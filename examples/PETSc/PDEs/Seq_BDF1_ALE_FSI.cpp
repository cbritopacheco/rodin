/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @example Time-dependent BDF1 ALE FSI with explicit solid-fluid coupling
 *
 * This example is intentionally small.  It couples:
 *   - a transient Oseen fluid solve on a moving ALE mesh, and
 *   - a compressible NeoHookean solid strip.
 *
 * The geometry consists of two matching triangular strip meshes:
 *
 *   solid:  [0, L] x [H_f, H_f + H_s]
 *   fluid:  [0, L] x [0,   H_f]
 *
 * The interface is the common line y = H_f.  Matching interface faces are kept
 * in an explicit map.  Interface interpolation is done by forwarding a point to
 * the matching face on the other mesh:
 *
 *   Point q(p);
 *   q.setPolytope(mappedFace);
 *
 * The time loop is a staggered, explicit FSI split:
 *
 *   1. Use the old solid displacement / velocity to update the fluid ALE mesh.
 *   2. Solve one BDF1 Oseen step on the moving fluid mesh.
 *   3. Transfer the current fluid pressure to the solid interface.
 *   4. Solve one BDF1 hyperelastic solid step.
 *
 * This is not a production FSI algorithm: no subiterations, no conservative
 * mortar transfer, and only a pressure load is applied to the solid.  The point
 * is to show the simplest Rodin mechanics for matching-element maps, ALE mesh
 * motion, Oseen flow, and hyperelasticity in one place.
 *
 * Suggested command:
 *
 * ./examples/PETSc/PDEs/PETSc_Seq_BDF1_ALE_FSI \
 *   -ksp_type preonly -pc_type lu \
 *   -pc_factor_shift_type nonzero -pc_factor_shift_amount 1e-10
 *
 * The option -fsi_pressure_scale controls the toy pressure load transferred to
 * the solid.  It defaults to a small nondimensional value so the explicit
 * partitioned step remains gentle on coarse meshes.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include <Rodin/Types.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Solid.h>
#include <Rodin/Solver.h>
#include <Rodin/Solver/NewtonSolver.h>
#include <Rodin/Solver/SparseLU.h>
#include <Rodin/Variational.h>

#include <Rodin/PETSc.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Math;
using namespace Rodin::Solver;
using namespace Rodin::Variational;

namespace
{
  using LocalMesh = Mesh<Context::Local>;

  struct FluidBoundary
  {
    static constexpr Attribute Inlet    = 1;
    static constexpr Attribute Outlet   = 2;
    static constexpr Attribute Wall     = 3;
    static constexpr Attribute FSI      = 4;
  };

  struct SolidBoundary
  {
    static constexpr Attribute FSI        = 11;
    static constexpr Attribute ClampLeft  = 12;
    static constexpr Attribute ClampRight = 13;
    static constexpr Attribute Free       = 14;
  };

  struct InterfaceSegment
  {
    Index fluid;
    Index solid;
    Real xmin;
    Real xmax;
  };

  struct InterfaceMap
  {
    std::unordered_map<Index, Index> fluidToSolid;
    std::unordered_map<Index, Index> solidToFluid;
    std::vector<InterfaceSegment> segments;
  };

  static Real pulse(const Real t)
  {
    const Real pi = Constants::pi();
    return 0.25 * std::sin(2.0 * pi * t / 0.4);
  }

  static Math::SpatialPoint centroid(const LocalMesh& mesh, const Polytope& polytope)
  {
    Math::SpatialPoint c(mesh.getSpaceDimension());
    c.setZero();

    const auto& vertices = polytope.getVertices();
    for (const auto& v : vertices)
      c += mesh.getVertexCoordinates(v);
    c /= static_cast<Real>(vertices.size());

    return c;
  }

  static std::pair<Real, Real> segmentXRange(const LocalMesh& mesh, const Polytope& segment)
  {
    const auto& vertices = segment.getVertices();
    const Real x0 = mesh.getVertexCoordinates(vertices[0]).x();
    const Real x1 = mesh.getVertexCoordinates(vertices[1]).x();
    return { std::min(x0, x1), std::max(x0, x1) };
  }

  static Real segmentReferenceCoordinate(const LocalMesh& mesh, const Polytope& segment, Real x)
  {
    const auto& vertices = segment.getVertices();
    const Real x0 = mesh.getVertexCoordinates(vertices[0]).x();
    const Real x1 = mesh.getVertexCoordinates(vertices[1]).x();

    if (std::abs(x1 - x0) < 1e-14)
      return 0.5;

    return std::clamp((x - x0) / (x1 - x0), 0.0, 1.0);
  }

  static void mapMeshToBox(LocalMesh& mesh, Real L, Real H, Real y0)
  {
    Real xmin =  std::numeric_limits<Real>::max();
    Real xmax = -std::numeric_limits<Real>::max();
    Real ymin =  std::numeric_limits<Real>::max();
    Real ymax = -std::numeric_limits<Real>::max();

    for (auto it = mesh.getVertex(); it; ++it)
    {
      const auto x = mesh.getVertexCoordinates(it->getIndex());
      xmin = std::min(xmin, x.x());
      xmax = std::max(xmax, x.x());
      ymin = std::min(ymin, x.y());
      ymax = std::max(ymax, x.y());
    }

    for (auto it = mesh.getVertex(); it; ++it)
    {
      auto x = mesh.getVertexCoordinates(it->getIndex());
      x(0) = L * (x.x() - xmin) / (xmax - xmin);
      x(1) = y0 + H * (x.y() - ymin) / (ymax - ymin);
      mesh.setVertexCoordinates(it->getIndex(), x);
    }
    mesh.flush();
  }

  static LocalMesh makeStripMesh(size_t nx, size_t ny, Real L, Real H, Real y0)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { nx + 1, ny + 1 });
    mapMeshToBox(mesh, L, H, y0);
    mesh.getConnectivity().compute(1, 2);
    return mesh;
  }

  static void labelFluidBoundaries(LocalMesh& mesh, Real L, Real H)
  {
    constexpr Real eps = 1e-10;
    for (auto it = mesh.getBoundary(); it; ++it)
    {
      const auto c = centroid(mesh, *it);
      if (std::abs(c.x()) < eps)
        mesh.setAttribute(it.key(), FluidBoundary::Inlet);
      else if (std::abs(c.x() - L) < eps)
        mesh.setAttribute(it.key(), FluidBoundary::Outlet);
      else if (std::abs(c.y()) < eps)
        mesh.setAttribute(it.key(), FluidBoundary::Wall);
      else if (std::abs(c.y() - H) < eps)
        mesh.setAttribute(it.key(), FluidBoundary::FSI);
    }
  }

  static void labelSolidBoundaries(LocalMesh& mesh, Real L, Real Hf, Real Hs)
  {
    constexpr Real eps = 1e-10;
    const Real yTop = Hf + Hs;

    for (auto it = mesh.getBoundary(); it; ++it)
    {
      const auto c = centroid(mesh, *it);
      if (std::abs(c.y() - Hf) < eps)
        mesh.setAttribute(it.key(), SolidBoundary::FSI);
      else if (std::abs(c.x()) < eps)
        mesh.setAttribute(it.key(), SolidBoundary::ClampLeft);
      else if (std::abs(c.x() - L) < eps)
        mesh.setAttribute(it.key(), SolidBoundary::ClampRight);
      else if (std::abs(c.y() - yTop) < eps)
        mesh.setAttribute(it.key(), SolidBoundary::Free);
    }
  }

  static InterfaceMap buildInterfaceMap(
      const LocalMesh& fluidMesh,
      const LocalMesh& solidMesh)
  {
    InterfaceMap map;
    std::map<long long, Index> solidByMidpoint;

    auto key = [](Real x) -> long long
    {
      return static_cast<long long>(std::llround(1e12 * x));
    };

    for (auto it = solidMesh.getBoundary(); it; ++it)
    {
      if (it->getAttribute() != SolidBoundary::FSI)
        continue;
      solidByMidpoint.emplace(key(centroid(solidMesh, *it).x()), it->getIndex());
    }

    for (auto it = fluidMesh.getBoundary(); it; ++it)
    {
      if (it->getAttribute() != FluidBoundary::FSI)
        continue;

      const auto c = centroid(fluidMesh, *it);
      const auto found = solidByMidpoint.find(key(c.x()));
      if (found == solidByMidpoint.end())
        throw std::runtime_error("Failed to match a fluid interface segment.");

      const Index fluidFace = it->getIndex();
      const Index solidFace = found->second;

      map.fluidToSolid.emplace(fluidFace, solidFace);
      map.solidToFluid.emplace(solidFace, fluidFace);

      const auto [xmin, xmax] = segmentXRange(fluidMesh, *it);
      map.segments.push_back({ fluidFace, solidFace, xmin, xmax });
    }

    std::sort(
        map.segments.begin(),
        map.segments.end(),
        [](const InterfaceSegment& a, const InterfaceSegment& b)
        {
          return a.xmin < b.xmin;
        });

    return map;
  }

  static Point forwardFluidPointToSolid(
      const Point& p,
      const LocalMesh& solidMesh,
      const InterfaceMap& map)
  {
    const auto found = map.fluidToSolid.find(p.getPolytope().getIndex());
    if (found == map.fluidToSolid.end())
      throw std::runtime_error("Fluid point is not on a mapped FSI face.");

    auto solidFace = solidMesh.getFace(found->second);
    Point q(p);
    q.setPolytope(*solidFace);
    return q;
  }

  static Point forwardSolidPointToFluid(
      const Point& p,
      const LocalMesh& fluidMesh,
      const InterfaceMap& map)
  {
    const auto found = map.solidToFluid.find(p.getPolytope().getIndex());
    if (found == map.solidToFluid.end())
      throw std::runtime_error("Solid point is not on a mapped FSI face.");

    auto fluidFace = fluidMesh.getFace(found->second);
    Point q(p);
    q.setPolytope(*fluidFace);
    return q;
  }

  static Point pointOnFluidInterface(
      const LocalMesh& fluidMesh,
      const InterfaceMap& map,
      const Real x)
  {
    constexpr Real eps = 1e-12;
    for (const auto& seg : map.segments)
    {
      if (x + eps < seg.xmin || x - eps > seg.xmax)
        continue;

      auto face = fluidMesh.getFace(seg.fluid);
      const Real r = segmentReferenceCoordinate(fluidMesh, *face, x);
      return Point(*face, Math::SpatialPoint{ r });
    }
    throw std::runtime_error("Point abscissa is outside the mapped FSI interface.");
  }

  template <class SolidGridFunction>
  static Math::SpatialVector<Real> evaluateSolidOnInterface(
      const SolidGridFunction& field,
      const LocalMesh& fluidMesh,
      const LocalMesh& solidMesh,
      const InterfaceMap& map,
      const Real x)
  {
    Point pf = pointOnFluidInterface(fluidMesh, map, x);
    Point ps = forwardFluidPointToSolid(pf, solidMesh, map);
    return field(ps);
  }

  template <class SolidGridFunction>
  static void updateFluidALEMesh(
      LocalMesh& fluidMesh,
      const std::vector<Math::SpatialPoint>& fluidReferenceVertices,
      const LocalMesh& solidMesh,
      const InterfaceMap& map,
      const SolidGridFunction& solidDisplacement,
      const Real Hf)
  {
    for (auto it = fluidMesh.getVertex(); it; ++it)
    {
      const Index vi = it->getIndex();
      const auto& X = fluidReferenceVertices[vi];

      auto x = X;
      const Real eta = std::clamp(X.y() / Hf, 0.0, 1.0);
      const auto dGamma =
        evaluateSolidOnInterface(solidDisplacement, fluidMesh, solidMesh, map, X.x());

      x += eta * dGamma;
      fluidMesh.setVertexCoordinates(vi, x);
    }
    fluidMesh.flush();
  }

  static bool isFinite(const Math::SpatialVector<Real>& x)
  {
    for (Index i = 0; i < x.size(); ++i)
    {
      if (!std::isfinite(x(i)))
        return false;
    }
    return true;
  }
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  PetscInt nxOpt = 24;
  PetscInt fluidNyOpt = 6;
  PetscInt solidNyOpt = 3;
  PetscInt ntOpt = 25;
  PetscReal finalTimeOpt = 0.25;
  PetscReal pressureScaleOpt = 1e-12;

  PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fsi_nx", &nxOpt, PETSC_NULLPTR);
  PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fsi_fluid_ny", &fluidNyOpt, PETSC_NULLPTR);
  PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fsi_solid_ny", &solidNyOpt, PETSC_NULLPTR);
  PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fsi_Nt", &ntOpt, PETSC_NULLPTR);
  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fsi_T", &finalTimeOpt, PETSC_NULLPTR);
  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fsi_pressure_scale", &pressureScaleOpt, PETSC_NULLPTR);

  const size_t nx = static_cast<size_t>(std::max<PetscInt>(2, nxOpt));
  const size_t fluidNy = static_cast<size_t>(std::max<PetscInt>(2, fluidNyOpt));
  const size_t solidNy = static_cast<size_t>(std::max<PetscInt>(1, solidNyOpt));
  const Index Nt = static_cast<Index>(std::max<PetscInt>(1, ntOpt));

  const Real L  = 4.0;
  const Real Hf = 1.0;
  const Real Hs = 0.25;
  const Real T  = static_cast<Real>(finalTimeOpt);
  const Real dt = T / static_cast<Real>(Nt);
  const Real pressureScale = static_cast<Real>(pressureScaleOpt);

  {
  LocalMesh fluidMesh = makeStripMesh(nx, fluidNy, L, Hf, 0.0);
  LocalMesh solidMesh = makeStripMesh(nx, solidNy, L, Hs, Hf);

  labelFluidBoundaries(fluidMesh, L, Hf);
  labelSolidBoundaries(solidMesh, L, Hf, Hs);

  const InterfaceMap interfaceMap = buildInterfaceMap(fluidMesh, solidMesh);

  std::vector<Math::SpatialPoint> fluidReferenceVertices(fluidMesh.getVertexCount());
  for (auto it = fluidMesh.getVertex(); it; ++it)
    fluidReferenceVertices[it->getIndex()] = fluidMesh.getVertexCoordinates(it->getIndex());

  const size_t dim = fluidMesh.getSpaceDimension();

  H1 uh(std::integral_constant<size_t, 2>{}, fluidMesh, dim);
  H1 ph(std::integral_constant<size_t, 1>{}, fluidMesh);
  P0g lh(fluidMesh);
  P1 solidVh(solidMesh, dim);

  PETSc::Variational::TrialFunction u(uh); u.setName("FluidVelocity");
  PETSc::Variational::TrialFunction p(ph); p.setName("FluidPressure");
  PETSc::Variational::TrialFunction l(lh); l.setName("FluidPressureMean");
  PETSc::Variational::TestFunction  v(uh);
  PETSc::Variational::TestFunction  q(ph);
  PETSc::Variational::TestFunction  m(lh);

  PETSc::Variational::GridFunction uOld(uh);
  PETSc::Variational::GridFunction pOld(ph);
  PETSc::Variational::GridFunction meshVelocity(uh);
  PETSc::Variational::GridFunction transportVelocity(uh);
  PETSc::Variational::GridFunction inletVelocity(uh);

  uOld.setName("FluidVelocityOld");
  pOld.setName("FluidPressureOld");
  meshVelocity.setName("ALEMeshVelocity");
  transportVelocity.setName("ALETransportVelocity");

  uOld = Zero(dim);
  pOld = 0.0;
  meshVelocity = Zero(dim);
  transportVelocity = Zero(dim);
  inletVelocity = Zero(dim);

  GridFunction solidDisplacement(solidVh);
  GridFunction solidDisplacementOld(solidVh);
  GridFunction solidVelocity(solidVh);
  GridFunction solidVelocityOld(solidVh);
  GridFunction solidFluidTraction(solidVh);

  solidDisplacement.setName("SolidDisplacement");
  solidVelocity.setName("SolidVelocity");
  solidFluidTraction.setName("FluidPressureTraction");

  const auto zero = Zero(dim);
  solidDisplacement = zero;
  solidDisplacementOld = zero;
  solidVelocity = zero;
  solidVelocityOld = zero;
  solidFluidTraction = zero;

  IO::XDMF xdmf("BDF1_ALE_FSI");
  auto fluidGrid = xdmf.grid("Fluid");
  fluidGrid.setMesh(fluidMesh, IO::XDMF::MeshPolicy::Transient);
  fluidGrid.add("velocity", u.getSolution());
  fluidGrid.add("pressure", p.getSolution());
  fluidGrid.add("mesh_velocity", meshVelocity);

  auto solidGrid = xdmf.grid("Solid");
  solidGrid.setMesh(solidMesh);
  solidGrid.add(solidDisplacement);
  solidGrid.add(solidVelocity);
  solidGrid.add(solidFluidTraction);
  xdmf.write(0.0);

  const Real rhoF = 1.0;
  const Real muF  = 0.02;

  const Real solidE = 8.0e2;
  const Real solidNu = 0.3;
  const Real solidLambda =
    solidE * solidNu / ((1.0 + solidNu) * (1.0 - 2.0 * solidNu));
  const Real solidMu = solidE / (2.0 * (1.0 + solidNu));
  const Real rhoS = 5.0;
  const Real dampingS = 0.25;
  const Real solidMassCoeff = rhoS / (dt * dt) + dampingS / dt;

  Solid::NeoHookean law(solidLambda, solidMu);
  auto nFluid = BoundaryNormal(fluidMesh);

  Real t = 0.0;
  for (Index step = 1; step <= Nt; ++step)
  {
    t += dt;

    updateFluidALEMesh(
        fluidMesh,
        fluidReferenceVertices,
        solidMesh,
        interfaceMap,
        solidDisplacementOld,
        Hf);

    meshVelocity = VectorFunction(dim, [&](const Point& x)
    {
      const Real eta = std::clamp(x.y() / Hf, 0.0, 1.0);
      return eta * evaluateSolidOnInterface(
          solidVelocityOld, fluidMesh, solidMesh, interfaceMap, x.x());
    });

    inletVelocity = VectorFunction(dim, [&](const Point& x)
    {
      const Real y = std::clamp(x.y(), 0.0, Hf);
      const Real profile = 4.0 * y * (Hf - y) / (Hf * Hf);
      return Math::Vector<Real>{{ pulse(t) * profile, 0.0 }};
    });

    auto interfaceVelocity = VectorFunction(dim, [&](const Point& x)
    {
      const Point xs = forwardFluidPointToSolid(x, solidMesh, interfaceMap);
      return solidVelocityOld(xs);
    });

    transportVelocity = uOld - meshVelocity;

    const auto convU = Mult(Jacobian(u), transportVelocity);
    const auto divTransport = Div(transportVelocity);
    const auto beta = Max(-Dot(transportVelocity, nFluid), 0.0);

    Problem flow(u, p, l, v, q, m);
    flow =
          (rhoF / dt) * Integral(u, v)
        - (rhoF / dt) * Integral(uOld, v)
        + rhoF * Integral(Dot(convU, v))
        + 0.5 * rhoF * Integral(divTransport * Dot(u, v))
        + muF * Integral(Jacobian(u), Jacobian(v))
        - Integral(p, Div(v))
        + Integral(Div(u), q)
        + Integral(l, q)
        + Integral(p, m)
        + BoundaryIntegral(0.5 * rhoF * beta * Dot(u, v)).over(FluidBoundary::Outlet)
        + 1e-10 * Integral(p, q)
        + 1e-10 * Integral(l, m)
        + DirichletBC(u, inletVelocity).on(FluidBoundary::Inlet)
        + DirichletBC(u, Zero(dim)).on(FluidBoundary::Wall)
        + DirichletBC(u, interfaceVelocity).on(FluidBoundary::FSI);

    Alert::Info() << "Fluid BDF1 ALE/Oseen step " << step << " / " << Nt << Alert::Raise;
    flow.assemble().setFieldSplits();
    Solver::KSP(flow).solve();

    uOld.setData(u.getSolution().getData());
    pOld.setData(p.getSolution().getData());

    TrialFunction du(solidVh);
    TestFunction  w(solidVh);

    Solid::MaterialTangent tangent(law, du, w);
    tangent.setDisplacement(solidDisplacement);

    Solid::InternalForce internal(law, w);
    internal.setDisplacement(solidDisplacement);

    auto fluidPressureLoad = VectorFunction(dim, [&](const Point& xs)
    {
      Math::SpatialVector<Real> traction(dim);
      traction.setZero();
      if (pressureScale == 0.0)
        return traction;

      const Point xf = forwardSolidPointToFluid(xs, fluidMesh, interfaceMap);
      const Real pressure = pOld(xf);
      const auto normal = nFluid(xf);
      if (!std::isfinite(pressure) || !isFinite(normal))
        return traction;

      traction = pressureScale * pressure * normal;
      return traction;
    });

    solidFluidTraction = zero;
    solidFluidTraction.project(
        Region::Boundary,
        fluidPressureLoad,
        SolidBoundary::FSI);

    Problem solid(du, w);
    solid =
          tangent
        + solidMassCoeff * Integral(du, w)
        + internal
        + solidMassCoeff * Integral(solidDisplacement, w)
        - solidMassCoeff * Integral(solidDisplacementOld, w)
        - (rhoS / dt) * Integral(solidVelocityOld, w)
        - BoundaryIntegral(solidFluidTraction, w).over(SolidBoundary::FSI)
        + DirichletBC(du, Zero(dim)).on(SolidBoundary::ClampLeft)
        + DirichletBC(du, Zero(dim)).on(SolidBoundary::ClampRight);

    SparseLU solidLinearSolver(solid);
    NewtonSolver solidSolver(solidLinearSolver);
    solidSolver.setMaxIterations(25)
               .setDampingFactor(1.0)
               .setAbsoluteTolerance(1e-10)
               .setRelativeTolerance(1e-8);

    Alert::Info() << "Solid BDF1 hyperelastic step " << step << " / " << Nt << Alert::Raise;
    solidSolver.solve(solidDisplacement);

    solidVelocity = (1.0 / dt) * (solidDisplacement - solidDisplacementOld);
    solidDisplacementOld = solidDisplacement;
    solidVelocityOld = solidVelocity;

    xdmf.write(t).flush();
  }

  xdmf.close();
  }

  PetscFinalize();

  return 0;
}
