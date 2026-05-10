/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @example Time-dependent BDF1 ALE FSI with explicit solid-fluid coupling
 *
 * Geometry convention:
 *
 *   solidMesh:
 *     reference solid mesh, used for hyperelastic assembly.
 *
 *   fluidMesh:
 *     current ALE fluid mesh, used for fluid assembly.
 *
 *   solidOutputMesh:
 *     moved solid mesh, used only for XDMF output.
 *
 * The fluid ALE displacement d_f is computed by harmonic extension:
 *
 *   -Δ d_f = 0       in Ω_f
 *      d_f = d_s     on Γ_FSI
 *      d_f = 0       on inlet, outlet, and wall
 *
 * Then the fluid mesh is moved by:
 *
 *   x_f = X_f + d_f.
 *
 * The solid output mesh is moved by:
 *
 *   x_s = X_s + d_s.
 *
 * Important visualization note:
 *
 *   The solid mesh written to XDMF is already geometrically moved. Do not
 *   additionally apply ParaView "Warp By Vector" to SolidDisplacement unless
 *   you intentionally want to see X_s + 2 d_s.
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

  static Real pulseTrain(
      const Real t,
      const Real amplitude,
      const Real duration,
      const Real gap,
      const Index count)
  {
    if (t < 0.0 || duration <= 0.0 || count == 0)
      return 0.0;

    const Real period = duration + std::max<Real>(gap, 0.0);
    const Index pulse = static_cast<Index>(std::floor(t / period));

    if (count > 0 && pulse >= count)
      return 0.0;

    const Real localTime = t - static_cast<Real>(pulse) * period;
    if (localTime > duration)
      return 0.0;

    const Real pi = Constants::pi();

    return amplitude * std::sin(pi * localTime / duration);
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
      const LocalMesh& fluidReferenceMesh,
      const LocalMesh& solidReferenceMesh)
  {
    InterfaceMap map;
    std::map<long long, Index> solidByMidpoint;

    auto key = [](Real x) -> long long
    {
      return static_cast<long long>(std::llround(1e12 * x));
    };

    for (auto it = solidReferenceMesh.getBoundary(); it; ++it)
    {
      if (it->getAttribute() != SolidBoundary::FSI)
        continue;

      solidByMidpoint.emplace(
          key(centroid(solidReferenceMesh, *it).x()),
          it->getIndex());
    }

    for (auto it = fluidReferenceMesh.getBoundary(); it; ++it)
    {
      if (it->getAttribute() != FluidBoundary::FSI)
        continue;

      const auto c = centroid(fluidReferenceMesh, *it);
      const auto found = solidByMidpoint.find(key(c.x()));

      if (found == solidByMidpoint.end())
        throw std::runtime_error("Failed to match a fluid interface segment.");

      const Index fluidFace = it->getIndex();
      const Index solidFace = found->second;

      map.fluidToSolid.emplace(fluidFace, solidFace);
      map.solidToFluid.emplace(solidFace, fluidFace);

      const auto [xmin, xmax] = segmentXRange(fluidReferenceMesh, *it);
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
      const LocalMesh& solidReferenceMesh,
      const InterfaceMap& map)
  {
    const auto found = map.fluidToSolid.find(p.getPolytope().getIndex());

    if (found == map.fluidToSolid.end())
      throw std::runtime_error("Fluid point is not on a mapped FSI face.");

    auto solidFace = solidReferenceMesh.getFace(found->second);

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

  static Point pointOnVertex(const LocalMesh& mesh, Index vi)
  {
    auto v = mesh.getVertex(vi);

    Math::SpatialPoint rc(1);
    rc.setZero();

    return Point(*v, rc);
  }

  static void saveReferenceVertices(
      const LocalMesh& mesh,
      std::vector<Math::SpatialPoint>& vertices)
  {
    vertices.resize(mesh.getVertexCount());

    for (auto it = mesh.getVertex(); it; ++it)
      vertices[it->getIndex()] = mesh.getVertexCoordinates(it->getIndex());
  }

  static void restoreMesh(
      LocalMesh& mesh,
      const std::vector<Math::SpatialPoint>& referenceVertices)
  {
    for (auto it = mesh.getVertex(); it; ++it)
    {
      const Index vi = it->getIndex();
      mesh.setVertexCoordinates(vi, referenceVertices[vi]);
    }

    mesh.flush();
  }

  template <class GridFunctionType>
  static void moveMeshWithDisplacement(
      LocalMesh& mesh,
      const LocalMesh& referenceMesh,
      const std::vector<Math::SpatialPoint>& referenceVertices,
      const GridFunctionType& displacement)
  {
    for (auto it = mesh.getVertex(); it; ++it)
    {
      const Index vi = it->getIndex();

      auto x = referenceVertices[vi];

      const Point p = pointOnVertex(referenceMesh, vi);
      x += displacement(p);

      mesh.setVertexCoordinates(vi, x);
    }

    mesh.flush();
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

  static bool isFinite(const Math::SpatialMatrix<Real>& A)
  {
    for (int i = 0; i < A.rows(); ++i)
    {
      for (int j = 0; j < A.cols(); ++j)
      {
        if (!std::isfinite(A(i, j)))
          return false;
      }
    }

    return true;
  }
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  PetscInt nxOpt = 48;
  PetscInt fluidNyOpt = 12;
  PetscInt solidNyOpt = 12;
  PetscInt ntOpt = 100;

  PetscReal finalTimeOpt = 1.0;

  PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fsi_nx", &nxOpt, PETSC_NULLPTR);
  PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fsi_fluid_ny", &fluidNyOpt, PETSC_NULLPTR);
  PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fsi_solid_ny", &solidNyOpt, PETSC_NULLPTR);
  PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fsi_Nt", &ntOpt, PETSC_NULLPTR);
  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fsi_T", &finalTimeOpt, PETSC_NULLPTR);

  PetscReal pressureScaleOpt = 0.05;
  PetscReal viscousScaleOpt = 1.0;
  PetscReal tractionScaleOpt = 1.0;
  PetscReal backflowStabilizationScaleOpt = 1.0;
  PetscReal inletAmplitudeOpt = 1.25;
  PetscReal pulseDurationOpt = 0.18;
  PetscReal pulseGapOpt = 0.35;
  PetscInt pulseCountOpt = 1;
  PetscReal solidEOpt = 2.0e2;
  PetscReal solidNuOpt = 0.3;
  PetscReal solidDensityOpt = 5.0e-1;
  PetscReal solidRayleighAlphaOpt = 5.0e-2;
  PetscReal solidRayleighBetaOpt = 2.5e-4;
  PetscReal newtonDampingOpt = 1.0;
  PetscBool checkNoSlipOpt = PETSC_TRUE;

  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
      "-fsi_pressure_scale", &pressureScaleOpt, PETSC_NULLPTR);

  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
      "-fsi_viscous_scale", &viscousScaleOpt, PETSC_NULLPTR);

  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
      "-fsi_traction_scale", &tractionScaleOpt, PETSC_NULLPTR);

  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
      "-fsi_backflow_stabilization_scale",
      &backflowStabilizationScaleOpt,
      PETSC_NULLPTR);

  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
      "-fsi_inlet_amplitude", &inletAmplitudeOpt, PETSC_NULLPTR);

  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
      "-fsi_pulse_duration", &pulseDurationOpt, PETSC_NULLPTR);

  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
      "-fsi_pulse_gap", &pulseGapOpt, PETSC_NULLPTR);

  PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR,
      "-fsi_pulse_count", &pulseCountOpt, PETSC_NULLPTR);

  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
      "-fsi_solid_E", &solidEOpt, PETSC_NULLPTR);

  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
      "-fsi_solid_nu", &solidNuOpt, PETSC_NULLPTR);

  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
      "-fsi_solid_density", &solidDensityOpt, PETSC_NULLPTR);

  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
      "-fsi_solid_rayleigh_alpha", &solidRayleighAlphaOpt, PETSC_NULLPTR);

  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
      "-fsi_solid_damping", &solidRayleighAlphaOpt, PETSC_NULLPTR);

  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
      "-fsi_solid_rayleigh_beta", &solidRayleighBetaOpt, PETSC_NULLPTR);

  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
      "-fsi_solid_newton_damping", &newtonDampingOpt, PETSC_NULLPTR);

  PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR,
      "-fsi_check_no_slip", &checkNoSlipOpt, PETSC_NULLPTR);

  const size_t nx = static_cast<size_t>(std::max<PetscInt>(2, nxOpt));
  const size_t fluidNy = static_cast<size_t>(std::max<PetscInt>(2, fluidNyOpt));
  const size_t solidNy = static_cast<size_t>(std::max<PetscInt>(1, solidNyOpt));
  const Index Nt = static_cast<Index>(std::max<PetscInt>(1, ntOpt));

  const Real L  = 4.0;
  const Real Hf = 1.0;
  const Real Hs = 0.25;

  const Real T  = static_cast<Real>(finalTimeOpt);
  const Real dt = T / static_cast<Real>(Nt);

  {
    /*
     * Meshes.
     *
     * fluidMesh is mutable/current.
     * solidMesh remains reference.
     * solidOutputMesh is only for XDMF output.
     */
    LocalMesh fluidMesh = makeStripMesh(nx, fluidNy, L, Hf, 0.0);
    LocalMesh solidMesh = makeStripMesh(nx, solidNy, L, Hs, Hf);
    LocalMesh solidOutputMesh = solidMesh;

    labelFluidBoundaries(fluidMesh, L, Hf);
    labelSolidBoundaries(solidMesh, L, Hf, Hs);
    labelSolidBoundaries(solidOutputMesh, L, Hf, Hs);

    std::vector<Math::SpatialPoint> fluidReferenceVertices;
    std::vector<Math::SpatialPoint> solidReferenceVertices;

    saveReferenceVertices(fluidMesh, fluidReferenceVertices);
    saveReferenceVertices(solidMesh, solidReferenceVertices);

    const InterfaceMap interfaceMap = buildInterfaceMap(fluidMesh, solidMesh);

    const size_t dim = fluidMesh.getSpaceDimension();

    /*
     * Function spaces.
     */
    H1 uh(std::integral_constant<size_t, 2>{}, fluidMesh, dim);
    H1 ph(std::integral_constant<size_t, 1>{}, fluidMesh);
    P0g lh(fluidMesh);

    P1 solidVh(solidMesh, dim);
    P1 solidOutputVh(solidOutputMesh, dim);

    /*
     * Fluid unknowns.
     */
    PETSc::Variational::TrialFunction u(uh);
    PETSc::Variational::TrialFunction p(ph);
    PETSc::Variational::TrialFunction l(lh);

    u.setName("FluidVelocity");
    p.setName("FluidPressure");
    l.setName("FluidPressureMean");

    PETSc::Variational::TestFunction v(uh);
    PETSc::Variational::TestFunction q(ph);
    PETSc::Variational::TestFunction m(lh);

    PETSc::Variational::GridFunction uOld(uh);
    PETSc::Variational::GridFunction pOld(ph);
    PETSc::Variational::GridFunction meshVelocity(uh);
    PETSc::Variational::GridFunction transportVelocity(uh);
    PETSc::Variational::GridFunction inletVelocity(uh);

    PETSc::Variational::GridFunction fluidDisplacement(uh);
    PETSc::Variational::GridFunction fluidDisplacementOld(uh);

    uOld.setName("FluidVelocityOld");
    pOld.setName("FluidPressureOld");
    meshVelocity.setName("ALEMeshVelocity");
    transportVelocity.setName("ALETransportVelocity");
    fluidDisplacement.setName("FluidALEDisplacement");

    uOld = Zero(dim);
    pOld = 0.0;
    meshVelocity = Zero(dim);
    transportVelocity = Zero(dim);
    inletVelocity = Zero(dim);
    fluidDisplacement = Zero(dim);
    fluidDisplacementOld = Zero(dim);

    /*
     * Solid assembly fields on reference solidMesh.
     */
    GridFunction solidDisplacement(solidVh);
    GridFunction solidDisplacementOld(solidVh);
    GridFunction solidVelocity(solidVh);
    GridFunction solidVelocityOld(solidVh);
    GridFunction solidFluidTraction(solidVh);

    solidDisplacement.setName("SolidDisplacement");
    solidVelocity.setName("SolidVelocity");
    solidFluidTraction.setName("FluidCauchyTraction");

    const auto zero = Zero(dim);

    solidDisplacement = zero;
    solidDisplacementOld = zero;
    solidVelocity = zero;
    solidVelocityOld = zero;
    solidFluidTraction = zero;

    /*
     * Solid output fields on solidOutputMesh.
     *
     * These are necessary because XDMF requires the attribute mesh to match
     * the grid mesh.
     */
    GridFunction solidOutputDisplacement(solidOutputVh);
    GridFunction solidOutputVelocity(solidOutputVh);
    GridFunction solidOutputFluidTraction(solidOutputVh);

    solidOutputDisplacement.setName("SolidDisplacement");
    solidOutputVelocity.setName("SolidVelocity");
    solidOutputFluidTraction.setName("FluidCauchyTraction");

    solidOutputDisplacement = zero;
    solidOutputVelocity = zero;
    solidOutputFluidTraction = zero;

    /*
     * XDMF output.
     */
    IO::XDMF xdmf("out/BDF1_ALE_FSI_HarmonicALE");

    auto fluidGrid = xdmf.grid("Fluid");
    fluidGrid.setMesh(fluidMesh, IO::XDMF::MeshPolicy::Transient);
    fluidGrid.add("velocity", u.getSolution());
    fluidGrid.add("pressure", p.getSolution());
    fluidGrid.add("mesh_velocity", meshVelocity);
    fluidGrid.add("ale_displacement", fluidDisplacement);

    auto solidGrid = xdmf.grid("Solid");
    solidGrid.setMesh(solidOutputMesh, IO::XDMF::MeshPolicy::Transient);
    solidGrid.add(solidOutputDisplacement);
    solidGrid.add(solidOutputVelocity);
    solidGrid.add(solidOutputFluidTraction);

    solidOutputDisplacement.setData(solidDisplacement.getData());
    solidOutputVelocity.setData(solidVelocity.getData());
    solidOutputFluidTraction.setData(solidFluidTraction.getData());

    xdmf.write(0.0);

    /*
     * Parameters.
     */
    const Real rhoF = 1.0;
    const Real muF  = 0.02;

    const Real solidE = static_cast<Real>(solidEOpt);
    const Real solidNu = static_cast<Real>(solidNuOpt);

    const Real solidLambda =
      solidE * solidNu / ((1.0 + solidNu) * (1.0 - 2.0 * solidNu));

    const Real solidMu =
      solidE / (2.0 * (1.0 + solidNu));

    const Real rhoS = static_cast<Real>(solidDensityOpt);
    const Real solidRayleighAlpha =
      static_cast<Real>(solidRayleighAlphaOpt);

    const Real solidRayleighBeta =
      static_cast<Real>(solidRayleighBetaOpt);

    const Real solidMassCoeff =
      rhoS / (dt * dt) + solidRayleighAlpha * rhoS / dt;

    Solid::NeoHookean law(solidLambda, solidMu);

    auto nFluid = BoundaryNormal(fluidMesh);

    /*
     * Harmonic ALE solve.
     *
     * This restores the fluid mesh to the reference coordinates, solves the
     * harmonic extension, stores it in fluidDisplacement, then moves fluidMesh.
     */
    auto solveHarmonicALE =
      [&](const auto& interfaceSolidDisplacement)
      {
        restoreMesh(fluidMesh, fluidReferenceVertices);

        PETSc::Variational::TrialFunction d(uh);
        PETSc::Variational::TestFunction  z(uh);

        d.setName("ALEDisplacementTrial");

        auto interfaceDisplacement = VectorFunction(dim, [&](const Point& xf)
        {
          const Point xs =
            forwardFluidPointToSolid(xf, solidMesh, interfaceMap);

          Math::SpatialVector<Real> value(dim);
          value.setZero();

          const auto ds = interfaceSolidDisplacement(xs);
          for (Index i = 0; i < static_cast<Index>(dim); ++i)
            value(i) = ds(i);

          return value;
        });

        Problem ale(d, z);

        ale =
            Integral(Jacobian(d), Jacobian(z))
          + DirichletBC(d, Zero(dim)).on(FluidBoundary::Inlet)
          + DirichletBC(d, Zero(dim)).on(FluidBoundary::Wall)
          + DirichletBC(d, interfaceDisplacement).on(FluidBoundary::FSI);

        Alert::Info()
          << "Solving harmonic ALE extension."
          << Alert::Raise;

        ale.assemble();
        Solver::KSP(ale).solve();

        fluidDisplacement.setData(d.getSolution().getData());

        moveMeshWithDisplacement(
            fluidMesh,
            fluidMesh,
            fluidReferenceVertices,
            fluidDisplacement);
      };

    Real t = 0.0;

    for (Index step = 1; step <= Nt; ++step)
    {
      t += dt;

      /*
       * Explicit ALE update for the fluid step:
       *
       * The fluid is solved on the mesh induced by the old solid displacement.
       */
      fluidDisplacementOld.setData(fluidDisplacement.getData());
      solveHarmonicALE(solidDisplacementOld);

      meshVelocity =
        (1.0 / dt) * (fluidDisplacement - fluidDisplacementOld);

      inletVelocity = VectorFunction(dim, [&](const Point& x)
      {
        const Real y = std::clamp(x.y(), 0.0, Hf);
        const Real profile = 4.0 * y * (Hf - y) / (Hf * Hf);

        return Math::Vector<Real>{{
            pulseTrain(
              t,
              static_cast<Real>(inletAmplitudeOpt),
              static_cast<Real>(pulseDurationOpt),
              static_cast<Real>(pulseGapOpt),
              static_cast<Index>(pulseCountOpt)) * profile,
            0.0 }};
      });

      auto interfaceVelocity = VectorFunction(dim, [&](const Point& xf)
      {
        const Point xs =
          forwardFluidPointToSolid(xf, solidMesh, interfaceMap);

        Math::SpatialVector<Real> value(dim);
        value.setZero();

        const auto vs = solidVelocityOld(xs);
        for (Index i = 0; i < static_cast<Index>(dim); ++i)
          value(i) = vs(i);

        return value;
      });

      transportVelocity = uOld - meshVelocity;

      const auto convU = Mult(Jacobian(u), transportVelocity);
      const auto divTransport = Div(transportVelocity);
      const auto beta = Max(-Dot(transportVelocity, nFluid), 0.0);
      const Real backflowStabilizationScale =
        static_cast<Real>(backflowStabilizationScaleOpt);

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
          + backflowStabilizationScale *
            BoundaryIntegral(
                0.5 * rhoF * beta * Dot(u, v)).over(FluidBoundary::Outlet)
          + 1e-10 * Integral(p, q)
          + 1e-10 * Integral(l, m)
          + DirichletBC(u, inletVelocity).on(FluidBoundary::Inlet)
          + DirichletBC(u, Zero(dim)).on(FluidBoundary::Wall)
          + DirichletBC(u, interfaceVelocity).on(FluidBoundary::FSI);

      Alert::Info()
        << "Fluid BDF1 ALE/Oseen step " << step << " / " << Nt
        << Alert::Raise;

      flow.assemble().setFieldSplits();
      Solver::KSP(flow).solve();

      if (checkNoSlipOpt)
      {
        auto fsiNoSlip =
          DirichletBC(u, interfaceVelocity).on(FluidBoundary::FSI);
        fsiNoSlip.assemble();

        const auto& fluidVelocityData = u.getSolution().getData();
        Real maxNoSlipResidual = 0.0;
        Real maxFluidFSIVelocity = 0.0;
        Real maxWallFSIVelocity = 0.0;

        for (const auto& [local, target] : fsiNoSlip.getDOFs())
        {
          const PetscInt idx = static_cast<PetscInt>(local);
          PetscScalar actual = 0.0;

          const PetscErrorCode ierr =
            VecGetValues(fluidVelocityData, 1, &idx, &actual);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;

          maxNoSlipResidual =
            std::max<Real>(
                maxNoSlipResidual,
                std::abs(static_cast<Real>(actual - target)));

          maxFluidFSIVelocity =
            std::max<Real>(
                maxFluidFSIVelocity,
                std::abs(static_cast<Real>(actual)));

          maxWallFSIVelocity =
            std::max<Real>(
                maxWallFSIVelocity,
                std::abs(static_cast<Real>(target)));
        }

        Alert::Info()
          << "  fluid FSI no-slip DOF max |u - u_FSI|: "
          << maxNoSlipResidual
          << " (max |u|=" << maxFluidFSIVelocity
          << ", max |u_FSI|=" << maxWallFSIVelocity << ")"
          << Alert::Raise;
      }

      uOld.setData(u.getSolution().getData());
      pOld.setData(p.getSolution().getData());

      /*
       * Solid tangent and internal force.
       */
      TrialFunction du(solidVh);
      TestFunction  w(solidVh);

      Solid::MaterialTangent tangent(law, du, w);
      tangent.setDisplacement(solidDisplacement);

      Solid::InternalForce internal(law, w);
      internal.setDisplacement(solidDisplacement);

      Solid::InternalForce internalOld(law, w);
      internalOld.setDisplacement(solidDisplacementOld);

      /*
       * Transfer fluid traction to the solid interface.
       */
      static bool tractionNaNReported = false;

      const Real pressureScale = static_cast<Real>(pressureScaleOpt);
      const Real viscousScale = static_cast<Real>(viscousScaleOpt);
      const Real tractionScale = static_cast<Real>(tractionScaleOpt);

      auto fluidTractionLoad = VectorFunction(dim, [&](const Point& xs)
      {
        Math::SpatialVector<Real> traction(dim);
        traction.setZero();

        const Point xf =
          forwardSolidPointToFluid(xs, fluidMesh, interfaceMap);

        const Real pressure = pOld(xf);
        const auto J = Jacobian(uOld)(xf);
        const auto normal = nFluid(xf);

        const bool pOk = std::isfinite(pressure);
        const bool nOk = isFinite(normal);
        const bool JOk =
          (J.rows() == static_cast<int>(dim)) &&
          (J.cols() == static_cast<int>(dim)) &&
          isFinite(J);

        if (!pOk || !nOk || !JOk)
        {
          if (!tractionNaNReported)
          {
            tractionNaNReported = true;

            const auto& pc = xs.getPhysicalCoordinates();

            Alert::Warning()
              << "fluidTractionLoad: non-finite at xs=("
              << pc.transpose() << ")"
              << " pressure=" << pressure
              << " normalOk=" << nOk
              << " J(" << J.rows() << "x" << J.cols() << ")Ok=" << JOk
              << Alert::Raise;
          }

          return traction;
        }

        const auto I = Math::SpatialMatrix<Real>::Identity(dim, dim);

        const auto sigmaF =
            -pressureScale * pressure * I
            + viscousScale * muF * (J + J.transpose());

        /*
         * nFluid is the outward normal of the fluid.
         * The traction exerted by the fluid on the solid is -sigma_f n_f.
         */
        traction = -tractionScale * sigmaF * normal;

        return traction;
      });

      solidFluidTraction = zero;

      solidFluidTraction.project(
          Region::Boundary,
          fluidTractionLoad,
          SolidBoundary::FSI);

      /*
       * Solid problem.
       *
       * Important:
       *   solidMesh is the reference mesh. It is not moved.
       */
      Problem solid(du, w);

      solid =
          tangent
        + solidMassCoeff * Integral(du, w)
        + (solidRayleighBeta / dt) * tangent

        + internal
        + solidMassCoeff * Integral(solidDisplacement, w)
        - solidMassCoeff * Integral(solidDisplacementOld, w)
        - (rhoS / dt) * Integral(solidVelocityOld, w)
        + (solidRayleighBeta / dt) * internal
        - (solidRayleighBeta / dt) * internalOld
        - BoundaryIntegral(solidFluidTraction, w).over(SolidBoundary::FSI)

        + DirichletBC(du, Zero(dim)).on(SolidBoundary::ClampLeft)
        + DirichletBC(du, Zero(dim)).on(SolidBoundary::ClampRight);

      Alert::Info()
        << "Solid BDF1 hyperelastic step " << step << " / " << Nt
        << Alert::Raise;

      Alert::Info()
        << "  solidDisplacement norm before Newton: "
        << solidDisplacement.getData().norm()
        << Alert::Raise;

      Alert::Info()
        << "  solidFluidTraction norm: "
        << solidFluidTraction.getData().norm()
        << Alert::Raise;

      const auto solidDisplacementInitial =
        solidDisplacement.getData();

      const Real tightSolidTolerance = 1e-10;
      const Real basinSolidTolerance = 1e-8;

      auto solveSolidNewton =
        [&](const char* label,
            Real damping,
            Real absoluteTolerance,
            Real relativeTolerance)
        {
          SparseLU solidLinearSolver(solid);
          NewtonSolver solidSolver(solidLinearSolver);

          solidSolver.setMaxIterations(75)
                     .setDampingFactor(damping)
                     .setAbsoluteTolerance(absoluteTolerance)
                     .setRelativeTolerance(relativeTolerance);

          Alert::Info()
            << "  " << label
            << " alpha=" << damping
            << " atol=" << absoluteTolerance
            << Alert::Raise;

          solidSolver.setMonitor(
              [](const auto& rep)
              {
                Alert::Info()
                  << "    it=" << rep.iterations
                  << " r=" << rep.final_residual
                  << " dx=" << rep.final_step_norm
                  << Alert::Raise;
              });

          solidSolver.solve(solidDisplacement);

          return solidSolver.converged();
        };

      solidDisplacement.setData(solidDisplacementInitial);

      bool solidNewtonConverged =
        solveSolidNewton(
            "full Newton attempt",
            1.0,
            tightSolidTolerance,
            1e-8);

      Real solidNewtonDamping =
        std::min<Real>(newtonDampingOpt, 0.5);

      for (Index attempt = 0; !solidNewtonConverged && attempt < 8; ++attempt)
      {
        solidDisplacement.setData(solidDisplacementInitial);

        Alert::Warning()
          << "  full Newton did not converge; seeking basin with alpha="
          << solidNewtonDamping
          << Alert::Raise;

        const bool reachedBasin =
          solveSolidNewton(
              "damped basin attempt",
              solidNewtonDamping,
              basinSolidTolerance,
              0.0);

        if (reachedBasin)
        {
          solidNewtonConverged =
            solveSolidNewton(
                "full Newton final",
                1.0,
                tightSolidTolerance,
                1e-8);

          if (solidNewtonConverged)
            break;

          Alert::Warning()
            << "  full Newton failed from damped iterate; reducing alpha."
            << Alert::Raise;
        }

        solidNewtonDamping *= 0.5;
      }

      if (!solidNewtonConverged)
      {
        throw std::runtime_error(
            "Solid Newton failed to converge after basin damping retries.");
      }

      /*
       * Commit solid state.
       */
      solidVelocity =
        (1.0 / dt) * (solidDisplacement - solidDisplacementOld);

      solidDisplacementOld = solidDisplacement;
      solidVelocityOld = solidVelocity;

      /*
       * Output geometry.
       *
       * Recompute harmonic ALE with the current solid displacement so that the
       * written fluid and solid meshes are geometrically compatible.
       *
       * Numerically, the fluid solve above was performed on the explicit old
       * ALE mesh. This recomputation is for visualization consistency.
       */
      fluidDisplacementOld.setData(fluidDisplacement.getData());
      solveHarmonicALE(solidDisplacement);

      meshVelocity =
        (1.0 / dt) * (fluidDisplacement - fluidDisplacementOld);

      restoreMesh(solidOutputMesh, solidReferenceVertices);

      moveMeshWithDisplacement(
          solidOutputMesh,
          solidMesh,
          solidReferenceVertices,
          solidDisplacement);

      /*
       * Copy solid assembly fields to output fields.
       *
       * The data layouts should match because solidMesh and solidOutputMesh
       * have the same topology and vertex ordering.
       */
      solidOutputDisplacement.setData(solidDisplacement.getData());
      solidOutputVelocity.setData(solidVelocity.getData());
      solidOutputFluidTraction.setData(solidFluidTraction.getData());

      xdmf.write(t).flush();
    }

    xdmf.close();
  }

  PetscFinalize();

  return 0;
}
