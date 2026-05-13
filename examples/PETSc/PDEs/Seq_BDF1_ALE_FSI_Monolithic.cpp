/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @example Monolithic BDF1 ALE FSI on a single conforming mesh
 *
 * This example is the monolithic counterpart of the staggered
 * Seq_BDF1_ALE_FSI.cpp example.
 *
 * The model is intentionally small.  The computational mesh is a single
 * conforming strip split into a fluid subdomain and a solid subdomain:
 *
 *   fluid:  [0, L] x [0,   Hf]
 *   solid:  [0, L] x [Hf, Hf + Hs]
 *
 * The interface is the internal facet set y = Hf.  The mesh carries cell
 * attributes Fluid/Solid and a facet attribute FSI on the common interface.
 *
 * Unknowns per Newton step are corrections:
 *
 *   du       fluid velocity correction
 *   dp       fluid pressure correction
 *   dlambda  pressure gauge correction
 *   deta     displacement increment correction
 *
 * The nonlinear state is:
 *
 *   uState      fluid velocity at t^{n+1}
 *   pState      pressure at t^{n+1}
 *   lambdaState pressure gauge variable
 *   etaState    displacement increment d^{n+1} - d^n
 *   dState      total displacement d^n + etaState
 *
 * The monolithic FSI kinematic constraint is imposed at correction level by
 *
 *   du = deta / dt  on Gamma_FSI.
 *
 * This requires the variational identification support from PR #276:
 *
 *   DirichletBC(du, (1.0 / dt) * deta).on(Boundary::FSI)
 *
 * Because this is a homogeneous correction constraint, the initial Newton
 * guess at each time step is chosen so that
 *
 *   uState = etaState / dt
 *
 * already holds on Gamma_FSI.  The correction constraint then preserves the
 * nonlinear kinematic relation.
 *
 * Notes:
 *
 * 1. This example assumes Integral(...).over(Attribute) and solid form
 *    restrictions such as solidTangent.over(Volume::Solid) are available.
 *    If Solid::MaterialTangent / Solid::InternalForce do not currently support
 *    .over(attr), add domain restriction to those forms or assemble them on a
 *    solid submesh.
 *
 * 2. DirichletBC identification is attribute-driven here: without `.on(...)`
 *    it scans the exterior boundary, while `.on(Boundary::FSI)` selects the
 *    tagged internal interface facets.
 *
 * 3. Velocity and pressure spaces below are global on the full mesh while their
 *    physical equations are restricted to the fluid.  The small inactive-region
 *    regularization terms avoid null rows in the solid part.  A cleaner future
 *    implementation should use subdomain-restricted spaces.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

#include <petscsys.h>

#include <Rodin/Types.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/PETSc.h>
#include <Rodin/Solid.h>
#include <Rodin/Solver.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Math;
using namespace Rodin::Solver;
using namespace Rodin::Variational;

namespace
{
  using LocalMesh = Mesh<Context::Local>;

  struct Volume
  {
    static constexpr Attribute Fluid = 1;
    static constexpr Attribute Solid = 2;
  };

  struct Boundary
  {
    static constexpr Attribute Inlet          = 10;
    static constexpr Attribute Outlet         = 11;
    static constexpr Attribute FluidWall      = 12;
    static constexpr Attribute SolidClampLeft = 13;
    static constexpr Attribute SolidClampRight = 14;
    static constexpr Attribute SolidFree      = 15;
    static constexpr Attribute FSI            = 16;
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

    return amplitude * std::sin(Constants::pi() * localTime / duration);
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

  /**
   * @brief Maps a uniform logical rectangle to a two-layer FSI strip.
   *
   * The interface row is forced to be exactly at y = Hf.  This avoids the
   * common bug where y = Hf does not coincide with a mesh row because
   * Hf / (Hf + Hs) differs from nfy / (nfy + nsy).
   */
  static void mapMeshToFSIStrip(
      LocalMesh& mesh,
      size_t nfy,
      size_t nsy,
      Real L,
      Real Hf,
      Real Hs)
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

    const Real rInterface =
      static_cast<Real>(nfy) / static_cast<Real>(nfy + nsy);

    for (auto it = mesh.getVertex(); it; ++it)
    {
      auto x = mesh.getVertexCoordinates(it->getIndex());

      const Real rx = (x.x() - xmin) / (xmax - xmin);
      const Real ry = (x.y() - ymin) / (ymax - ymin);

      x(0) = L * rx;

      if (ry <= rInterface)
      {
        x(1) = Hf * ry / rInterface;
      }
      else
      {
        x(1) = Hf + Hs * (ry - rInterface) / (Real(1) - rInterface);
      }

      mesh.setVertexCoordinates(it->getIndex(), x);
    }

    mesh.flush();
  }

  static LocalMesh makeFSIMesh(size_t nx, size_t nfy, size_t nsy, Real L, Real Hf, Real Hs)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { nx + 1, nfy + nsy + 1 });

    mapMeshToFSIStrip(mesh, nfy, nsy, L, Hf, Hs);

    mesh.getConnectivity().compute(1, 2);
    mesh.getConnectivity().compute(2, 1);

    constexpr Real eps = 1e-10;
    const Real yTop = Hf + Hs;

    for (auto it = mesh.getCell(); it; ++it)
    {
      const auto c = centroid(mesh, *it);
      if (c.y() < Hf)
        mesh.setAttribute(it.key(), Volume::Fluid);
      else
        mesh.setAttribute(it.key(), Volume::Solid);
    }

    /* Exterior boundary attributes. */
    for (auto it = mesh.getBoundary(); it; ++it)
    {
      const auto c = centroid(mesh, *it);

      if (std::abs(c.y()) < eps)
      {
        mesh.setAttribute(it.key(), Boundary::FluidWall);
      }
      else if (std::abs(c.x()) < eps && c.y() < Hf - eps)
      {
        mesh.setAttribute(it.key(), Boundary::Inlet);
      }
      else if (std::abs(c.x() - L) < eps && c.y() < Hf - eps)
      {
        mesh.setAttribute(it.key(), Boundary::Outlet);
      }
      else if (std::abs(c.x()) < eps && c.y() > Hf + eps)
      {
        mesh.setAttribute(it.key(), Boundary::SolidClampLeft);
      }
      else if (std::abs(c.x() - L) < eps && c.y() > Hf + eps)
      {
        mesh.setAttribute(it.key(), Boundary::SolidClampRight);
      }
      else if (std::abs(c.y() - yTop) < eps)
      {
        mesh.setAttribute(it.key(), Boundary::SolidFree);
      }
    }

    /* Internal interface facets.
     *
     * Faces touching x = 0 or x = L are left untagged. Their endpoint DOFs are
     * already governed by exterior essential conditions, so adding an FSI
     * identification there would over-constrain the correction system.
     */
    const size_t faceDim = mesh.getDimension() - 1;
    for (auto it = mesh.getFace(); it; ++it)
    {
      const auto c = centroid(mesh, *it);
      bool touchesExterior = false;
      for (const auto vertex : it->getVertices())
      {
        const auto x = mesh.getVertexCoordinates(vertex);
        touchesExterior =
          touchesExterior || std::abs(x.x()) < eps || std::abs(x.x() - L) < eps;
      }

      if (!touchesExterior && std::abs(c.y() - Hf) < eps)
        mesh.setAttribute({ faceDim, it->getIndex() }, Boundary::FSI);
    }

    return mesh;
  }
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  {
  PetscInt nxOpt = 48;
  PetscInt fluidNyOpt = 12;
  PetscInt solidNyOpt = 4;
  PetscInt ntOpt = 100;

  PetscReal finalTimeOpt = 1.0;
  PetscReal inletAmplitudeOpt = 1.25;
  PetscReal pulseDurationOpt = 0.18;
  PetscReal pulseGapOpt = 0.35;
  PetscInt pulseCountOpt = 1;

  PetscReal solidEOpt = 2.0e2;
  PetscReal solidNuOpt = 0.3;
  PetscReal solidDensityOpt = 5.0e-1;
  PetscReal solidRayleighAlphaOpt = 5.0e-2;
  PetscReal inactiveOpt = 1e-8;

  PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fsi_nx", &nxOpt, PETSC_NULLPTR);
  PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fsi_fluid_ny", &fluidNyOpt, PETSC_NULLPTR);
  PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fsi_solid_ny", &solidNyOpt, PETSC_NULLPTR);
  PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fsi_Nt", &ntOpt, PETSC_NULLPTR);
  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fsi_T", &finalTimeOpt, PETSC_NULLPTR);

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
      "-fsi_inactive_regularization", &inactiveOpt, PETSC_NULLPTR);

  const size_t nx = static_cast<size_t>(std::max<PetscInt>(3, nxOpt));
  const size_t nfy = static_cast<size_t>(std::max<PetscInt>(2, fluidNyOpt));
  const size_t nsy = static_cast<size_t>(std::max<PetscInt>(1, solidNyOpt));
  const Index Nt = static_cast<Index>(std::max<PetscInt>(1, ntOpt));

  const Real L = 4.0;
  const Real Hf = 1.0;
  const Real Hs = 0.25;
  const Real T = static_cast<Real>(finalTimeOpt);
  const Real dt = T / static_cast<Real>(Nt);

  LocalMesh mesh = makeFSIMesh(nx, nfy, nsy, L, Hf, Hs);
  const size_t dim = mesh.getSpaceDimension();

  /* Global spaces.
   *
   * Velocity/pressure are regularized in the inactive solid region below.
   * Displacement is physical in the solid and ALE mesh displacement in the
   * fluid.
   */
  H1 uh(std::integral_constant<size_t, 2>{}, mesh, dim);
  H1 ph(std::integral_constant<size_t, 1>{}, mesh);
  P0g lh(mesh);
  H1 dh(std::integral_constant<size_t, 2>{}, mesh, dim);

  PETSc::Variational::TrialFunction du(uh);
  PETSc::Variational::TrialFunction dp(ph);
  PETSc::Variational::TrialFunction dlambda(lh);
  PETSc::Variational::TrialFunction deta(dh);

  PETSc::Variational::TestFunction v(uh);
  PETSc::Variational::TestFunction q(ph);
  PETSc::Variational::TestFunction m(lh);
  PETSc::Variational::TestFunction z(dh);

  du.setName("du");
  dp.setName("dp");
  dlambda.setName("dlambda");
  deta.setName("deta");

  PETSc::Variational::GridFunction uState(uh);
  PETSc::Variational::GridFunction pState(ph);
  PETSc::Variational::GridFunction lambdaState(lh);
  PETSc::Variational::GridFunction etaState(dh);
  PETSc::Variational::GridFunction dState(dh);

  PETSc::Variational::GridFunction uOld(uh);
  PETSc::Variational::GridFunction pOld(ph);
  PETSc::Variational::GridFunction etaOld(dh);
  PETSc::Variational::GridFunction dOld(dh);

  PETSc::Variational::GridFunction meshVelocity(dh);
  PETSc::Variational::GridFunction inletVelocity(uh);

  uState.setName("FluidVelocity");
  pState.setName("FluidPressure");
  lambdaState.setName("PressureGauge");
  etaState.setName("DisplacementIncrement");
  dState.setName("Displacement");
  meshVelocity.setName("ALEMeshVelocity");

  uState = Zero(dim);
  pState = 0.0;
  lambdaState = 0.0;
  etaState = Zero(dim);
  dState = Zero(dim);
  uOld = Zero(dim);
  pOld = 0.0;
  etaOld = Zero(dim);
  dOld = Zero(dim);
  meshVelocity = Zero(dim);
  inletVelocity = Zero(dim);

  IO::XDMF xdmf("out/Seq_BDF1_Monolithic_ALE_FSI");
  auto grid = xdmf.grid("FSI");
  grid.setMesh(mesh, IO::XDMF::MeshPolicy::Static);
  grid.add("velocity", uState);
  grid.add("pressure", pState);
  grid.add("displacement", dState);
  grid.add("mesh_velocity", meshVelocity);
  xdmf.write(0.0);

  const Real rhoF = 1.0;
  const Real muF = 0.02;

  const Real solidE = static_cast<Real>(solidEOpt);
  const Real solidNu = static_cast<Real>(solidNuOpt);
  const Real solidLambda =
    solidE * solidNu / ((1.0 + solidNu) * (1.0 - 2.0 * solidNu));
  const Real solidMu = solidE / (2.0 * (1.0 + solidNu));
  const Real rhoS = static_cast<Real>(solidDensityOpt);
  const Real solidRayleighAlpha = static_cast<Real>(solidRayleighAlphaOpt);
  const Real inactive = static_cast<Real>(inactiveOpt);

  Solid::NeoHookean law(solidLambda, solidMu);

  Solid::MaterialTangent solidTangent(law, deta, z);
  solidTangent.setDisplacement(dState);

  Solid::InternalForce solidInternal(law, z);
  solidInternal.setDisplacement(dState);

  auto n = BoundaryNormal(mesh);

  /* Fully nonlinear BDF1 ALE Navier--Stokes fluid block.
   *
   * The ALE transport field is
   *
   *   beta(u, eta) = u - eta / dt,
   *
   * where eta is the displacement increment from t^n to t^{n+1}.  Therefore
   * the Newton correction of the transport field is
   *
   *   delta beta = du - deta / dt.
   *
   * The convective residual is written in skew-symmetric ALE form:
   *
   *   rho (beta . grad) u + 0.5 rho div(beta) u.
   *
   * Its Newton tangent is
   *
   *   rho (beta . grad) du
   * + rho (delta beta . grad) u
   * + 0.5 rho div(beta) du
   * + 0.5 rho div(delta beta) u.
   */
  const auto meshVelocityState = (1.0 / dt) * etaState;
  const auto transportVelocity = uState - meshVelocityState;

  const auto convState = Mult(Jacobian(uState), transportVelocity);
  const auto convTangentVelocity = Mult(Jacobian(du), transportVelocity);
  const auto convTangentTransportU = Mult(Jacobian(uState), du);
  const auto convTangentTransportEta = Mult(Jacobian(uState), deta);

  /* Div does not currently have a deduction guide for Sum expressions such as
   * Div(uState - etaState / dt), so write the divergence algebraically.
  */
  const auto divTransport = Div(uState) - (1.0 / dt) * Div(etaState);
  const auto beta = Max(-Dot(transportVelocity, n), 0.0);

  Problem fsi(du, dp, dlambda, deta, v, q, m, z);

  fsi =
      /* Fully nonlinear ALE Navier--Stokes tangent over fluid cells. */
        (rhoF / dt) * Integral(du, v).over(Volume::Fluid)
      + rhoF * Integral(Dot(convTangentVelocity, v)).over(Volume::Fluid)
      + rhoF * Integral(Dot(convTangentTransportU, v)).over(Volume::Fluid)
      - (rhoF / dt) * Integral(Dot(convTangentTransportEta, v)).over(Volume::Fluid)
      + 0.5 * rhoF * Integral(divTransport * Dot(du, v)).over(Volume::Fluid)
      + 0.5 * rhoF * Integral(Dot(Div(du) * uState, v)).over(Volume::Fluid)
      - (0.5 * rhoF / dt) * Integral(Dot(Div(deta) * uState, v)).over(Volume::Fluid)
      + muF * Integral(Jacobian(du), Jacobian(v)).over(Volume::Fluid)
      - Integral(dp, Div(v)).over(Volume::Fluid)
      + Integral(Div(du), q).over(Volume::Fluid)
      + Integral(dlambda, q).over(Volume::Fluid)
      + Integral(dp, m).over(Volume::Fluid)
      + BoundaryIntegral(0.5 * rhoF * beta * Dot(du, v)).over(Boundary::Outlet)
      + 1e-10 * Integral(dp, q).over(Volume::Fluid)
      + 1e-10 * Integral(dlambda, m)

      /* Fully nonlinear ALE Navier--Stokes residual over fluid cells. */
      + (rhoF / dt) * Integral(uState - uOld, v).over(Volume::Fluid)
      + rhoF * Integral(Dot(convState, v)).over(Volume::Fluid)
      + 0.5 * rhoF * Integral(Dot(divTransport * uState, v)).over(Volume::Fluid)
      + muF * Integral(Jacobian(uState), Jacobian(v)).over(Volume::Fluid)
      - Integral(pState, Div(v)).over(Volume::Fluid)
      + Integral(Div(uState), q).over(Volume::Fluid)
      + Integral(lambdaState, q).over(Volume::Fluid)
      + Integral(pState, m).over(Volume::Fluid)
      + BoundaryIntegral(0.5 * rhoF * Dot(beta * uState, v)).over(Boundary::Outlet)
      + 1e-10 * Integral(pState, q).over(Volume::Fluid)
      + 1e-10 * Integral(lambdaState, m)

      /* Harmonic ALE displacement in the fluid. */
      + Integral(Jacobian(deta), Jacobian(z)).over(Volume::Fluid)
      + Integral(Jacobian(etaState), Jacobian(z)).over(Volume::Fluid)

      /* Solid BDF1/Newmark-like displacement increment equation. */
      + (rhoS / (dt * dt)) * Integral(deta, z).over(Volume::Solid)
      + (solidRayleighAlpha * rhoS / dt) * Integral(deta, z).over(Volume::Solid)
      + solidTangent.over(Volume::Solid)

      + (rhoS / (dt * dt)) * Integral(etaState - etaOld, z).over(Volume::Solid)
      + (solidRayleighAlpha * rhoS / dt) * Integral(etaState, z).over(Volume::Solid)
      + solidInternal.over(Volume::Solid)

      /* Temporary inactive-region regularization for globally-defined fluid
       * variables.  Replace by subdomain-restricted spaces when available.
       */
      + inactive * Integral(du, v).over(Volume::Solid)
      + inactive * Integral(uState, v).over(Volume::Solid)
      + inactive * Integral(dp, q).over(Volume::Solid)
      + inactive * Integral(pState, q).over(Volume::Solid)

      /* Boundary and interface constraints. */
      + DirichletBC(du, inletVelocity - uState).on(Boundary::Inlet)
      + DirichletBC(du, -uState).on(Boundary::FluidWall)
      + DirichletBC(deta, -etaState).on(Boundary::Inlet, Boundary::Outlet, Boundary::FluidWall)
      + DirichletBC(deta, -etaState).on(Boundary::SolidClampLeft, Boundary::SolidClampRight)
      + DirichletBC(du, (1.0 / dt) * deta).on(Boundary::FSI);

  fsi.assemble().setFieldSplits();

  Solver::KSP ksp(fsi);
  Solver::SNES snes(ksp);

  snes.setTolerances(1e-10, 1e-8, 1e-10, 30, 10000);

  snes.setStateUpdate([&](const PETSc::Math::Vector& state)
  {
    const size_t uOffset = 0;
    const size_t pOffset = uOffset + uh.getSize();
    const size_t lambdaOffset = pOffset + ph.getSize();
    const size_t etaOffset = lambdaOffset + lh.getSize();

    uState.setData(state, uOffset);
    pState.setData(state, pOffset);
    lambdaState.setData(state, lambdaOffset);
    etaState.setData(state, etaOffset);

    dState = dOld + etaState;
    meshVelocity = (1.0 / dt) * etaState;
  });

  Real t = 0.0;

  for (Index step = 1; step <= Nt; ++step)
  {
    t += dt;

    inletVelocity = VectorFunction(dim, [&](const Point& x)
    {
      const Real y = std::clamp(x.y(), 0.0, Hf);
      const Real profile = 4.0 * y * (Hf - y) / (Hf * Hf);
      const Real amplitude = pulseTrain(
          t,
          static_cast<Real>(inletAmplitudeOpt),
          static_cast<Real>(pulseDurationOpt),
          static_cast<Real>(pulseGapOpt),
          static_cast<Index>(pulseCountOpt));

      Math::SpatialVector<Real> value(dim);
      value.setZero();
      value(0) = amplitude * profile;
      return value;
    });

    /* Initial Newton guess.
     *
     * The homogeneous correction constraint du = deta / dt preserves the
     * nonlinear relation only if the initial state already satisfies it on the
     * interface.  Use etaState = dt * uOld as a simple consistent interface
     * guess, then Newton corrections preserve du = deta / dt.
     */
    uState = uOld;
    pState = pOld;
    lambdaState = 0.0;
    etaState = dt * uOld;
    dState = dOld + etaState;
    meshVelocity = (1.0 / dt) * etaState;

    Alert::Info()
      << "Monolithic BDF1 ALE FSI step " << step << " / " << Nt
      << Alert::Raise;

    snes.solve();

    if (!snes.converged())
    {
      Alert::MemberFunctionException("Seq_BDF1_Monolithic_ALE_FSI", "main")
        << "SNES failed to converge at step " << step << Alert::Raise;
    }

    uOld.setData(uState.getData());
    pOld.setData(pState.getData());
    etaOld.setData(etaState.getData());
    dOld.setData(dState.getData());

    xdmf.write(t);
    xdmf.flush();
  }
  }

  PetscFinalize();
  return 0;
}
