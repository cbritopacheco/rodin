/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @example Transient 2D nonlinear lid-driven cavity with Taylor-Hood elements and PETSc SNES
 *
 * This example solves the transient incompressible Navier-Stokes equations in
 * the unit square with a fully implicit nonlinear solve at every time step.
 * It mirrors the Oseen lid-driven cavity example, but PETSc SNES solves the
 * nonlinear backward-Euler problem instead of freezing the transport velocity.
 *
 * At each time step, given u^n, SNES solves for (u^{n+1}, p^{n+1}):
 *
 *   (rho / dt) (u^{n+1} - u^n, v)
 * + rho ((u^{n+1} . grad) u^{n+1}, v)
 * + mu (grad u^{n+1}, grad v)
 * - (p^{n+1}, div v)
 * + (div u^{n+1}, q)
 * + eps_p (p^{n+1}, q)
 * = 0.
 *
 * The nonlinear convection Jacobian includes both Newton terms:
 *
 *   ((u . grad) du, v) + ((du . grad) u, v).
 *
 * Suggested quick run:
 *
 *   ./examples/PETSc/PDEs/PETSc_Seq_NavierStokes \
 *     -ns_n 12 -ns_nt 3 -ns_T 0.03 \
 *     -snes_monitor -snes_converged_reason \
 *     -ksp_type preonly -pc_type lu
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

#include <petscksp.h>
#include <petscsnes.h>

#include <Rodin/PETSc.h>

#include <Rodin/Alert.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Solver.h>
#include <Rodin/Types.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Math;
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  void labelUnitSquareBoundary(LocalMesh& mesh)
  {
    constexpr Attribute left   = 1;
    constexpr Attribute right  = 2;
    constexpr Attribute bottom = 3;
    constexpr Attribute top    = 4;
    constexpr Real eps = 1e-6;

    for (auto it = mesh.getFace(); it; ++it)
    {
      if (!mesh.isBoundary(it->getIndex()))
        continue;

      const auto& verts = it->getVertices();
      Math::SpatialPoint xc{0.0, 0.0};
      for (const auto& vtx : verts)
        xc += mesh.getVertexCoordinates(vtx);
      xc /= static_cast<Real>(verts.size());

      Attribute attr = 0;
      if (std::abs(xc.x()) < eps)
        attr = left;
      else if (std::abs(xc.x() - 1.0) < eps)
        attr = right;
      else if (std::abs(xc.y()) < eps)
        attr = bottom;
      else if (std::abs(xc.y() - 1.0) < eps)
        attr = top;
      else
        continue;

      mesh.setAttribute(it.key(), attr);
    }
  }

  void copyIntoStateVector(::Vec state, ::Vec values, PetscInt offset)
  {
    PetscErrorCode ierr;

    PetscInt n = 0;
    ierr = VecGetSize(values, &n);
    assert(ierr == PETSC_SUCCESS);

    ::IS is = PETSC_NULLPTR;
    ierr = ISCreateStride(PETSC_COMM_SELF, n, offset, 1, &is);
    assert(ierr == PETSC_SUCCESS);

    ::Vec sub = PETSC_NULLPTR;
    ierr = VecGetSubVector(state, is, &sub);
    assert(ierr == PETSC_SUCCESS);

    ierr = VecCopy(values, sub);
    assert(ierr == PETSC_SUCCESS);

    ierr = VecRestoreSubVector(state, is, &sub);
    assert(ierr == PETSC_SUCCESS);

    ierr = ISDestroy(&is);
    assert(ierr == PETSC_SUCCESS);
  }
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  constexpr Attribute left   = 1;
  constexpr Attribute right  = 2;
  constexpr Attribute bottom = 3;
  constexpr Attribute top    = 4;

  PetscInt nOpt = 64;
  PetscInt ntOpt = 1000;
  PetscReal finalTime = 10.0;
  PetscReal rhoOpt = 1.0;
  PetscReal muOpt = 0.01;
  PetscReal lidSpeedOpt = 1.0;

  PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-ns_n", &nOpt, PETSC_NULLPTR);
  PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-ns_nt", &ntOpt, PETSC_NULLPTR);
  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-ns_T", &finalTime, PETSC_NULLPTR);
  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-ns_rho", &rhoOpt, PETSC_NULLPTR);
  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-ns_mu", &muOpt, PETSC_NULLPTR);
  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-ns_lid_speed", &lidSpeedOpt, PETSC_NULLPTR);

  const size_t n = static_cast<size_t>(std::max<PetscInt>(nOpt, 2));
  const Index nt = static_cast<Index>(std::max<PetscInt>(ntOpt, 1));
  const Real T = static_cast<Real>(finalTime);
  const Real dt = T / static_cast<Real>(nt);
  const Real rho = static_cast<Real>(rhoOpt);
  const Real mu = static_cast<Real>(muOpt);
  const Real lidSpeed = static_cast<Real>(lidSpeedOpt);
  constexpr Real epsP = 1e-12;

  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Quadrilateral, { n, n });
  mesh.getConnectivity().compute(1, 2);
  mesh.scale(1.0 / (n - 1));
  labelUnitSquareBoundary(mesh);

  Alert::Info()
    << "Setting up nonlinear NavierStokes cavity with n = " << n
    << ", Nt = " << nt << ", dt = " << dt
    << Alert::Raise;

  IO::XDMF xdmf("NavierStokes");
  xdmf.setMesh(mesh);

  const size_t dim = mesh.getSpaceDimension();
  H1 uh(std::integral_constant<size_t, 2>{}, mesh, dim);
  H1 ph(std::integral_constant<size_t, 1>{}, mesh);

  int exitCode = 0;

  {
    PETSc::Variational::GridFunction uState(uh);
    PETSc::Variational::GridFunction pState(ph);
    PETSc::Variational::GridFunction uOld(uh);
    PETSc::Variational::GridFunction pOld(ph);
    PETSc::Variational::GridFunction lidVelocity(uh);

    uState = Zero(dim);
    pState = 0.0;
    uOld = Zero(dim);
    pOld = 0.0;
    lidVelocity = VectorFunction{lidSpeed, 0.0};

    xdmf.add("velocity", uState);
    xdmf.add("pressure", pState);
    xdmf.write(0.0).flush();

    PETSc::Variational::TrialFunction du(uh); du.setName("u");
    PETSc::Variational::TrialFunction dp(ph); dp.setName("p");
    PETSc::Variational::TestFunction  v(uh);
    PETSc::Variational::TestFunction  q(ph);

    const auto stateConvection = Mult(Jacobian(uState), uState);
    const auto newtonConvection =
      Mult(Jacobian(du), uState) + Mult(Jacobian(uState), du);

    Problem navierStokes(du, dp, v, q);
    navierStokes =
        (rho / dt) * Integral(du, v)
      + rho * Integral(Dot(newtonConvection, v))
      + mu * Integral(Jacobian(du), Jacobian(v))
      - Integral(dp, Div(v))
      + Integral(Div(du), q)
      + epsP * Integral(dp, q)

      + (rho / dt) * Integral(uState, v)
      - (rho / dt) * Integral(uOld, v)
      + rho * Integral(Dot(stateConvection, v))
      + mu * Integral(Jacobian(uState), Jacobian(v))
      - Integral(pState, Div(v))
      + Integral(Div(uState), q)
      + epsP * Integral(pState, q)

      + DirichletBC(du, -uState).on(left)
      + DirichletBC(du, -uState).on(right)
      + DirichletBC(du, -uState).on(bottom)
      + DirichletBC(du, lidVelocity - uState).on(top);

    navierStokes.assemble().setFieldSplits();

    ::Vec x = PETSC_NULLPTR;
    PetscErrorCode ierr = VecDuplicate(navierStokes.getLinearSystem().getSolution(), &x);
    assert(ierr == PETSC_SUCCESS);

    auto syncState = [&](::Vec state)
    {
      uState.setData(state, 0);
      pState.setData(state, uh.getSize());
    };

    Solver::KSP ksp(navierStokes);
    ksp.setType(KSPPREONLY);

    ::PC pc = PETSC_NULLPTR;
    ierr = KSPGetPC(ksp.getHandle(), &pc);
    assert(ierr == PETSC_SUCCESS);
    ierr = PCSetType(pc, PCLU);
    assert(ierr == PETSC_SUCCESS);
    ierr = PCFactorSetShiftType(pc, MAT_SHIFT_NONZERO);
    assert(ierr == PETSC_SUCCESS);
    ierr = PCFactorSetShiftAmount(pc, 1e-10);
    assert(ierr == PETSC_SUCCESS);

    Solver::SNES snes(ksp);
    snes.setType(SNESNEWTONLS)
        .setTolerances(1e-10, 1e-8, 1e-10, 25, 1000)
        .setStateUpdate(syncState);

    for (Index step = 1; step <= nt; ++step)
    {
      const Real t = static_cast<Real>(step) * dt;
      const Real ramp = std::min<Real>(1.0, t / std::max<Real>(dt, 0.1 * T));
      lidVelocity = VectorFunction{lidSpeed * ramp, 0.0};

      copyIntoStateVector(x, uOld.getData(), 0);
      copyIntoStateVector(x, pOld.getData(), static_cast<PetscInt>(uh.getSize()));
      syncState(x);

      Alert::Info()
        << "SNES solve for nonlinear Navier-Stokes step "
        << step << " / " << nt << " at t = " << t
        << Alert::Raise;

      snes.solve(x);

      PetscInt its = 0;
      ::SNESConvergedReason reason;
      ierr = SNESGetIterationNumber(snes.getHandle(), &its);
      assert(ierr == PETSC_SUCCESS);
      ierr = SNESGetConvergedReason(snes.getHandle(), &reason);
      assert(ierr == PETSC_SUCCESS);

      if (reason < 0)
      {
        std::cerr << "SNES failed at step " << step
                  << " with reason " << reason << "\n";
        exitCode = 1;
        break;
      }

      Alert::Success()
        << "Step " << step << " converged in " << its
        << " SNES iterations."
        << Alert::Raise;

      uOld.setData(uState.getData());
      pOld.setData(pState.getData());

      xdmf.write(t).flush();
    }

    if (exitCode == 0)
    {
      mesh.save("NavierStokes.mesh", IO::FileFormat::MFEM);
      uState.save("NavierStokes_velocity.gf", IO::FileFormat::MFEM);
      pState.save("NavierStokes_pressure.gf", IO::FileFormat::MFEM);
    }

    ierr = VecDestroy(&x);
    assert(ierr == PETSC_SUCCESS);
  }

  xdmf.close();
  PetscFinalize();
  return exitCode;
}
