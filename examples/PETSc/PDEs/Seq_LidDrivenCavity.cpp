/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @example Transient 2D lid-driven cavity with Taylor–Hood elements and PETSc
 *
 * This example solves the transient incompressible Navier–Stokes equations in
 * the unit square using:
 *   - a mixed finite element discretization,
 *   - Taylor–Hood spaces (quadratic velocity, linear pressure),
 *   - backward Euler time stepping,
 *   - a semi-implicit Oseen / Picard linearization of the convective term,
 *   - PETSc for linear algebra and iterative solution.
 *
 * Problem solved
 * --------------
 * Let Ω = [0,1]^2. At each time step, given the previous velocity u^n, the code
 * solves for (u^{n+1}, p^{n+1}) such that, for all test functions (v, q),
 *
 *     (rho / dt) (u^{n+1}, v)
 *   + rho ((u^n · ∇) u^{n+1}, v)
 *   + (rho / 2) ((∇ · u^n) u^{n+1}, v)
 *   + mu (∇u^{n+1}, ∇v)
 *   - (p^{n+1}, ∇ · v)
 *   + (∇ · u^{n+1}, q)
 *   + eps_p (p^{n+1}, q)
 *
 *   = (rho / dt) (u^n, v),
 *
 * where:
 *   - rho is the density,
 *   - mu is the viscosity,
 *   - dt is the time step,
 *   - the convective velocity is frozen at the previous time step u^n,
 *   - the additional (∇·u^n) term yields the skew-symmetric Oseen form,
 *   - eps_p is a tiny pressure-block diagonal filler used only to improve
 *     robustness of sparse direct factorization for the mixed system.
 *
 * Boundary conditions
 * -------------------
 * The cavity is driven by a moving lid:
 *   - top boundary:    u = (1, 0),
 *   - left boundary:   u = (0, 0),
 *   - right boundary:  u = (0, 0),
 *   - bottom boundary: u = (0, 0).
 *
 * This is the standard lid-driven cavity benchmark for incompressible flow.
 * Even at moderate Reynolds numbers, the moving top wall generates a primary
 * recirculation zone and nontrivial vorticity in the cavity.
 *
 * Physical interpretation
 * -----------------------
 * The code models viscous incompressible flow in a rigid square cavity with:
 *   - no-slip walls,
 *   - one moving boundary acting as the flow driver,
 *   - no inlet/outlet boundaries,
 *   - no external body force.
 *
 * Numerically, this is not a fully nonlinear implicit Navier–Stokes solve at
 * each time step. Instead, one Oseen solve is performed per time level using
 * the previous velocity as the transport field.
 *
 * Discretization
 * --------------
 *   - Velocity space: H1 of degree 2
 *   - Pressure space: H1 of degree 1
 *
 * This Taylor–Hood pair is inf-sup stable, so no PSPG/SUPG stabilization is
 * required for pressure stability in this example.
 *
 * Parameters
 * ----------
 * The default parameters are:
 *   - rho = 1
 *   - mu  = 0.01
 *   - lid speed = 1
 *   - T   = 5
 *   - Nt  = 500
 *
 * With cavity length scale L = 1 and lid speed U = 1, the kinematic viscosity
 * is nu = mu / rho = 0.01, so the nominal Reynolds number is:
 *
 *   Re = U L / nu = 100.
 *
 * This is a standard moderate-Reynolds-number regime for a first cavity test.
 *
 * Notes
 * -----
 *   - The pressure field is only defined up to a constant at the continuous
 *     level. The tiny term eps_p (p, q) is included here only as a numerical
 *     aid for sparse direct factorization in PETSc.
 *   - The boundary attribute numbering used below may depend on the conventions
 *     of the generated mesh. If the lid or wall conditions appear to be applied
 *     to the wrong boundaries, check the attribute labels of the mesh.
 *
 * Output
 * ------
 * The program writes:
 *
 *   LidDrivenCavity.mesh
 *   LidDrivenCavity_velocity_XXXXXX.gf
 *   LidDrivenCavity_pressure_XXXXXX.gf
 *   cavity_centerline.txt
 *
 * Recommended PETSc options
 * -------------------------
 * Direct-solver diagnostic path:
 *
 *   -ksp_type preonly
 *   -pc_type lu
 *   -pc_factor_shift_type nonzero
 *   -pc_factor_shift_amount 1e-10
 *
 * This is useful to verify that the assembled mixed system is solvable and that
 * any remaining issues are more likely due to iterative preconditioning than to
 * the variational formulation itself.
 *
 * Iterative solver path for mixed systems:
 *
 *   -ksp_type fgmres
 *   -pc_type fieldsplit
 *   -pc_fieldsplit_type schur
 *   -pc_fieldsplit_schur_fact_type lower
 *   -pc_fieldsplit_schur_precondition selfp
 *   -fieldsplit_u_ksp_type preonly
 *   -fieldsplit_u_pc_type gamg
 *   -fieldsplit_p_ksp_type preonly
 *   -fieldsplit_p_pc_type jacobi
 *   -ksp_rtol 1e-6
 *   -ksp_monitor
 *   -ksp_converged_reason
 *
 * The direct-solver path is recommended first when debugging a new example.
 * Once that works, the Schur-complement fieldsplit path is the more natural
 * iterative strategy for the mixed velocity-pressure system.
 */

#include "Rodin/Alert/Info.h"
#include "Rodin/Alert/Raise.h"
#include "Rodin/Math/SpatialVector.h"
#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/PETSc.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace Rodin;
using namespace Rodin::Math;
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  // Boundary tags for the generated unit-square mesh.
  //
  // Important:
  // The actual numbering may depend on the mesh-generation convention used by
  // UniformGrid. If the lid or wall conditions are applied to the wrong sides,
  // inspect the boundary attributes of the generated mesh and adjust these tags.
  constexpr Attribute left   = 1;
  constexpr Attribute right  = 2;
  constexpr Attribute bottom = 3;
  constexpr Attribute top    = 4;

  size_t n = 64;
  const Real eps = 1e-6;

  // Physical parameters.
  //
  // rho = density
  // mu  = dynamic viscosity
  //
  // With U_lid = 1 and cavity size L = 1, the Reynolds number is approximately
  // Re = rho * U * L / mu = 100 for the chosen values below.
  const Real rho = 1.0;
  const Real mu  = 0.01;

  // Time-stepping parameters.
  const Real T   = 5.0;
  const Index Nt = 500;
  const Real dt  = T / Nt;

  // Structured unit-square mesh, triangulated.
  //
  // The mesh is first generated on an integer grid and then scaled to [0,1]^2.
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
  mesh.getConnectivity().compute(1, 2);
  mesh.scale(1.0 / (n - 1));

  for (auto it = mesh.getFace(); it; ++it)
  {
    // Skip interior faces: only boundary faces should be labeled.
    if (!mesh.isBoundary(it->getIndex()))
      continue;

    const auto& verts = it->getVertices();

    Math::SpatialPoint xc{0.0, 0.0};
    for (const auto& v : verts)
      xc += mesh.getVertexCoordinates(v);

    xc /= static_cast<Real>(verts.size());

    Attribute attr = 0;

    if (std::abs(xc.x() - 0.0) < eps)
      attr = left;
    else if (std::abs(xc.x() - 1.0) < eps)
      attr = right;
    else if (std::abs(xc.y() - 0.0) < eps)
      attr = bottom;
    else if (std::abs(xc.y() - 1.0) < eps)
      attr = top;
    else
      continue;

    mesh.setAttribute(it.key(), attr);
  }

  mesh.save("LidDrivenCavity.mesh", IO::FileFormat::MEDIT);

  const size_t dim = mesh.getSpaceDimension();

  // Taylor–Hood spaces:
  //   - velocity: H1 degree 2, vector-valued
  //   - pressure: H1 degree 1, scalar-valued
  H1 uh(std::integral_constant<size_t, 2>{}, mesh, dim);
  H1 ph(std::integral_constant<size_t, 1>{}, mesh);

  // Unknowns and test functions.
  PETSc::Variational::TrialFunction u(uh); u.setName("u");
  PETSc::Variational::TrialFunction p(ph); p.setName("p");
  PETSc::Variational::TestFunction  v(uh);
  PETSc::Variational::TestFunction  q(ph);

  // Solution fields from the previous time step.
  PETSc::Variational::GridFunction u_old(uh);
  PETSc::Variational::GridFunction p_old(ph);

  // Initial condition: fluid initially at rest, zero pressure.
  u_old = Math::Vector<Real>{{0.0, 0.0}};
  p_old = 0.0;

  // Prescribed lid velocity.
  const VectorFunction lidVelocity{1.0, 0.0};

  // Time loop.
  for (Index k = 0; k < Nt; k++)
  {
    // Oseen / Picard linearization:
    //   (u_old · ∇)u
    const auto conv_u = Mult(Jacobian(u), u_old);

    // Divergence of the frozen transport velocity.
    //
    // Included in the skew-symmetric linearization:
    //   ((u_old · ∇)u, v) + 0.5 ((div u_old) u, v)
    const auto div_u_old = Div(u_old);

    // Mixed transient Oseen problem at the current time step.
    Problem flow(u, p, v, q);
    flow =
        // Backward Euler time derivative:
        //   (rho / dt) (u - u_old, v)
        (rho / dt) * Integral(u, v)
      - (rho / dt) * Integral(u_old, v)

        // Linearized convection:
        //   rho ((u_old · ∇)u, v)
      + rho * Integral(Dot(conv_u, v))

        // Skew-symmetric correction:
        //   (rho / 2) ((div u_old) u, v)
      + 0.5 * rho * Integral(div_u_old * Dot(u, v))

        // Viscous term:
        //   mu (∇u, ∇v)
      + mu * Integral(Jacobian(u), Jacobian(v))

        // Pressure-velocity coupling:
        //   -(p, div v)
      - Integral(p, Div(v))

        // Incompressibility constraint:
        //   (div u, q)
      + Integral(Div(u), q)

        // Tiny pressure-block diagonal filler.
        //
        // This is not pressure stabilization in the inf-sup sense. It is only a
        // numerical aid for sparse direct factorization of the mixed system.
      + 1e-12 * Integral(p, q)

        // No-slip boundary conditions on the cavity walls.
      + DirichletBC(u, Zero(dim)).on(left)
      + DirichletBC(u, Zero(dim)).on(right)
      + DirichletBC(u, Zero(dim)).on(bottom)

        // Moving lid on the top boundary.
      + DirichletBC(u, lidVelocity).on(top);

    Alert::Info() << "Assembling time step " << k + 1 << " / " << Nt << "..."
      << Alert::Raise;
    flow.assemble().setFieldSplits();

    Alert::Info() << "Solving time step " << k + 1 << " / " << Nt << "..."
      << Alert::Raise;
    Solver::KSP(flow).solve();

    // Advance the time history.
    u_old.setData(u.getSolution().getData());
    p_old.setData(p.getSolution().getData());

    // Save solution snapshots.
    mesh.save("LidDrivenCavity." + std::to_string(k) + ".mesh", IO::FileFormat::MEDIT);
    u_old.save("LidDrivenCavity." + std::to_string(k) + ".sol", IO::FileFormat::MEDIT);
    // p_old.save("LidDrivenCavity_pressure_" + std::to_string(k) + ".gf", IO::FileFormat::MEDIT);
  }

  PetscFinalize();
  return 0;
}
