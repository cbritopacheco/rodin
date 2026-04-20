/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * Recommended PETSc options
 * -------------------------
 * Direct-solver diagnostic path:
 *
 *  -ksp_type preonly
 *  -pc_type lu
 *  -pc_factor_shift_type nonzero
 *  -pc_factor_shift_amount 1e-10
 *
 * This is useful to verify that the assembled mixed system is solvable and that
 * any remaining issues are more likely due to iterative preconditioning than to
 * the variational formulation itself.
 *
 * Performant direct solve:
 *
 *  -ksp_type preonly
 *  -pc_type lu
 *  -pc_factor_mat_solver_type mumps
 *  -mat_mumps_icntl_7 7
 *  -ksp_converged_reason
 *  -pc_factor_shift_type nonzero
 *  -pc_factor_shift_amount 1e-10
 *
 * Iterative solver path for mixed systems:
 *
 *  -ksp_type fgmres
 *  -pc_type fieldsplit
 *  -pc_fieldsplit_type schur
 *  -pc_fieldsplit_schur_fact_type lower
 *  -pc_fieldsplit_schur_precondition selfp
 *  -fieldsplit_u_ksp_type preonly
 *  -fieldsplit_u_pc_type gamg
 *  -fieldsplit_p_ksp_type preonly
 *  -fieldsplit_p_pc_type jacobi
 *  -ksp_rtol 1e-6
 *  -ksp_monitor
 *  -ksp_converged_reason
 *
 * The direct-solver path is recommended first when debugging a new example.
 * Once that works, the Schur-complement fieldsplit path is the more natural
 * iterative strategy for the mixed velocity-pressure system.
 */

#include "Rodin/Geometry/Types.h"
#include "Rodin/IO/ForwardDecls.h"
#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Variational.h>

#include <Rodin/PETSc.h>

using namespace Rodin;
using namespace Rodin::Math;
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  Attribute wall = 2;
  Attribute inlet = 3;
  FlatSet<Attribute> outlet{ 4, 5, 6, 7, 8, 9 };

  const Real eps = 1e-12;

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

  Mesh mesh;
  mesh.load("CoronaryArtery.mesh", IO::FileFormat::MEDIT);

  Alert::Info() << "Computing connectivity for CoronaryArtery.mesh ..." << Alert::Raise;

  mesh.getConnectivity().compute(3, 2);
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);

  mesh.getConnectivity().compute(2, 3);

  Alert::Info() << "Setting up CoronaryArtery.xdmf ..." << Alert::Raise;

  IO::XDMF xdmf("CoronaryArtery");
  xdmf.setMesh(mesh);

  const size_t dim = mesh.getSpaceDimension();

  H1 uh(std::integral_constant<size_t, 2>{}, mesh, dim);
  H1 ph(std::integral_constant<size_t, 1>{}, mesh);

  {
    // Unknowns and test functions.
    PETSc::Variational::TrialFunction u(uh); u.setName("u");
    PETSc::Variational::TrialFunction p(ph); p.setName("p");
    PETSc::Variational::TestFunction  v(uh);
    PETSc::Variational::TestFunction  q(ph);

    // Solution fields from the previous time step.
    PETSc::Variational::GridFunction u_old(uh);
    PETSc::Variational::GridFunction p_old(ph);

    xdmf.add("velocity", u.getSolution());
    xdmf.add("pressure", p.getSolution());

    // Initial condition: fluid initially at rest, zero pressure.
    u_old = Math::Vector<Real>{{0.0, 0.0, 0.0}};
    p_old = 0.0;

    const auto n = BoundaryNormal(mesh);

    const auto conv_u = Mult(Jacobian(u), u_old);
    const auto div_u_old = Div(u_old);

    // Time loop.
    Real t = 0.0;

    for (Index k = 0; k < Nt; k++)
    {
      const auto pulse = Max(Sin(RealFunction(2 * M_PI * t * 50.0 / 60.0)), 0.0);
      const auto in = -pulse * n;

      Problem flow(u, p, v, q);
      flow =
          (rho / dt) * Integral(u, v)
        - (rho / dt) * Integral(u_old, v)
        + rho * Integral(Dot(conv_u, v))
        + 0.5 * rho * Integral(div_u_old * Dot(u, v))
        + mu * Integral(Jacobian(u), Jacobian(v))
        - Integral(p, Div(v))
        + Integral(Div(u), q)
        + eps * Integral(p, q)

        + BoundaryIntegral(0.5 * rho * Max(-Dot(u_old, n), 0.0) * Dot(u, v)).over(outlet)
        + DirichletBC(u, Zero(dim)).on(wall)
        + DirichletBC(u, in).on(inlet);

      Alert::Info() << "Assembling time step " << k + 1 << " / " << Nt << "..."
        << Alert::Raise;

      flow.assemble().setFieldSplits();

      Alert::Info() << "Solving time step " << k + 1 << " / " << Nt << "..."
        << Alert::Raise;

      Solver::KSP(flow).solve();

      // Advance the time history.
      u_old.setData(u.getSolution().getData());
      p_old.setData(p.getSolution().getData());

      xdmf.write(t).flush();

      t += dt;
    }
  }

  PetscFinalize();
  return 0;
}

