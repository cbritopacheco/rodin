/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Geometry/Polytope.h"
#include <Rodin/PETSc.h>

#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Math;
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, { 16, 16, 16 });
  mesh.scale(1.0 / (16 - 1));
  mesh.getConnectivity().compute(2, 3);
  mesh.getConnectivity().compute(3, 2);
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);

  RealFunction f = 1;

  H1 uh(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
  P1 ph(mesh);

  P0g p0g(mesh);

  {
    VectorFunction u_exact{ F::y, F::z, F::x };

    auto p_exact = Zero();

    VectorFunction f{ Zero(), Zero(), Zero() };

    // Define problem
    PETSc::Variational::TrialFunction u(uh);
    PETSc::Variational::TrialFunction p(ph);
    PETSc::Variational::TestFunction  v(uh);
    PETSc::Variational::TestFunction  q(ph);

    PETSc::Variational::TrialFunction lambda(p0g);
    PETSc::Variational::TestFunction  mu(p0g);

    Problem stokes(u, p, v, q, lambda, mu);
    stokes = Integral(Jacobian(u), Jacobian(v))
           - Integral(p, Div(v))
           + Integral(Div(u), q)
           + Integral(lambda, q)
           + Integral(p, mu)
           - Integral(f, v)
           + DirichletBC(u, u_exact);

    GMRES gmres(stokes);
    gmres.solve();

    // Save solution
    u.getSolution().save("Stokes.gf");
    mesh.save("Stokes.mesh");
  }

  PetscFinalize();

  return 0;
}


