/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int, char**)
{
  // Build a mesh
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 16, 16 });
  mesh.getConnectivity().compute(1, 2);

  // Functions
  P1<Complex> vh(mesh);
  GridFunction gf(vh);

  TrialFunction u(vh);
  TestFunction  v(vh);

  Problem helmholtz(u, v);
  helmholtz = Integral(Grad(u), Grad(v))
            - Integral(u, v)
            + DirichletBC(u, Zero());

  CG(helmholtz).solve();
  // poisson.solve(solver);

  // Save solution
  // u.getSolution().save("Poisson.gf");
  // mesh.save("Poisson.mesh");

  return 0;
}

