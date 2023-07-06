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
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int, char**)
{
  // Build a mesh
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Geometry::Triangle, 256, 256);
  mesh.getConnectivity().compute(1, 2);

  // Functions
  P1 vh(mesh);

  TrialFunction u(vh);
  TestFunction  v(vh);

  // Define problem
  ScalarFunction f = 1.0;
  ScalarFunction g = 5.0;

  Problem poisson(u, v);
  poisson = Integral(f * Grad(u), Grad(v))
          - Integral(v)
          + DirichletBC(u, g);
  poisson.assemble();

  // Solve the problem
  Solver::SparseLU solver;
  poisson.solve(solver);

  // Save solution
  u.getSolution().save("Poisson.gf");
  mesh.save("Poisson.mesh");

  return 0;
}
