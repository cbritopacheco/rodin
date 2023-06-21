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
  // Load mesh
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Geometry::Triangle, 256, 256);
  mesh.getConnectivity().compute(1, 2);

  // Functions
  P1 vh(mesh);

  TrialFunction u(vh);
  TestFunction  v(vh);

  // Define problem
  ScalarFunction f = 1.0;
  ScalarFunction g = 0.0;

  // std::cout << "lf\n";
  // LinearForm lf(v);
  // lf = Integral(f, v);
  // lf.assemble();

  // std::cout << "bf\n";
  // BilinearForm bf(u, v);
  // bf = Integral(Grad(u), Grad(v));
  // bf.assemble();

  Problem poisson(u, v);
  poisson = Integral(Grad(u), Grad(v))
          - Integral(f, v)
          + DirichletBC(u, g);

  // // for (int i = 0; i < 25; i++)
  // poisson.assemble();


  // Solve the problem
  Solver::CG solver;
  solver.setMaxIterations(200).setTolerance(1e-12);
  poisson.solve(solver);

  // std::cout << poisson.getStiffnessOperator() << std::endl;

  // Save solution
  u.getSolution().save("Poisson.gf");
  mesh.save("Poisson.mesh");

  return 0;
}
