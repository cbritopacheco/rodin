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
  const char* meshFile = "../resources/mfem/StarSquare.mfem.mesh";

  // Define boundary attributes
  int Gamma = 1;

  // Load mesh
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Geometry::Triangle, 16, 16);
  mesh.scale(1. / 15);
  mesh.getConnectivity().compute(1, 2);

  for (auto it = mesh.getBoundary(); !it.end(); ++it)
    mesh.setAttribute(it->getDimension(), it->getIndex(), 1);

  // Functions
  P1 vh(mesh);

  TrialFunction u(vh);
  TestFunction  v(vh);

  // Define problem
  ScalarFunction f = 1.0;
  ScalarFunction g = 0.0;

  Problem poisson(u, v);
  poisson = Integral(Grad(u), Grad(v))
          - Integral(f, v)
          + DirichletBC(u, g);

  // u.getSolution().projectOnBoundary(g);
  // u.getSolution().save("poisson.gf");

  Solver::CG solver;
  poisson.solve(solver);

  std::cout << poisson.getStiffnessOperator() << std::endl;
  std::cout << poisson.getMassVector() << std::endl;

  // std::cout << "weights: \n" << u.getSolution().getWeights().value() << std::endl;
  // std::cout << "guess: \n" << u.getSolution().getData() << std::endl;

  // Save solution
  u.getSolution().save("Poisson.gf");
  mesh.save("Poisson.mesh");

  return 0;
}
