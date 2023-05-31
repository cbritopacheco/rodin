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
  mesh.load(meshFile);

  // Functions
  P1 vh(mesh);

  TrialFunction u(vh);
  TestFunction  v(vh);

  // Define problem
  ScalarFunction f = 1.0;
  ScalarFunction g = 0.0;

 // Solver::CG solver;
 // solver.printIterations(true);

 // Problem poisson(u, v);
 // poisson = Integral(Grad(u), Grad(v))
 //         - Integral(f, v)
 //         + DirichletBC(u, g).on(Gamma);

 // poisson.solve(solver);

 // // Save solution
 // u.getSolution().getHandle().Save("u.gf");
 // Omega.save("Omega.mesh");

  return 0;
}
