/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

#include <Rodin/PETSc.h>
#include <Rodin/PETSc/Assembly/Sequential.h>

using namespace Rodin;
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);
  PetscPrintf(PETSC_COMM_WORLD,"Hello World\n");

  // Build a mesh
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Quadrilateral, { 32, 32 });

  ScalarFunction f(1.0);

  P1 vh(mesh);

  TrialFunction u(vh);
  TestFunction  v(vh);

  // PETSc::Matrix x;
  // BilinearForm bf(u, v, std::move(x));
  // bf = Integral(Grad(u), Grad(v));
  // bf.assemble();

  // PETSc::Vector b;
  // LinearForm lf(v, b);
  // lf = Integral(f, v);
  // lf.assemble();
  // lf.getVector();

  // // Define problem
  // Problem poisson(u, v);
  // poisson = Integral(Grad(u), Grad(v))
  //         - Integral(f, v)
  //         + DirichletBC(u, Zero());

  // // Solve
  // CG(poisson).solve();

  // // Save solution
  // u.getSolution().save("Poisson.gf");
  // mesh.save("Poisson.mesh");

  PetscFinalize();
  return 0;
}

