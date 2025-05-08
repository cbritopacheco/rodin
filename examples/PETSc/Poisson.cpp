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
using namespace Rodin::Math;
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


  // BilinearForm bf(u, v, x);
  // bf = Integral(Grad(u), Grad(v));
  // bf.assemble();
  // bf.getOperator();

  // LinearForm lf(v, b);
  // lf = Integral(f, v);
  // lf.assemble();

  // Mat a;
  // MatCreate(PETSC_COMM_SELF, &a);

  // Vec x;
  // VecCreate(PETSC_COMM_SELF, &x);

  // Vec b;
  // VecCreate(PETSC_COMM_SELF, &b);

  // LinearSystem axb(a, x, b);

  // // Define problem
  // Problem poisson(u, v, axb);
  // poisson = Integral(Grad(u), Grad(v))
  //         - Integral(f, v)
  //         + DirichletBC(u, Zero());

  // // // Solve
  // // CG(poisson).solve();

  // // // Save solution
  // // u.getSolution().save("Poisson.gf");
  // // mesh.save("Poisson.mesh");

  // MatDestroy(&a);
  // VecDestroy(&x);
  // VecDestroy(&b);

  PetscFinalize();
  return 0;
}

