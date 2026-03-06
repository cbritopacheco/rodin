/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
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

  // Build a mesh
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Hexahedron, { 16, 16, 16 });
  mesh.scale(1.0 / 127.0);
  mesh.getConnectivity().compute(2, 3);

  RealFunction f = 1;

  P1 vh(mesh);

  {
    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TestFunction  v(vh);

    // Define problem
    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());
    std::cout << "Assembling..." << std::endl;
    poisson.assemble();

    std::cout << "Solving..." << std::endl;
    CG(poisson).solve();


    // Save solution
    u.getSolution().save("Poisson.gf");
    mesh.save("Poisson.mesh");
  }

  PetscFinalize();

  return 0;
}

