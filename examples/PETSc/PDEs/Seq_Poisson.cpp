/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <chrono>

#include <Rodin/PETSc.h>

#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Alert.h>

using namespace Rodin;
using namespace Rodin::Math;
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int argc, char** argv)
{
  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);
  assert(ierr == PETSC_SUCCESS);

  size_t n = 16;

  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Hexahedron, { n, n, n });

  Alert::Info() << "Number of cells in the mesh: " << mesh.getCellCount() << Alert::Raise;
  Alert::Info() << "Number of vertices in the mesh: " << mesh.getVertexCount() << Alert::Raise;

  mesh.scale(1.0 / (n - 1));

  Alert::Info() << "Computing mesh connectivity..." << Alert::Raise;
  auto t0 = std::chrono::high_resolution_clock::now();

  const size_t cellDim = mesh.getDimension();
  mesh.getConnectivity().compute(cellDim - 1, cellDim);
  mesh.getConnectivity().compute(cellDim, cellDim);

  auto t1 = std::chrono::high_resolution_clock::now();

  Alert::Success()
    << "Computed mesh connectivity in "
    << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
    << " ms."
    << Alert::Raise;

  RealFunction f = 1;

  Alert::Info() << "Constructing finite element space..." << Alert::Raise;
  P1 vh(mesh);
  Alert::Success() << "Constructed finite element space." << Alert::Raise;

  {
    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());

    Alert::Info() << "Assembling linear system..." << Alert::Raise;

    t0 = std::chrono::high_resolution_clock::now();
    poisson.assemble();
    t1 = std::chrono::high_resolution_clock::now();

    Alert::Success()
      << "Assembled linear system in "
      << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
      << " ms."
      << Alert::Raise;

    Alert::Info() << "Solving linear system..." << Alert::Raise;

    t0 = std::chrono::high_resolution_clock::now();
    CG(poisson).solve();
    t1 = std::chrono::high_resolution_clock::now();

    Alert::Success()
      << "Solved linear system in "
      << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
      << " ms."
      << Alert::Raise;

    Alert::Info() << "Saving solution and mesh..." << Alert::Raise;

    u.getSolution().save("Poisson.gf", IO::FileFormat::MFEM);
    mesh.save("Poisson.mesh", IO::FileFormat::MFEM);

    Alert::Success() << "Saved solution and mesh to Poisson.gf and Poisson.mesh." << Alert::Raise;
  }

  (void) ierr;
  PetscFinalize();
  return 0;
}
