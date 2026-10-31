/*
 *          Copyright Carlos BRITO PACHECO 2026.
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
#include <Rodin/Kokkos/PETSc.h>

#ifdef KOKKOS_ENABLE_CUDA
using Backend = Rodin::Context::CUDA;
#else
using Backend = ::Kokkos::DefaultHostExecutionSpace;
#endif

// A global using directive on Rodin would make Kokkos ambiguous in the Kokkos headers.
using DeviceContext = Rodin::Context::Kokkos<Rodin::Context::Local, Backend>;

/**
 * @brief Elapsed milliseconds since a time point.
 * @param t0 Start time.
 * @return Milliseconds since @p t0.
 */
static auto since(std::chrono::high_resolution_clock::time_point t0)
{
  return std::chrono::duration_cast<std::chrono::milliseconds>(
    std::chrono::high_resolution_clock::now() - t0)
    .count();
}

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;

  ::Kokkos::ScopeGuard runtime(argc, argv);
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);
  {
    const size_t n = argc > 1 ? std::stoul(argv[1]) : 16;
    Mesh host = Mesh().UniformGrid(Polytope::Type::Hexahedron, {n, n, n});
    host.scale(1.0 / (n - 1));
    Alert::Info() << "Number of cells in the mesh: " << host.getCellCount()
                  << Alert::Raise;
    Alert::Info() << "Number of vertices in the mesh: " << host.getVertexCount()
                  << Alert::Raise;

    auto t0 = std::chrono::high_resolution_clock::now();
    const size_t cellDim = host.getDimension();
    host.getConnectivity().compute(cellDim - 1, cellDim);
    const Mesh<DeviceContext> mesh(DeviceContext{}, host);
    Alert::Success() << "Computed the connectivity and copied the mesh to the device in "
                     << since(t0) << " ms." << Alert::Raise;

    RealFunction f = 1;
    P1 vh(mesh);
    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TestFunction v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v)) - Integral(f, v) + DirichletBC(u, Zero());

    t0 = std::chrono::high_resolution_clock::now();
    poisson.assemble();
    Alert::Success() << "Assembled the linear system on the device in " << since(t0)
                     << " ms." << Alert::Raise;

    t0 = std::chrono::high_resolution_clock::now();
    Solver::CG(poisson).solve();
    Alert::Success() << "Solved the linear system in " << since(t0) << " ms."
                     << Alert::Raise;

    u.getSolution().save("Poisson.gf", IO::FileFormat::MFEM);
    host.save("Poisson.mesh", IO::FileFormat::MFEM);
    Alert::Success() << "Saved solution and mesh to Poisson.gf and Poisson.mesh."
                     << Alert::Raise;
  }
  PetscFinalize();
  return 0;
}
