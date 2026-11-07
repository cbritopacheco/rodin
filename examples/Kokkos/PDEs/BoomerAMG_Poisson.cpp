/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file BoomerAMG_Poisson.cpp
 * @brief Poisson problem on a hexahedral grid, assembled on the device and
 * solved by the conjugate gradient of hypre preconditioned by BoomerAMG,
 * without PETSc.
 *
 * The matrix stays on the device from the assembly to the solve. The
 * arguments are the number of grid nodes per side, the strength threshold,
 * the levels of aggressive coarsening and the maximum number of elements per
 * row of the interpolation. The hierarchy is built once and reused by the
 * solve that follows a second assembly. Settings suited to three dimensional problems:
 * @code
 *   KokkosBoomerAMG_Poisson 100 0.7 2 4
 * @endcode
 */
#include <chrono>

#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Alert.h>
#include <Rodin/Kokkos.h>
#include <Rodin/Kokkos/Solver/BoomerAMG.h>

#if defined(KOKKOS_ENABLE_CUDA)
using Backend = Rodin::Context::CUDA;

#elif defined(KOKKOS_ENABLE_HIP)
using Backend = Rodin::Context::HIP;
#elif defined(KOKKOS_ENABLE_SYCL)
using Backend = Rodin::Context::SYCL;
#else
using Backend = ::Kokkos::DefaultHostExecutionSpace;
#endif

// A global using directive on Rodin would make Kokkos ambiguous in the Kokkos headers.
using DeviceContext = Rodin::Context::Kokkos<Rodin::Context::Local, Backend>;

/// @brief System whose matrix stays on the device.
using DeviceSystem =
  Rodin::Math::LinearSystem<Rodin::Kokkos::SparseMatrix<Rodin::Real, Backend>,
    Rodin::Math::Vector<Rodin::Real>>;

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
  {
    const size_t n = argc > 1 ? std::stoul(argv[1]) : 16;
    const Real threshold = argc > 2 ? std::stod(argv[2]) : 0.7;
    const size_t aggressive = argc > 3 ? std::stoul(argv[3]) : 2;
    const size_t interpolation = argc > 4 ? std::stoul(argv[4]) : 4;
    Mesh host = Mesh().UniformGrid(Polytope::Type::Hexahedron, {n, n, n});
    host.scale(1.0 / (n - 1));
    Alert::Info() << "Number of vertices in the mesh: " << host.getVertexCount()
                  << Alert::Raise;

    const size_t cellDim = host.getDimension();
    host.getConnectivity().compute(cellDim - 1, cellDim);
    const Mesh<DeviceContext> mesh(DeviceContext{}, host);

    RealFunction f = 1;
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);

    Problem<DeviceSystem, decltype(u), decltype(v)> poisson(u, v);
    poisson = Integral(Grad(u), Grad(v)) - Integral(f, v) + DirichletBC(u, Zero());

    auto t0 = std::chrono::high_resolution_clock::now();
    poisson.assemble();
    Alert::Success() << "Assembled the linear system on the device in " << since(t0)
                     << " ms." << Alert::Raise;

    Solver::BoomerAMG solver(poisson);
    solver.setStrongThreshold(threshold)
      .setAggressiveLevels(aggressive)
      .setMaxInterpolationElements(interpolation);
    t0 = std::chrono::high_resolution_clock::now();
    solver.setup(poisson.getLinearSystem());
    Alert::Success() << "Built the hierarchy in " << since(t0) << " ms." << Alert::Raise;

    // The second solve, after a new assembly, reuses the hierarchy, as a time
    // step would.
    for (size_t step = 1; step <= 2; ++step)
    {
      if (step > 1)
        poisson.assemble();
      t0 = std::chrono::high_resolution_clock::now();
      solver.solve(poisson.getLinearSystem());
      Alert::Success() << "Solved linear system " << step << " in " << since(t0)
                       << " ms, " << solver.getIterationNumber()
                       << " iterations, relative residual " << solver.getError() << "."
                       << Alert::Raise;
    }
  }
  // The solver initialized hypre and MPI, whose finalization is left to the
  // application.
  HYPRE_Finalize();
  MPI_Finalize();
  return 0;
}
