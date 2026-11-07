/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file H1_Assembly.cpp
 * @brief Times the device assembly of H1 forms on a grid.
 *
 * The first argument is the number of grid nodes per side, the second the
 * number of assemblies of each form, the third the degree, 2 or 3. The scalar
 * form is the stiffness of an H1 space, the vector form the Jacobian product
 * of its vector counterpart, the hyperelastic form the Neo-Hookean tangent
 * of the vector space at a small displacement. The first assembly of each form also builds the sparsity, so
 * the later ones measure the integration kernel and the matrix update. A
 * fourth argument `tetrahedra` meshes the grid with tetrahedra instead of
 * hexahedra.
 * @code
 *   KokkosH1_Assembly 30 3 2
 *   KokkosH1_Assembly 30 3 2 tetrahedra
 * @endcode
 */
#include <chrono>
#include <cmath>
#include <string>

#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Alert.h>
#include <Rodin/Solid.h>
#include <Rodin/Kokkos.h>

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

/**
 * @brief Assembles a problem several times and reports each duration.
 * @param name Name of the form.
 * @param problem Problem to assemble.
 * @param repetitions Number of assemblies.
 */
template <class Problem>
static void time(const char* name, Problem& problem, size_t repetitions)
{
  for (size_t i = 0; i < repetitions; ++i)
  {
    const auto t0 = std::chrono::high_resolution_clock::now();
    problem.assemble();
    ::Kokkos::fence();
    Rodin::Alert::Success() << name << " assembly " << i + 1 << ": " << since(t0)
                            << " ms." << Rodin::Alert::Raise;
  }
}

/**
 * @brief Times the scalar and vector forms of degree K.
 * @param mesh Device mesh.
 * @param repetitions Number of assemblies of each form.
 */
template <size_t K, class Mesh>
static void run(const Mesh& mesh, size_t repetitions)
{
  using namespace Rodin;
  using namespace Rodin::Variational;
  constexpr std::integral_constant<size_t, K> degree{};
  {
    H1 vh(degree, mesh);
    Alert::Info() << "Scalar degrees of freedom: " << vh.getSize() << Alert::Raise;
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem<DeviceSystem, decltype(u), decltype(v)> problem(u, v);
    problem = Integral(Grad(u), Grad(v));
    time("Scalar", problem, repetitions);
  }
  {
    H1 vh(degree, mesh, mesh.getSpaceDimension());
    Alert::Info() << "Vector degrees of freedom: " << vh.getSize() << Alert::Raise;
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem<DeviceSystem, decltype(u), decltype(v)> problem(u, v);
    problem = Integral(Jacobian(u), Jacobian(v));
    time("Vector", problem, repetitions);
  }
  {
    H1 vh(degree, mesh, mesh.getSpaceDimension());
    GridFunction w(vh);
    auto& data = w.getData();
    for (Index i = 0; i < static_cast<Index>(data.size()); ++i)
      data.coeffRef(i) = Real(0.005) * std::sin(Real(0.7) * i);
    TrialFunction u(vh);
    TestFunction v(vh);
    const auto work = Solid::InternalVirtualWork(Solid::NeoHookean(1.2, 0.8), w);
    Problem<DeviceSystem, decltype(u), decltype(v)> problem(u, v);
    problem = work.Tangent(u, v);
    time("Hyperelastic", problem, repetitions);
  }
}

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Geometry;

  ::Kokkos::ScopeGuard runtime(argc, argv);
  {
    const size_t n = argc > 1 ? std::stoul(argv[1]) : 16;
    const size_t repetitions = argc > 2 ? std::stoul(argv[2]) : 3;
    const size_t degree = argc > 3 ? std::stoul(argv[3]) : 2;
    const bool tetrahedra = argc > 4 && std::string(argv[4]) == "tetrahedra";
    Mesh host = Mesh().UniformGrid(
      tetrahedra ? Polytope::Type::Tetrahedron : Polytope::Type::Hexahedron, {n, n, n});
    host.scale(1.0 / (n - 1));
    Alert::Info() << "Number of cells in the mesh: " << host.getCellCount()
                  << Alert::Raise;

    const size_t cellDim = host.getDimension();
    host.getConnectivity().compute(cellDim - 1, cellDim);
    const Mesh<DeviceContext> mesh(DeviceContext{}, host);
    if (degree == 3)
      run<3>(mesh, repetitions);
    else
      run<2>(mesh, repetitions);
  }
  return 0;
}
