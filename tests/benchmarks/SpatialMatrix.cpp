/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Native matrix assembly versus nine independent scalar assemblies. */
#include <benchmark/benchmark.h>
#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"
using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
namespace
{
  template <size_t K>
  void matrixAssembly(benchmark::State& state)
  {
    const bool threeD = state.range(0) == 3;
    auto mesh = threeD ? LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3})
                       : LocalMesh::UniformGrid(Polytope::Type::Triangle, {8, 8});
    for (size_t d = 1; d <= mesh.getDimension(); ++d)
      for (size_t l = 0; l < d; ++l)
        mesh.getConnectivity().compute(d, l);
    auto run = [&](const auto& fes, size_t repetitions) {
      TrialFunction u(fes);
      TestFunction v(fes);
      BilinearForm form(u, v);
      if (state.range(1))
        form = Integral(Grad(u), Grad(v));
      else
        form = Integral(u, v);
      for (auto iteration : state)
      {
        for (size_t repeat = 0; repeat < repetitions; ++repeat)
        {
          form.assemble();
          benchmark::DoNotOptimize(form.getOperator().valuePtr());
        }
        benchmark::ClobberMemory();
      }
      state.counters["matrix_dofs"] = fes.getSize() * (repetitions == 9 ? 9 : 1);
    };
    if (state.range(2))
      run(H1(std::integral_constant<size_t, K>{}, mesh), 9);
    else
      run(H1(std::integral_constant<size_t, K>{}, mesh, 3, 3), 1);
  }
  // Arguments: dimension, stiffness (otherwise mass), nine scalar solves (otherwise native matrix).
  BENCHMARK_TEMPLATE(matrixAssembly, 1)
    ->ArgsProduct({{2, 3}, {0, 1}, {0, 1}})
    ->UseRealTime();
  BENCHMARK_TEMPLATE(matrixAssembly, 2)
    ->ArgsProduct({{2, 3}, {0, 1}, {0, 1}})
    ->UseRealTime();
  BENCHMARK_TEMPLATE(matrixAssembly, 3)
    ->ArgsProduct({{2, 3}, {0, 1}, {0, 1}})
    ->UseRealTime();
}
