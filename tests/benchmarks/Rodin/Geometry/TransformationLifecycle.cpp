/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Measures geometric transformation construction and value lifetimes.
 */
#include <benchmark/benchmark.h>
#include "Rodin/Geometry/IdentityTransformation.h"

namespace Rodin::Tests::Benchmarks
{
  bool transformationLifecycleFailed = false;

  /**
   * @brief Measures identity-map construction without geometric evaluation.
   *
   * Construct measures one map. ConstructCopy and ConstructMove each construct
   * a source and a destination; source construction is deliberately included,
   * rather than hidden behind per-iteration timing pauses. Both addresses
   * escape through mutable pointer barriers. No heap allocation, mesh, formula,
   * assembly or solver is included. Dimensions and the identity Jacobian are
   * qualified outside timing, including the empty zero-dimensional map.
   *
   * @tparam Operation Zero for construction, one for construction/copy, two for
   * construction/move.
   * @param state Runtime map dimension and repetition protocol.
   */
  template <size_t Operation>
  void transformationLifecycle(benchmark::State& state)
  {
    const size_t dimension = state.range(0);
    const auto valid = [&]() {
      Geometry::IdentityTransformation map(dimension);
      Math::SpatialPoint reference(dimension);
      reference.setZero();
      Math::SpatialMatrix<Real> jacobian;
      map.jacobian(jacobian, reference);
      if (map.getReferenceDimension() != dimension ||
        map.getPhysicalDimension() != dimension || jacobian.rows() != dimension ||
        jacobian.cols() != dimension)
        return false;
      for (size_t i = 0; i < dimension; ++i)
      {
        for (size_t j = 0; j < dimension; ++j)
        {
          if (jacobian(i, j) != Real(i == j))
            return false;
        }
      }
      return true;
    };
    if (!valid())
    {
      transformationLifecycleFailed = true;
      state.SkipWithError("Identity-map dimension/Jacobian qualification failed");
      return;
    }
    for (auto _ : state)
    {
      Geometry::IdentityTransformation source(dimension);
      auto* address = &source;
      benchmark::DoNotOptimize(address);
      if constexpr (Operation != 0)
      {
        Geometry::IdentityTransformation destination = [&]() {
          if constexpr (Operation == 1)
            return Geometry::IdentityTransformation(source);
          else
            return Geometry::IdentityTransformation(std::move(source));
        }();
        auto* destinationAddress = &destination;
        benchmark::DoNotOptimize(destinationAddress);
      }
      benchmark::ClobberMemory();
    }
    if (!valid())
    {
      transformationLifecycleFailed = true;
      state.SkipWithError("Identity-map qualification changed after timing");
    }
    state.counters["dimension"] = dimension;
    state.counters["constructions_per_iteration"] = Operation == 0 ? 1 : 2;
    state.counters["object_bytes"] = sizeof(Geometry::IdentityTransformation);
  }
}

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  benchmark::RegisterBenchmark(
    "TransformationLifecycle/Identity/Construct", &transformationLifecycle<0>)
    ->DenseRange(0, 3)
    ->UseRealTime();
  benchmark::RegisterBenchmark(
    "TransformationLifecycle/Identity/ConstructCopy", &transformationLifecycle<1>)
    ->DenseRange(0, 3)
    ->UseRealTime();
  benchmark::RegisterBenchmark(
    "TransformationLifecycle/Identity/ConstructMove", &transformationLifecycle<2>)
    ->DenseRange(0, 3)
    ->UseRealTime();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  const auto cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || transformationLifecycleFailed ? 1 : 0;
}
