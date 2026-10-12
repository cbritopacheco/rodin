/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Separates fixed-capacity spatial-vector construction and storage costs.
 */
#include <benchmark/benchmark.h>
#include <array>
#include <string>
#include <vector>
#include "Rodin/Math/SpatialVector.h"

namespace Rodin::Tests::Benchmarks
{
  bool spatialVectorStorageFailed = false;

  enum class SpatialVectorStorageOperation
  {
    ZeroConstruct,
    ZeroResize,
    CopyAppend,
    TemporaryFillAppend,
    DirectFillAppend,
    CopyAssign
  };

  /**
   * @brief Measures fixed-capacity values inside a reserved owning array.
   *
   * Each sweep creates or assigns a batch of spatial vectors. Array allocation
   * and source preparation are untimed; append scopes include clearing the
   * previous batch. TemporaryFillAppend constructs a zeroed local vector,
   * fills its active components and moves it into the array. DirectFillAppend
   * constructs the zeroed vector in the final array slot before filling it.
   * CopyAppend uses preconstructed immutable vectors; CopyAssign retains both
   * the size and capacity of the destination array throughout timing.
   * ZeroResize constructs dimension zero and then changes the logical size.
   * These are operation scopes, not isolated instruction or allocator timings.
   *
   * ## Exact guard
   *
   * Every component and logical size is checked exactly before and after the
   * timed sweep, including the unused fixed capacity. All constructor-zeroed
   * inactive components must remain zero. A perturbed active or inactive
   * component must be rejected. Input pointers and memory barriers prevent
   * invariant-output hoisting. Dimensions zero through three, both scalar
   * fields and three batch sizes separate logical work from fixed storage.
   *
   * @tparam Scalar Real or complex component type.
   */
  template <class Scalar>
  class SpatialVectorStorage
  {
    public:
    /**
     * @brief Measures one storage operation and checks the resulting batch.
     * @param state Logical dimension, batch size and timing iterations.
     * @tparam Operation Storage scope, fixed before the timed loop.
     */
    template <SpatialVectorStorageOperation Operation>
    static void run(benchmark::State& state)
    {
      using Vector = Math::SpatialVector<Scalar>;
      const auto dim = static_cast<std::uint8_t>(state.range(0));
      const size_t count = state.range(1);
      constexpr bool zero = Operation == SpatialVectorStorageOperation::ZeroConstruct ||
        Operation == SpatialVectorStorageOperation::ZeroResize;
      std::vector<std::array<Scalar, Vector::MaxSize>> samples(count);
      std::vector<Vector> source;
      source.reserve(count);
      for (size_t a = 0; a < count; ++a)
      {
        auto& point = source.emplace_back(dim);
        for (size_t j = 0; j < dim; ++j)
        {
          if constexpr (std::is_same_v<Scalar, Complex>)
            samples[a][j] = Scalar(Real(a + j + 1) / 7, -Real(a + 2 * j + 1) / 11);
          else
            samples[a][j] = Real(a + j + 1) / 7;
          point(j) = samples[a][j];
        }
      }
      std::vector<Vector> output;
      output.reserve(count);
      if constexpr (Operation == SpatialVectorStorageOperation::CopyAssign)
        output.resize(count);
      const auto evaluate = [&] {
        benchmark::DoNotOptimize(samples.data());
        benchmark::DoNotOptimize(source.data());
        benchmark::ClobberMemory();
        if constexpr (Operation == SpatialVectorStorageOperation::CopyAssign)
        {
          for (size_t a = 0; a < count; ++a)
            output[a] = source[a];
        }
        else
        {
          output.clear();
          for (size_t a = 0; a < count; ++a)
          {
            switch (Operation)
            {
              case SpatialVectorStorageOperation::ZeroConstruct:
                output.emplace_back(dim);
                break;
              case SpatialVectorStorageOperation::ZeroResize:
                output.emplace_back().resize(dim);
                break;
              case SpatialVectorStorageOperation::CopyAppend:
                output.push_back(source[a]);
                break;
              case SpatialVectorStorageOperation::TemporaryFillAppend:
              {
                Vector point(dim);
                for (size_t j = 0; j < dim; ++j)
                  point(j) = samples[a][j];
                output.push_back(std::move(point));
                break;
              }
              case SpatialVectorStorageOperation::DirectFillAppend:
              {
                auto& point = output.emplace_back(dim);
                for (size_t j = 0; j < dim; ++j)
                  point(j) = samples[a][j];
                break;
              }
              case SpatialVectorStorageOperation::CopyAssign:
                break;
            }
          }
        }
        benchmark::DoNotOptimize(output.data());
        benchmark::ClobberMemory();
      };
      const auto valid = [&] {
        if (output.size() != count)
          return false;
        for (size_t a = 0; a < count; ++a)
        {
          if (output[a].size() != dim)
            return false;
          for (size_t j = 0; j < Vector::MaxSize; ++j)
          {
            const Scalar expected = zero || j >= dim ? Scalar(0) : samples[a][j];
            if (output[a].getData()(j) != expected)
              return false;
          }
        }
        return true;
      };
      const auto fail = [&](const char* message) {
        spatialVectorStorageFailed = true;
        state.SkipWithError(message);
      };
      evaluate();
      if (!valid())
      {
        fail("Spatial-vector storage differs from exact size/component oracle");
        return;
      }
      output.front().getData()(0) += Scalar(1);
      if (valid())
      {
        fail("Storage oracle accepted a perturbed component");
        return;
      }
      evaluate();
      for (auto _ : state)
        evaluate();
      if (!valid())
        fail("Timed storage sweep changed size, active values or zeroed capacity");
      state.counters["dimension"] = dim;
      state.counters["vectors_per_sweep"] = count;
      state.counters["active_entries"] = count * dim;
      state.counters["capacity_entries"] = count * Vector::MaxSize;
      state.counters["storage_bytes"] = count * sizeof(Vector);
    }

    /**
     * @brief Registers all logical dimensions and three batch sizes.
     * @param scalar Scalar-field label.
     */
    static void registerCases(const char* scalar)
    {
      const auto registerOperation = [&]<SpatialVectorStorageOperation Operation>(
                                       const char* label) {
        const std::string name =
          std::string("SpatialVectorStorage/") + label + "/" + scalar;
        auto* registration = benchmark::RegisterBenchmark(
          name.c_str(), [](auto& state) { run<Operation>(state); });
        for (const size_t dim : {0u, 1u, 2u, 3u})
        {
          for (const size_t count : {32u, 256u, 2048u})
            registration->Args({static_cast<int64_t>(dim), static_cast<int64_t>(count)});
        }
        registration->UseRealTime();
      };
      registerOperation.template operator()<SpatialVectorStorageOperation::ZeroConstruct>(
        "ZeroConstruct");
      registerOperation.template operator()<SpatialVectorStorageOperation::ZeroResize>(
        "ZeroResize");
      registerOperation.template operator()<SpatialVectorStorageOperation::CopyAppend>(
        "CopyAppend");
      registerOperation
        .template operator()<SpatialVectorStorageOperation::TemporaryFillAppend>(
          "TemporaryFillAppend");
      registerOperation.template
      operator()<SpatialVectorStorageOperation::DirectFillAppend>("DirectFillAppend");
      registerOperation.template operator()<SpatialVectorStorageOperation::CopyAssign>(
        "CopyAssign");
    }
  };
}

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Tests::Benchmarks;
  SpatialVectorStorage<Real>::registerCases("Real");
  SpatialVectorStorage<Complex>::registerCases("Complex");
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  benchmark::AddCustomContext(
    "scope", "Reserved spatial-value storage; no allocation or assembly");
  const size_t cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || spatialVectorStorageFailed ? 1 : 0;
}
