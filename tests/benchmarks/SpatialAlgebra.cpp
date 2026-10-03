/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Fixed-capacity spatial arithmetic and contraction timings. */
#include <array>
#include <limits>
#include <string>
#include <type_traits>
#include <benchmark/benchmark.h>
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/SpatialTensor.h"

namespace Rodin::Tests::Benchmarks
{
  /**
   * @brief Times spatial value operations with runtime extents and escaping results.
   *
   * Operand initialization and independent contraction checks precede timing.
   * Compiler barriers expose operands and results on every iteration so neither
   * constant folding nor dead-store removal can eliminate the measured operation.
   */
  template <class Scalar, size_t Rank>
  class SpatialAlgebraBenchmark
  {
    public:
      enum class Operation
      {
        Construct,
        Copy,
        Add,
        Scale,
        Dot
      };

      static auto make(const benchmark::State& state)
      {
        if constexpr (Rank == 1)
          return Math::SpatialVector<Scalar>(state.range(0));
        else if constexpr (Rank == 2)
          return Math::SpatialMatrix<Scalar>(state.range(0), state.range(1));
        else
        {
          std::array<size_t, Rank> extents;
          for (size_t axis = 0; axis < Rank; ++axis)
            extents[axis] = state.range(axis);
          return Math::SpatialTensor<Scalar, Rank>(extents);
        }
      }

      static Scalar coefficient(size_t i, size_t offset)
      {
        // Bounded, nonuniform values; complex inputs have nonzero imaginary parts.
        const Real real = Real(1 + (i + offset) % 7) / 8;
        if constexpr (std::is_same_v<Scalar, Complex>)
          return Scalar(real, Real(1 + (2 * i + offset) % 5) / 8);
        else
          return real;
      }

      template <class Value>
      static size_t count(const Value& value)
      {
        if constexpr (requires {
                        value.rows();
                        value.cols();
                      })
          return value.rows() * value.cols();
        else
          return value.size();
      }

      template <class Value>
      static Scalar& entry(Value& value, size_t i)
      {
        if constexpr (Rank == 2)
          return value(i / value.cols(), i % value.cols());
        else
          return value[i];
      }

      template <Operation Op>
      static void arithmetic(benchmark::State& state)
      {
        auto lhs = make(state);
        auto rhs = make(state);
        for (size_t i = 0; i < count(lhs); ++i)
        {
          entry(lhs, i) = coefficient(i, 0);
          entry(rhs, i) = coefficient(i, 1);
        }
        const Scalar factor = coefficient(0, 2);
        for (auto iteration : state)
        {
          benchmark::DoNotOptimize(&lhs);
          benchmark::DoNotOptimize(&rhs);
          if constexpr (Op == Operation::Construct)
          {
            auto result = make(state);
            benchmark::DoNotOptimize(&result);
            benchmark::ClobberMemory();
          }
          else if constexpr (Op == Operation::Copy)
          {
            auto result = lhs;
            benchmark::DoNotOptimize(&result);
            benchmark::ClobberMemory();
          }
          else if constexpr (Op == Operation::Add)
          {
            auto result = lhs + rhs;
            benchmark::DoNotOptimize(&result);
            benchmark::ClobberMemory();
          }
          else if constexpr (Op == Operation::Scale)
          {
            auto result = lhs * factor;
            benchmark::DoNotOptimize(&result);
            benchmark::ClobberMemory();
          }
          else
          {
            auto result = lhs.dot(rhs);
            benchmark::DoNotOptimize(result);
            benchmark::ClobberMemory();
          }
        }
        state.counters["entries"] = count(lhs);
        state.counters["object_bytes"] = sizeof(lhs);
        state.SetItemsProcessed(state.iterations() * count(lhs));
      }

      static void contraction(benchmark::State& state)
        requires(Rank >= 2)
      {
        auto lhs = make(state);
        for (size_t i = 0; i < count(lhs); ++i)
          entry(lhs, i) = coefficient(i, 0);
        auto rhs = [&] {
          if constexpr (Rank == 2)
            return Math::SpatialVector<Scalar>(lhs.cols());
          else if constexpr (Rank == 3)
            return Math::SpatialVector<Scalar>(lhs.getDimension(2));
          else
            return Math::SpatialMatrix<Scalar>(lhs.getDimension(2), lhs.getDimension(3));
        }();
        for (size_t i = 0; i < count(rhs); ++i)
        {
          if constexpr (Rank == 4)
            rhs(i / rhs.cols(), i % rhs.cols()) = coefficient(i, 1);
          else
            rhs[i] = coefficient(i, 1);
        }
        // Independently flatten only the reference calculation, outside timing.
        const auto checked = lhs * rhs;
        const size_t contracted = count(rhs);
        for (size_t i = 0; i < count(checked); ++i)
        {
          Scalar expected = 0;
          for (size_t j = 0; j < contracted; ++j)
          {
            if constexpr (Rank == 4)
              expected +=
                entry(lhs, i * contracted + j) * rhs(j / rhs.cols(), j % rhs.cols());
            else
              expected += entry(lhs, i * contracted + j) * rhs[j];
          }
          const Scalar actual = [&] {
            if constexpr (Rank == 2)
              return checked[i];
            else
              return checked(i / checked.cols(), i % checked.cols());
          }();
          // Sum of at most nine products; allow rounding across accumulation orders.
          constexpr Real tolerance = 32 * std::numeric_limits<Real>::epsilon();
          if (std::abs(actual - expected) > tolerance * (1 + std::abs(expected)))
          {
            state.SkipWithError("Spatial contraction disagrees with reference sum.");
            return;
          }
        }
        for (auto iteration : state)
        {
          benchmark::DoNotOptimize(&lhs);
          benchmark::DoNotOptimize(&rhs);
          auto result = lhs * rhs;
          benchmark::DoNotOptimize(&result);
          benchmark::ClobberMemory();
        }
        state.counters["entries"] = count(lhs);
        state.counters["object_bytes"] = sizeof(lhs);
        state.SetItemsProcessed(state.iterations() * count(lhs));
      }

      static void matrixProduct(benchmark::State& state)
        requires(Rank == 2)
      {
        Math::SpatialMatrix<Scalar> lhs(state.range(0), state.range(1));
        Math::SpatialMatrix<Scalar> rhs(state.range(1), state.range(2));
        for (size_t i = 0; i < count(lhs); ++i)
          entry(lhs, i) = coefficient(i, 0);
        for (size_t i = 0; i < count(rhs); ++i)
          entry(rhs, i) = coefficient(i, 1);
        const auto checked = lhs * rhs;
        for (size_t i = 0; i < checked.rows(); ++i)
          for (size_t j = 0; j < checked.cols(); ++j)
          {
            Scalar expected = 0;
            for (size_t k = 0; k < lhs.cols(); ++k)
              expected += lhs(i, k) * rhs(k, j);
            constexpr Real tolerance = 32 * std::numeric_limits<Real>::epsilon();
            if (std::abs(checked(i, j) - expected) > tolerance * (1 + std::abs(expected)))
            {
              state.SkipWithError("Spatial matrix product disagrees with reference sum.");
              return;
            }
          }
        for (auto iteration : state)
        {
          benchmark::DoNotOptimize(&lhs);
          benchmark::DoNotOptimize(&rhs);
          auto result = lhs * rhs;
          benchmark::DoNotOptimize(&result);
          benchmark::ClobberMemory();
        }
        state.SetItemsProcessed(state.iterations() * count(checked) * lhs.cols());
      }

      static void registerCases(const std::string& name)
      {
        using Op = Operation;
        auto shape = [](auto* benchmark) {
          if constexpr (Rank == 1)
            benchmark->Arg(0)->Arg(1)->Arg(2)->Arg(3);
          else if constexpr (Rank == 2)
            benchmark->Args({0, 0})->ArgsProduct({{1, 2, 3}, {1, 2, 3}});
          else if constexpr (Rank == 3)
            benchmark->Args({0, 0, 0})
              ->Args({1, 1, 1})
              ->Args({2, 2, 2})
              ->Args({3, 3, 3})
              ->Args({2, 3, 2})
              ->Args({3, 2, 3});
          else
            benchmark->Args({0, 0, 0, 0})
              ->Args({1, 1, 1, 1})
              ->Args({2, 2, 2, 2})
              ->Args({3, 3, 3, 3})
              ->Args({2, 3, 3, 2})
              ->Args({3, 2, 2, 3});
        };
        shape(benchmark::RegisterBenchmark(
          (name + "/Construct").c_str(), &arithmetic<Op::Construct>));
        shape(
          benchmark::RegisterBenchmark((name + "/Copy").c_str(), &arithmetic<Op::Copy>));
        shape(
          benchmark::RegisterBenchmark((name + "/Add").c_str(), &arithmetic<Op::Add>));
        shape(benchmark::RegisterBenchmark(
          (name + "/Scale").c_str(), &arithmetic<Op::Scale>));
        shape(
          benchmark::RegisterBenchmark((name + "/Dot").c_str(), &arithmetic<Op::Dot>));
        if constexpr (Rank >= 2)
          shape(benchmark::RegisterBenchmark((name + "/Contract").c_str(), &contraction));
        if constexpr (Rank == 2)
          benchmark::RegisterBenchmark((name + "/MatMat").c_str(), &matrixProduct)
            ->ArgsProduct({{1, 2, 3}, {1, 2, 3}, {1, 2, 3}});
      }
  };

  const bool RegisteredSpatialAlgebra = [] {
    SpatialAlgebraBenchmark<Real, 1>::registerCases("SpatialVector/Real");
    SpatialAlgebraBenchmark<Complex, 1>::registerCases("SpatialVector/Complex");
    SpatialAlgebraBenchmark<Real, 2>::registerCases("SpatialMatrix/Real");
    SpatialAlgebraBenchmark<Complex, 2>::registerCases("SpatialMatrix/Complex");
    SpatialAlgebraBenchmark<Real, 3>::registerCases("SpatialTensor3/Real");
    SpatialAlgebraBenchmark<Complex, 3>::registerCases("SpatialTensor3/Complex");
    SpatialAlgebraBenchmark<Real, 4>::registerCases("SpatialTensor4/Real");
    SpatialAlgebraBenchmark<Complex, 4>::registerCases("SpatialTensor4/Complex");
    return true;
  }();
}
