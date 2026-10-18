/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Scalar and one-component contraction representation measurements.
 */
#include <algorithm>
#include <benchmark/benchmark.h>
#include <complex>
#include <string>
#include <type_traits>
#include <vector>
#include "Rodin/Math/SpatialVector.h"

namespace Rodin::Tests::Benchmarks
{
  /// @brief Records untimed numerical qualification failures.
  inline bool scalarDotFailed = false;

  /// @brief Public contraction representation, including explicit zero seeds.
  enum class ScalarDotPrimitive
  {
    Scalar,
    ScalarSeeded,
    FixedOne,
    FixedOneSeeded,
    SpatialOne
  };

  /**
   * @brief Measures one-component pair contractions without finite-element setup.
   *
   * ## Mathematical foundation
   *
   * A complete batch computes @f$D_{ba}=z_b\overline{y_a}@f$ or replaces
   * a zeroed output by @f$D_{ba}\mathrel{+}=wz_b\overline{y_a}@f$.
   * Seeded scopes additionally form @f$0+z_b\overline{y_a}@f$ before
   * storing or weighting it. These represent scalar arithmetic, fixed-size
   * Eigen views and the public one-component spatial-vector operation. No
   * mathematical equality of their general floating-point algorithms is assumed.
   *
   * ## Architecture
   *
   * Raw scalar arrays are shared by the scalar and fixed-view protocols.
   * Spatial protocols instead retain arrays of one-component spatial values;
   * their object footprint is reported separately. Dyadic inputs give exact
   * finite reference values, including genuinely complex right operands and
   * real-valued right operands representing reference basis values. A perturbed
   * output must fail qualification. This does not certify arbitrary complex,
   * nonfinite or rounding-mode behavior across representations.
   *
   * Setup and qualification are untimed. Each iteration includes pair traversal
   * and output stores; weighted scopes also include output zeroing. Barriers
   * surround a batch rather than each contraction. No geometry, quadrature,
   * basis scaling, matrix scatter or solving is measured.
   *
   * @tparam Scalar Real or complex scalar type.
   * @tparam Primitive Public representation and zero-seed protocol.
   * @tparam Weighted Whether a weighted accumulation is measured.
   * @tparam RealTest Whether right inputs have zero imaginary part.
   */
  template <class Scalar, ScalarDotPrimitive Primitive, bool Weighted, bool RealTest>
  class ScalarDotBatch
  {
    public:
      /**
       * @brief Measures one complete pair batch.
       * @param state Operand count and iteration protocol.
       */
      static void run(benchmark::State& state)
      {
        constexpr bool Spatial = Primitive == ScalarDotPrimitive::SpatialOne;
        constexpr bool Fixed = Primitive == ScalarDotPrimitive::FixedOne ||
          Primitive == ScalarDotPrimitive::FixedOneSeeded;
        constexpr bool Seeded = Primitive == ScalarDotPrimitive::ScalarSeeded ||
          Primitive == ScalarDotPrimitive::FixedOneSeeded;
        using Operand = std::conditional_t<Spatial, Math::SpatialVector<Scalar>, Scalar>;
        using FixedVector = Eigen::Matrix<Scalar, 1, 1>;
        const size_t count = state.range(0);
        std::vector<Operand> lhs(count), rhs(count);
        for (size_t i = 0; i < count; ++i)
        {
          const Scalar left = [&]() {
            if constexpr (std::is_same_v<Scalar, Complex>)
              return Scalar(Real(i + 1) / 8, Real(i + 2) / 16);
            else
              return Scalar(Real(i + 1) / 8);
          }();
          const Scalar right = [&]() {
            if constexpr (std::is_same_v<Scalar, Complex> && !RealTest)
              return Scalar(Real(i + 1) / 16, -Real(i + 3) / 32);
            else
              return Scalar(Real(i + 1) / 16);
          }();
          if constexpr (Spatial)
          {
            lhs[i].resize(1);
            rhs[i].resize(1);
            lhs[i][0] = left;
            rhs[i][0] = right;
          }
          else
          {
            lhs[i] = left;
            rhs[i] = right;
          }
        }
        const auto scalar = [](const Operand& value) -> Scalar {
          if constexpr (Spatial)
            return value[0];
          else
            return value;
        };
        Scalar weight = Scalar(Real(3) / 4);
        std::vector<Scalar> expected(count * count), actual(count * count);
        for (size_t b = 0; b < count; ++b)
        {
          for (size_t a = 0; a < count; ++a)
          {
            const Scalar value = [&]() {
              if constexpr (std::is_same_v<Scalar, Complex>)
                return scalar(lhs[b]) * std::conj(scalar(rhs[a]));
              else
                return scalar(lhs[b]) * scalar(rhs[a]);
            }();
            expected[b * count + a] = Weighted ? Scalar(0) + weight * value : value;
          }
        }
        const auto compute = [&]() {
          if constexpr (Weighted)
            std::fill(actual.begin(), actual.end(), Scalar(0));
          for (size_t b = 0; b < count; ++b)
          {
            for (size_t a = 0; a < count; ++a)
            {
              const Scalar product = [&]() {
                if constexpr (Fixed)
                  return Math::dot(Eigen::Map<const FixedVector>(&lhs[b]),
                    Eigen::Map<const FixedVector>(&rhs[a]));
                else
                  return Math::dot(lhs[b], rhs[a]);
              }();
              const Scalar contraction = [&]() {
                if constexpr (Seeded)
                {
                  Scalar value = 0;
                  value += product;
                  return value;
                }
                else
                  return product;
              }();
              if constexpr (Weighted)
                actual[b * count + a] += weight * contraction;
              else
                actual[b * count + a] = contraction;
            }
          }
        };
        compute();
        if (actual != expected)
        {
          scalarDotFailed = true;
          state.SkipWithError(
            "One-component contraction differs from the dyadic reference");
          return;
        }
        actual[0] += Scalar(1);
        if (actual == expected)
        {
          scalarDotFailed = true;
          state.SkipWithError("One-component contraction accepted a perturbed entry");
          return;
        }
        compute();
        for (auto _ : state)
        {
          benchmark::DoNotOptimize(lhs.data());
          benchmark::DoNotOptimize(rhs.data());
          if constexpr (Weighted)
            benchmark::DoNotOptimize(weight);
          compute();
          benchmark::DoNotOptimize(actual.data());
          benchmark::ClobberMemory();
        }
        if (actual != expected)
        {
          scalarDotFailed = true;
          state.SkipWithError("Timed contraction changed a reference entry");
        }
        state.counters["components"] = 1;
        state.counters["operands_per_side"] = count;
        state.counters["pairs"] = count * count;
        state.counters["scalar_products"] = count * count;
        state.counters["input_object_bytes"] = 2 * count * sizeof(Operand);
        state.counters["output_bytes"] = count * count * sizeof(Scalar);
        state.counters["zeroed_entries"] = Weighted ? count * count : 0;
        state.counters["weighted"] = Weighted;
        state.counters["initial_zero"] = Seeded;
        state.counters["real_test_values"] = RealTest;
      }

      /**
       * @brief Registers three basis-sized operand counts.
       * @param label Scalar and right-input family name.
       */
      static void registerCases(const char* label)
      {
        const std::string name = std::string("ScalarDot/") +
          (Primitive == ScalarDotPrimitive::Scalar              ? "Scalar/"
              : Primitive == ScalarDotPrimitive::ScalarSeeded   ? "ScalarSeeded/"
              : Primitive == ScalarDotPrimitive::FixedOne       ? "FixedOne/"
              : Primitive == ScalarDotPrimitive::FixedOneSeeded ? "FixedOneSeeded/"
                                                                : "SpatialOne/") +
          label + (Weighted ? "/Weighted" : "/Dot");
        benchmark::RegisterBenchmark(name.c_str(), [](auto& state) { run(state); })
          ->Arg(8)
          ->Arg(27)
          ->Arg(343)
          ->UseRealTime();
      }
  };
}

/**
 * @brief Runs opt-in scalar and one-component contraction measurements.
 * @param argc Number of command-line arguments.
 * @param argv Benchmark selection and measurement arguments.
 * @returns Zero when selected measurements and numerical guards succeed.
 */
int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Tests::Benchmarks;
  const auto registerFamily = []<class Scalar, bool RealTest>(const char* label) {
    const auto primitive = [&]<ScalarDotPrimitive Kind>() {
      ScalarDotBatch<Scalar, Kind, false, RealTest>::registerCases(label);
      ScalarDotBatch<Scalar, Kind, true, RealTest>::registerCases(label);
    };
    primitive.template operator()<ScalarDotPrimitive::Scalar>();
    primitive.template operator()<ScalarDotPrimitive::ScalarSeeded>();
    primitive.template operator()<ScalarDotPrimitive::FixedOne>();
    primitive.template operator()<ScalarDotPrimitive::FixedOneSeeded>();
    primitive.template operator()<ScalarDotPrimitive::SpatialOne>();
  };
  registerFamily.template operator()<Real, true>("Real");
  registerFamily.template operator()<Complex, true>("ComplexRealTest");
  registerFamily.template operator()<Complex, false>("ComplexGeneral");
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  benchmark::AddCustomContext(
    "scope", "One-component pair contractions; no geometry or assembly");
  const size_t cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || scalarDotFailed ? 1 : 0;
}
