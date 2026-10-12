/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Paired public spatial-dot and original dynamic-reduction measurements.
 */
#include <benchmark/benchmark.h>
#include "Rodin/Math/SpatialVector.h"

namespace Rodin::Tests::Benchmarks
{
  bool spatialDotFailed = false;

  /**
   * @brief Measures the scalar contraction used inside diffusion kernels.
   *
   * The public operation computes @f$\sum_i u_i\overline{v_i}@f$ on the
   * active entries. The reference explicitly calls the original dynamic Eigen
   * reduction. Both operate on the same vectors; construction and exact result
   * guards are outside timing. This does not measure the differently defined
   * unconjugated SpatialVector member dot, mapped quadrature or full assembly.
   *
   * @tparam LHS Left operand's scalar type.
   * @tparam RHS Right operand's scalar type.
   * @tparam Original Whether to measure the original reduction.
   * @param state Logical dimension and timed iterations.
   */
  template <class LHS, class RHS, bool Original>
  void spatialDot(benchmark::State& state)
  {
    Math::SpatialVector<LHS> lhs(state.range(0));
    Math::SpatialVector<RHS> rhs(state.range(0));
    for (size_t i = 0; i < lhs.size(); ++i)
    {
      if constexpr (std::is_same_v<LHS, Complex>)
        lhs(i) = LHS(Real(i + 1) / 7, Real(i + 2) / 11);
      else
        lhs(i) = Real(i + 1) / 7;
      if constexpr (std::is_same_v<RHS, Complex>)
        rhs(i) = RHS(Real(i + 2) / 13, -Real(i + 1) / 17);
      else
        rhs(i) = Real(i + 2) / 13;
    }
    const auto reference = [&]() {
      const auto size = static_cast<Eigen::Index>(lhs.size());
      return Math::dot(lhs.getData().head(size), rhs.getData().head(size));
    };
    const auto valid = [&]() { return Math::dot(lhs, rhs) == reference(); };
    if (!valid())
    {
      spatialDotFailed = true;
      state.SkipWithError("Spatial dot differs from original reduction");
      return;
    }
    for (auto _ : state)
    {
      benchmark::DoNotOptimize(&lhs);
      benchmark::DoNotOptimize(&rhs);
      auto result = [&]() {
        if constexpr (Original)
          return reference();
        else
          return Math::dot(lhs, rhs);
      }();
      benchmark::DoNotOptimize(result);
    }
    if (!valid())
    {
      spatialDotFailed = true;
      state.SkipWithError("Timed spatial dot changed the result");
    }
    state.counters["dimension"] = lhs.size();
  }

  /**
   * @brief Registers both scalar orders without conflating mixed contractions.
   * @tparam LHS Left operand's scalar type.
   * @tparam RHS Right operand's scalar type.
   * @param label Scalar-pair label used in benchmark names.
   */
  template <class LHS, class RHS>
  void registerSpatialDot(const char* label)
  {
    const auto registerCase = [&]<bool Original>(const char* operation) {
      const std::string name = std::string("SpatialDot/") + operation + label;
      benchmark::RegisterBenchmark(
        name.c_str(), [](auto& state) { spatialDot<LHS, RHS, Original>(state); })
        ->Arg(0)
        ->Arg(1)
        ->Arg(2)
        ->Arg(3)
        ->UseRealTime();
    };
    registerCase.template operator()<false>("Public/");
    registerCase.template operator()<true>("Original/");
  }
}

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Tests::Benchmarks;
  registerSpatialDot<Real, Real>("RealReal");
  registerSpatialDot<Real, Complex>("RealComplex");
  registerSpatialDot<Complex, Real>("ComplexReal");
  registerSpatialDot<Complex, Complex>("ComplexComplex");
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  benchmark::AddCustomContext("scope", "Spatial dot only; no geometry or assembly");
  const size_t cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || spatialDotFailed ? 1 : 0;
}
