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
#include <limits>
#include <type_traits>
#include <vector>
#include "Rodin/Math/SpatialVector.h"

namespace Rodin::Tests::Benchmarks
{
  bool spatialDotFailed = false;

  /**
   * @brief Measures the scalar contraction used inside diffusion kernels.
   *
   * The public operation computes @f$\sum_i u_i\overline{v_i}@f$ on the
   * active entries. The reference explicitly calls the original dynamic Eigen
   * reduction through the original const Eigen views. Compile-time guards
   * retain that expression-type contract. Both operate on the same vectors;
   * construction and exact result guards are outside timing. This does not
   * measure the differently defined
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
    const auto& left = lhs;
    const auto& right = rhs;
    const auto reference = [&]() {
      static_assert(std::is_const_v<std::remove_reference_t<decltype(left.getData())>>);
      static_assert(std::is_const_v<std::remove_reference_t<decltype(right.getData())>>);
      const auto size = static_cast<Eigen::Index>(left.size());
      return Math::dot(left.getData().head(size), right.getData().head(size));
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

  /**
   * @brief Measures public spatial contractions over a retained pair batch.
   *
   * ## Architecture
   *
   * For @f$N@f$ left and right vectors of logical dimension @f$d@f$, stores
   * @f$D_{ba}=\sum_{j=0}^{d-1}u_{bj}\overline{v_{aj}}@f$. Full batches
   * evaluate all @f$N^2@f$ pairs; lower batches evaluate @f$N(N+1)/2@f$ pairs
   * without assuming that the arbitrary input vectors define a symmetric
   * operator. Barriers surround a complete batch, not each scalar contraction.
   * Construction and exact entry qualification remain outside timing.
   *
   * The original dynamic Eigen reduction through const views supplies every
   * reference entry. Compile-time guards reject mutable comparator views.
   * Inactive vector storage is poisoned, including all storage at dimension
   * zero. A changed output must fail qualification. This measures Math::dot
   * with pair-loop and output-store overhead, not physical gradient mapping,
   * weighted local-matrix accumulation or complete assembly.
   * RowSnapshot scopes copy the current left vector into a local value before
   * traversing its right-hand pairs. Ordinary Public and Original scopes keep
   * a borrowed row. This isolates value-ownership effects without changing the
   * reduction expression or adding a retained cache. The strategy is not an
   * assembly timing decomposition.
   *
   * @tparam LHS Left operand scalar type.
   * @tparam RHS Right operand scalar type.
   * @tparam Original Whether to time the original dynamic reduction.
   * @tparam Lower Whether only lower-index pairs are evaluated.
   * @tparam Snapshot Whether each left row is copied into a local value.
   */
  template <class LHS, class RHS, bool Original, bool Lower, bool Snapshot = false>
  class SpatialDotBatch
  {
    public:
    /**
     * @brief Measures one complete pair batch per iteration.
     * @param state Logical dimension, vector count and repetition protocol.
     */
    static void run(benchmark::State& state)
    {
      const size_t dimension = state.range(0);
      const size_t count = state.range(1);
      std::vector<Math::SpatialVector<LHS>> lhs(count);
      std::vector<Math::SpatialVector<RHS>> rhs(count);
      for (size_t b = 0; b < count; ++b)
      {
        lhs[b].resize(dimension);
        rhs[b].resize(dimension);
        lhs[b].getData().setConstant(LHS(std::numeric_limits<Real>::quiet_NaN()));
        rhs[b].getData().setConstant(RHS(std::numeric_limits<Real>::quiet_NaN()));
        for (size_t j = 0; j < dimension; ++j)
        {
          if constexpr (std::is_same_v<LHS, Complex>)
            lhs[b](j) = LHS(Real(b + j + 1) / 7, Real(b + 2 * j + 2) / 11);
          else
            lhs[b](j) = Real(b + j + 1) / 7;
          if constexpr (std::is_same_v<RHS, Complex>)
            rhs[b](j) = RHS(Real(2 * b + j + 2) / 13, -Real(b + j + 1) / 17);
          else
            rhs[b](j) = Real(2 * b + j + 2) / 13;
        }
      }
      const auto& left = lhs;
      const auto& right = rhs;
      static_assert(
        std::is_const_v<std::remove_reference_t<decltype(left[0].getData())>>);
      static_assert(
        std::is_const_v<std::remove_reference_t<decltype(right[0].getData())>>);
      using Value = decltype(Math::dot(left[0], right[0]));
      std::vector<Value> expected(count * count, Value(0));
      std::vector<Value> actual(count * count, Value(0));
      for (size_t b = 0; b < count; ++b)
      {
        const auto n = static_cast<Eigen::Index>(left[b].size());
        for (size_t a = 0; a < (Lower ? b + 1 : count); ++a)
        {
          expected[b * count + a] =
            Math::dot(left[b].getData().head(n), right[a].getData().head(n));
        }
      }
      const auto compute = [&]() {
        for (size_t b = 0; b < count; ++b)
        {
          const auto& row = [&]() -> decltype(auto) {
            if constexpr (Snapshot)
              return Math::SpatialVector<LHS>(left[b]);
            else
              return (left[b]);
          }();
          for (size_t a = 0; a < (Lower ? b + 1 : count); ++a)
          {
            if constexpr (Original)
            {
              const auto n = static_cast<Eigen::Index>(row.size());
              actual[b * count + a] =
                Math::dot(row.getData().head(n), right[a].getData().head(n));
            }
            else
              actual[b * count + a] = Math::dot(row, right[a]);
          }
        }
      };
      const auto valid = [&]() { return actual == expected; };
      compute();
      if (!valid())
      {
        spatialDotFailed = true;
        state.SkipWithError("Batched spatial dot differs from original entry reduction");
        return;
      }
      actual[0] += Value(1);
      if (valid())
      {
        spatialDotFailed = true;
        state.SkipWithError("Batched spatial dot failed to reject a changed entry");
        return;
      }
      actual[0] = expected[0];
      for (auto _ : state)
      {
        benchmark::DoNotOptimize(lhs.data());
        benchmark::DoNotOptimize(rhs.data());
        compute();
        benchmark::DoNotOptimize(actual.data());
        benchmark::ClobberMemory();
      }
      if (!valid())
      {
        spatialDotFailed = true;
        state.SkipWithError("Timed spatial dot batch changed a reference entry");
      }
      const size_t pairs = Lower ? count * (count + 1) / 2 : count * count;
      state.counters["dimension"] = dimension;
      state.counters["row_snapshot"] = Snapshot;
      state.counters["vectors_per_operand"] = count;
      state.counters["pairs"] = pairs;
      state.counters["scalar_products"] = pairs * dimension;
      state.counters["input_bytes"] = count * (sizeof(lhs[0]) + sizeof(rhs[0]));
      state.counters["output_bytes"] = actual.size() * sizeof(Value);
      state.counters["written_bytes"] = pairs * sizeof(Value);
    }

    /**
     * @brief Registers empty and nonempty dimensions at three batch sizes.
     * @param label Scalar-pair label used in benchmark names.
     */
    static void registerCases(const char* label)
    {
      const std::string name = std::string("SpatialDot/Batch/") +
        (Snapshot      ? "RowSnapshot/"
            : Original ? "Original/"
                       : "Public/") +
        (Lower ? "Lower/" : "Full/") + label;
      for (const size_t dimension : {size_t(0), size_t(1), size_t(2), size_t(3)})
      {
        for (const size_t count : {size_t(4), size_t(16), size_t(64)})
        {
          benchmark::RegisterBenchmark(name.c_str(), [](auto& state) { run(state); })
            ->Args({static_cast<int64_t>(dimension), static_cast<int64_t>(count)})
            ->UseRealTime();
        }
      }
    }
  };

  /**
   * @brief Registers paired implementations and both pair-loop shapes.
   * @tparam LHS Left operand scalar type.
   * @tparam RHS Right operand scalar type.
   * @param label Scalar-pair label used in benchmark names.
   */
  template <class LHS, class RHS>
  void registerSpatialDotBatches(const char* label)
  {
    SpatialDotBatch<LHS, RHS, false, false>::registerCases(label);
    SpatialDotBatch<LHS, RHS, true, false>::registerCases(label);
    SpatialDotBatch<LHS, RHS, false, true>::registerCases(label);
    SpatialDotBatch<LHS, RHS, true, true>::registerCases(label);
    SpatialDotBatch<LHS, RHS, false, false, true>::registerCases(label);
    SpatialDotBatch<LHS, RHS, false, true, true>::registerCases(label);
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
  registerSpatialDotBatches<Real, Real>("RealReal");
  registerSpatialDotBatches<Real, Complex>("RealComplex");
  registerSpatialDotBatches<Complex, Real>("ComplexReal");
  registerSpatialDotBatches<Complex, Complex>("ComplexComplex");
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  benchmark::AddCustomContext("scope", "Spatial dot only; no geometry or assembly");
  const size_t cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || spatialDotFailed ? 1 : 0;
}
