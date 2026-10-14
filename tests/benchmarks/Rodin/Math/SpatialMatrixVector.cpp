/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Runtime-shape spatial matrix-vector products used by gradient mapping.
 */
#include <benchmark/benchmark.h>
#include <cmath>
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/SpatialMatrix.h"

namespace Rodin::Tests::Benchmarks
{
  bool spatialMatrixVectorFailed = false;

  /**
   * @brief Measures the canonical fixed-capacity product with runtime extents.
   *
   * ## Architecture
   *
   * All rectangular shapes with zero through three rows/columns and all
   * real/complex operand pairs are measured. Input padding is poisoned; the
   * logical contraction must not consume it and inactive result storage must
   * be zero. Bounded dyadic operands give an exact independent reference sum,
   * including genuinely complex values and mixed scalar promotion.
   *
   * Timed calls include result construction and its storage initialization.
   * Escaping operand/result addresses and memory barriers prevent constant
   * folding and discarded padding stores. Guards and operand setup are untimed.
   * The owning Variational physical-gradient benchmark measures the consumer;
   * these primitive timings are not assembly or gradient-sweep timings.
   *
   * @tparam MatrixScalar Matrix entry field.
   * @tparam VectorScalar Vector entry field.
   */
  template <class MatrixScalar, class VectorScalar>
  class SpatialMatrixVector
  {
    public:
    /**
     * @brief Measures one qualified product shape.
     * @param state Runtime rows, contraction dimension and timing protocol.
     */
    static void run(benchmark::State& state)
    {
      using Out = typename FormLanguage::Mult<MatrixScalar, VectorScalar>::Type;
      const size_t rows = state.range(0), columns = state.range(1);
      Math::SpatialMatrix<MatrixScalar> matrix(rows, columns);
      Math::SpatialVector<VectorScalar> vector(columns);
      matrix.getData().setConstant(MatrixScalar(128));
      vector.getData().setConstant(VectorScalar(256));
      const auto coefficient = []<class Scalar>(size_t i, size_t j) {
        const Real real = Real(2 * i + j + 1) / 8;
        if constexpr (std::is_same_v<Scalar, Complex>)
          return Scalar(real, -Real(i + 2 * j + 1) / 16);
        else
          return real;
      };
      for (size_t i = 0; i < rows; ++i)
      {
        for (size_t j = 0; j < columns; ++j)
          matrix(i, j) = coefficient.template operator()<MatrixScalar>(i, j);
      }
      for (size_t j = 0; j < columns; ++j)
        vector(j) = coefficient.template operator()<VectorScalar>(j, 1);
      Math::SpatialVector<Out> expected(rows);
      for (size_t i = 0; i < rows; ++i)
      {
        for (size_t j = 0; j < columns; ++j)
          expected(i) += Out(matrix(i, j)) * Out(vector(j));
      }
      const auto valid = [&](const auto& result) {
        if (result.size() != rows)
          return false;
        for (size_t i = 0; i < 3; ++i)
        {
          const Out value = result.getData()(i);
          if (!std::isfinite(std::abs(value)) ||
            value != (i < rows ? expected(i) : Out(0)))
            return false;
        }
        return true;
      };
      auto result = matrix * vector;
      auto wrong = result;
      if (rows > 0)
        wrong(0) += Out(1);
      else
        wrong.resize(1);
      if (!valid(result) || valid(wrong))
      {
        spatialMatrixVectorFailed = true;
        state.SkipWithError("Product, dimensions, padding or sensitivity guard failed");
        return;
      }
      for (auto _ : state)
      {
        benchmark::DoNotOptimize(&matrix);
        benchmark::DoNotOptimize(&vector);
        auto product = matrix * vector;
        benchmark::DoNotOptimize(&product);
        benchmark::ClobberMemory();
      }
      result = matrix * vector;
      if (!valid(result))
      {
        spatialMatrixVectorFailed = true;
        state.SkipWithError("Repeated product changed a qualified component");
      }
      state.counters["rows"] = rows;
      state.counters["columns"] = columns;
      state.counters["scalar_products"] = rows * columns;
      state.counters["result_bytes"] = sizeof(result);
    }

    /// @brief Registers all runtime shapes for the operand scalar pair.
    static void registerCases()
    {
      const std::string name = std::string("SpatialMatrixVector/") +
        (std::is_same_v<MatrixScalar, Complex> ? "Complex" : "Real") + "By" +
        (std::is_same_v<VectorScalar, Complex> ? "Complex" : "Real");
      benchmark::RegisterBenchmark(name.c_str(), &run)
        ->ArgsProduct({{0, 1, 2, 3}, {0, 1, 2, 3}})
        ->UseRealTime();
    }
  };
}

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Tests::Benchmarks;
  SpatialMatrixVector<Real, Real>::registerCases();
  SpatialMatrixVector<Real, Complex>::registerCases();
  SpatialMatrixVector<Complex, Real>::registerCases();
  SpatialMatrixVector<Complex, Complex>::registerCases();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  benchmark::AddCustomContext(
    "scope", "Canonical product including result initialization");
  const auto cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || spatialMatrixVectorFailed ? 1 : 0;
}
