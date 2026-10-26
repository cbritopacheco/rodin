/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Spatial matrix products underlying vector-basis gradient transforms.
 */
#include <benchmark/benchmark.h>
#include <cmath>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include "Rodin/Math/SpatialMatrix.h"

namespace Rodin::Tests::Benchmarks
{
  bool spatialMatrixProductFailed = false;

  /**
   * @brief Measures the canonical matrix product over retained basis-sized arrays.
   *
   * ## Mathematical operation
   *
   * For @f$A_a\in\mathbb F^{r\times k}@f$ and
   * @f$B\in\mathbb G^{k\times c}@f$, the operation is @f$C_a=A_aB@f$.
   * There is no conjugation in this product. Reference Jacobians have
   * @f$r=m@f$, @f$k=d@f$ and the geometric inverse has @f$c=s@f$, where
   * @f$m@f$ is the field component count and @f$d,s@f$ are reference and
   * ambient dimensions. The real/real and complex/real selections therefore
   * isolate the product used by vector-basis pullbacks in assembly. The other
   * scalar pairs cover the same public arithmetic contract.
   *
   * ## Measurement boundary
   *
   * Each iteration multiplies every retained left matrix by the same right
   * matrix, including result construction and assignment to retained output.
   * Counts 6, 12, 30 and 81 model local vector-basis sizes; they do not identify
   * particular geometries. No finite-element evaluation, geometric inverse,
   * symmetric-gradient construction, coefficient weighting, pair contraction
   * or sparse insertion is included. This is a primitive measurement, not a
   * decomposition into additive fractions of a complete cell-kernel time.
   *
   * All dimensions zero through three and all real/complex scalar pairs are
   * qualified. Bounded dyadic inputs admit an exact independent scalar-loop
   * reference. Inactive input storage is poisoned, while inactive output must
   * remain zero. Wrong values or dimensions must be rejected outside timing.
   *
   * @tparam Left Left matrix scalar type.
   * @tparam Right Right matrix scalar type.
   */
  template <class Left, class Right>
  class SpatialMatrixProduct
  {
    public:
      /**
       * @brief Measures one runtime shape and array length.
       * @param state Rows, contraction dimension, columns and matrix count.
       */
      static void run(benchmark::State& state)
      {
        using Scalar = typename FormLanguage::Mult<Left, Right>::Type;
        using Matrix = Math::SpatialMatrix<Scalar>;
        const size_t rows = state.range(0), inner = state.range(1);
        const size_t columns = state.range(2), count = state.range(3);
        std::vector<Math::SpatialMatrix<Left>> left(count);
        Math::SpatialMatrix<Right> right(inner, columns);
        right.getData().setConstant(Right(std::numeric_limits<Real>::quiet_NaN()));
        const auto coefficient = []<class Entry>(size_t a, size_t i, size_t j) {
          const Real real = Real((a + 2 * i + j) % 7 + 1) / 8;
          if constexpr (std::is_same_v<Entry, Complex>)
            return Entry(real, -Real((a + i + 3 * j) % 5 + 1) / 16);
          else
            return real;
        };
        for (size_t a = 0; a < count; ++a)
        {
          left[a].resize(rows, inner);
          left[a].getData().setConstant(Left(std::numeric_limits<Real>::quiet_NaN()));
          for (size_t i = 0; i < rows; ++i)
          {
            for (size_t k = 0; k < inner; ++k)
              left[a](i, k) = coefficient.template operator()<Left>(a, i, k);
          }
        }
        for (size_t k = 0; k < inner; ++k)
        {
          for (size_t j = 0; j < columns; ++j)
            right(k, j) = coefficient.template operator()<Right>(1, k, j);
        }
        std::vector<Matrix> expected(count), output(count);
        for (size_t a = 0; a < count; ++a)
        {
          expected[a].resize(rows, columns);
          for (size_t i = 0; i < rows; ++i)
          {
            for (size_t j = 0; j < columns; ++j)
            {
              for (size_t k = 0; k < inner; ++k)
                expected[a](i, j) += Scalar(left[a](i, k)) * Scalar(right(k, j));
            }
          }
        }
        const auto valid = [&](const auto& matrices) {
          if (matrices.size() != count)
            return false;
          for (size_t a = 0; a < count; ++a)
          {
            const auto& result = matrices[a];
            if (result.rows() != rows || result.cols() != columns)
              return false;
            for (size_t i = 0; i < 3; ++i)
            {
              for (size_t j = 0; j < 3; ++j)
              {
                const Scalar value = result.getData()(i, j);
                const Scalar reference =
                  i < rows && j < columns ? expected[a](i, j) : Scalar(0);
                if (!std::isfinite(std::abs(value)) || value != reference)
                  return false;
              }
            }
          }
          return true;
        };
        const auto evaluate = [&] {
          for (size_t a = 0; a < count; ++a)
            output[a] = left[a] * right;
        };
        evaluate();
        auto wrong = output;
        if (rows > 0 && columns > 0)
          wrong.front()(0, 0) += Scalar(1);
        else
          wrong.front().resize(rows == 0 ? 1 : rows, columns == 0 ? 1 : columns);
        if (!valid(output) || valid(wrong))
        {
          spatialMatrixProductFailed = true;
          state.SkipWithError(
            "Matrix product, shape, padding or sensitivity guard failed");
          return;
        }
        for (auto _ : state)
        {
          benchmark::DoNotOptimize(left.data());
          benchmark::DoNotOptimize(&right);
          evaluate();
          benchmark::DoNotOptimize(output.data());
          benchmark::ClobberMemory();
        }
        if (!valid(output))
        {
          spatialMatrixProductFailed = true;
          state.SkipWithError("Repeated matrix products changed qualified values");
        }
        state.counters["rows"] = rows;
        state.counters["inner_dimension"] = inner;
        state.counters["columns"] = columns;
        state.counters["matrix_count"] = count;
        state.counters["scalar_products"] = rows * inner * columns * count;
        state.counters["output_bytes"] = count * sizeof(Matrix);
      }

      /// @brief Registers all runtime shapes and five array lengths.
      static void registerCases()
      {
        const std::string name = std::string("SpatialMatrixProduct/") +
          (std::is_same_v<Left, Complex> ? "Complex" : "Real") + "By" +
          (std::is_same_v<Right, Complex> ? "Complex" : "Real");
        benchmark::RegisterBenchmark(name.c_str(), &run)
          ->ArgsProduct({{0, 1, 2, 3}, {0, 1, 2, 3}, {0, 1, 2, 3}, {1, 6, 12, 30, 81}})
          ->ArgNames({"rows", "inner", "columns", "count"});
      }
  };
}

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Tests::Benchmarks;
  SpatialMatrixProduct<Real, Real>::registerCases();
  SpatialMatrixProduct<Real, Complex>::registerCases();
  SpatialMatrixProduct<Complex, Real>::registerCases();
  SpatialMatrixProduct<Complex, Complex>::registerCases();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  benchmark::AddCustomContext("scope", "Public matrix product into retained arrays");
  const auto cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || spatialMatrixProductFailed ? 1 : 0;
}
