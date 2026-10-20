/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Conjugating spatial-matrix contractions and weighted pair batches.
 */
#include <benchmark/benchmark.h>
#include <algorithm>
#include <limits>
#include <type_traits>
#include <vector>
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Common.h"

namespace Rodin::Tests::Benchmarks
{
  bool matrixContractionFailed = false;

  /// @brief Retained matrix occupancy and basis-sized pair-count protocols.
  enum class MatrixContractionInput
  {
    Dense,
    ComponentDense,
    ComponentSparse
  };

  /**
   * @brief Measures the public Frobenius product used by named elasticity.
   *
   * ## Mathematical operation
   *
   * For retained matrices with logical shape @f$r\times c@f$, the scalar
   * scope evaluates @f$\sum_{i,j} A_{ij}\overline{B_{ij}}@f$. The batch scope
   * evaluates @f$K_{ba}=w\sum_{i,j} A^a_{ij}\overline{B^b_{ij}}@f$ with
   * @f$w=3/4@f$ and @f$N^2@f$ pairs. Output zeroing and weighted accumulation
   * are timed. Setup, exact qualification and perturbation checks are untimed.
   *
   * ## Measurement boundary
   *
   * Original uses const dynamic Eigen blocks, preserving the existing public
   * reduction expression and reading dimensions from the operands. KnownColumns
   * selects three-column views before the batch pair loop; the row extent and
   * total size remain dynamic to retain the original reduction strategy.
   * Other shapes retain the public spatial-matrix dot call. This experimental
   * caller-known shape protocol models a local assembly consumer, not an
   * optimization of the generic spatial-matrix dot API.
   * All scopes retain runtime dimensions and poison
   * inactive storage. Bounded dyadic inputs admit an independent exact sum,
   * including mixed real/complex fields and empty rectangular shapes.
   * No operand symmetry is assumed. Counts 6, 12 and 30 cover representative
   * vector-element pair-loop sizes, not complete finite element kernels.
   * Geometry, coefficient evaluation, Jacobian mapping and sparse insertion
   * are excluded. The named-form benchmark supplies the assembly consumer.
   *
   * Component protocols retain @f$mN@f$ matrices for @f$m=rc@f$ components
   * and @f$N\in\{6,27,64\}@f$ scalar-basis-count proxies. ComponentSparse
   * stores one active component @f$a\bmod m@f$ per matrix. Its independent
   * oracle uses the component Kronecker delta, rather than reducing the input
   * matrices. ComponentDense retains the same dimensions and pair counts,
   * but fills every active entry. These are arithmetic inputs resembling
   * component bases, not evaluations of a particular element or geometry.
   * Empty shapes remain structural bounds, not physical refinement cases.
   *
   * @tparam Left Left matrix scalar type.
   * @tparam Right Right matrix scalar type.
   * @tparam Original Whether to use the original dynamic-block expression.
   * @tparam Batch Whether to measure weighted pair accumulation.
   * @tparam KnownColumns Whether the caller dispatches three columns once per batch.
   * @tparam Input Matrix occupancy and pair-count protocol.
   */
  template <class Left, class Right, bool Original, bool Batch, bool KnownColumns = false,
    MatrixContractionInput Input = MatrixContractionInput::Dense>
  class SpatialMatrixContraction
  {
    public:
      /**
       * @brief Runs one qualified shape and pair-count protocol.
       * @param state Rows, columns, pair count and timing repetitions.
       */
      static void run(benchmark::State& state)
      {
        using Scalar = typename FormLanguage::Mult<Left, Right>::Type;
        const size_t rows = state.range(0), columns = state.range(1);
        const size_t count = state.range(2);
        const size_t components = rows * columns;
        constexpr Real Weight = Real(3) / 4; ///< Dimensionless accumulation weight.
        std::vector<Math::SpatialMatrix<Left>> left(count);
        std::vector<Math::SpatialMatrix<Right>> right(count);
        const auto fill = [&]<class Entry>(
                            Math::SpatialMatrix<Entry>& matrix, size_t a, size_t side) {
          matrix.resize(rows, columns);
          matrix.getData().setConstant(Entry(std::numeric_limits<Real>::quiet_NaN()));
          for (size_t j = 0; j < columns; ++j)
          {
            for (size_t i = 0; i < rows; ++i)
            {
              if constexpr (Input == MatrixContractionInput::ComponentSparse)
              {
                if (components == 0 || i * columns + j != a % components)
                {
                  matrix(i, j) = Entry(0);
                  continue;
                }
              }
              const Real real = Real(a + 2 * i + 3 * j + side + 1) / 16;
              if constexpr (std::is_same_v<Entry, Complex>)
                matrix(i, j) = Entry(real, Real(i + a + 1) / (side == 0 ? 32 : -32));
              else
                matrix(i, j) = real;
            }
          }
        };
        for (size_t a = 0; a < count; ++a)
        {
          fill(left[a], a, 0);
          fill(right[a], a, 1);
        }
        const auto original = [](const auto& lhs, const auto& rhs) {
          static_assert(
            std::is_const_v<std::remove_reference_t<decltype(lhs.getData())>>);
          static_assert(
            std::is_const_v<std::remove_reference_t<decltype(rhs.getData())>>);
          const auto activeRows = static_cast<Eigen::Index>(lhs.rows());
          const auto activeColumns = static_cast<Eigen::Index>(lhs.cols());
          return Math::dot(lhs.getData().topLeftCorner(activeRows, activeColumns),
            rhs.getData().topLeftCorner(activeRows, activeColumns));
        };
        const auto product = [&]<int Columns = Eigen::Dynamic>(
                               const auto& lhs, const auto& rhs) {
          if constexpr (Original)
            return original(lhs, rhs);
          else if constexpr (KnownColumns && Columns != Eigen::Dynamic)
          {
            return Math::dot(
              lhs.getData().template topLeftCorner<Eigen::Dynamic, Columns>(
                lhs.rows(), lhs.cols()),
              rhs.getData().template topLeftCorner<Eigen::Dynamic, Columns>(
                rhs.rows(), rhs.cols()));
          }
          else
            return Math::dot(lhs, rhs);
        };
        Math::Matrix<Scalar> expected(count, count), output(count, count);
        for (size_t b = 0; b < count; ++b)
        {
          for (size_t a = 0; a < count; ++a)
          {
            Scalar sum = 0;
            if constexpr (Input == MatrixContractionInput::ComponentSparse)
            {
              if (components != 0 && a % components == b % components)
              {
                const size_t component = a % components;
                const size_t i = component / columns, j = component % columns;
                const auto lhs = [&]() {
                  const Real real = Real(a + 2 * i + 3 * j + 1) / 16;
                  if constexpr (std::is_same_v<Left, Complex>)
                    return Left(real, Real(i + a + 1) / 32);
                  else
                    return Left(real);
                }();
                const auto rhs = [&]() {
                  const Real real = Real(b + 2 * i + 3 * j + 2) / 16;
                  if constexpr (std::is_same_v<Right, Complex>)
                    return Right(real, -Real(i + b + 1) / 32);
                  else
                    return Right(real);
                }();
                sum = lhs * Math::conj(rhs);
              }
            }
            else
            {
              for (size_t j = 0; j < columns; ++j)
              {
                for (size_t i = 0; i < rows; ++i)
                  sum += left[a](i, j) * Math::conj(right[b](i, j));
              }
            }
            if (Math::dot(left[a], right[b]) != sum || original(left[a], right[b]) != sum)
            {
              matrixContractionFailed = true;
              state.SkipWithError("Matrix contraction differs from exact dyadic sum");
              return;
            }
            expected(b, a) = Weight * sum;
          }
        }
        const auto evaluateShape = [&]<int Columns>() {
          output.setZero();
          for (size_t b = 0; b < count; ++b)
          {
            for (size_t a = 0; a < count; ++a)
            {
              output(b, a) +=
                Weight * product.template operator()<Columns>(left[a], right[b]);
            }
          }
        };
        const auto evaluate = [&] {
          if constexpr (KnownColumns)
          {
            if (columns == 3)
              evaluateShape.template operator()<3>();
            else
              evaluateShape.template operator()<Eigen::Dynamic>();
          }
          else
            evaluateShape.template operator()<Eigen::Dynamic>();
        };
        const auto valid = [&] { return (output.array() == expected.array()).all(); };
        evaluate();
        if (!valid())
        {
          matrixContractionFailed = true;
          state.SkipWithError("Weighted batch differs from exact entries");
          return;
        }
        output(0, 0) += Scalar(1);
        if (valid())
        {
          matrixContractionFailed = true;
          state.SkipWithError("Batch oracle accepted a perturbed output");
          return;
        }
        for (auto _ : state)
        {
          benchmark::DoNotOptimize(left.data());
          benchmark::DoNotOptimize(right.data());
          if constexpr (Batch)
          {
            evaluate();
            benchmark::DoNotOptimize(output.data());
          }
          else
          {
            auto result = [&] {
              if constexpr (KnownColumns)
              {
                if (columns == 3)
                  return product.template operator()<3>(left[0], right[0]);
              }
              return product(left[0], right[0]);
            }();
            benchmark::DoNotOptimize(result);
          }
          benchmark::ClobberMemory();
        }
        if constexpr (!Batch)
          evaluate();
        if (!valid())
        {
          matrixContractionFailed = true;
          state.SkipWithError("Timed batch or inputs no longer match exact entries");
        }
        state.counters["rows"] = rows;
        state.counters["columns"] = columns;
        state.counters["pairs"] = Batch ? count * count : 1;
        state.counters["entries_per_pair"] = rows * columns;
        if constexpr (Input != MatrixContractionInput::Dense)
        {
          state.counters["field_components"] = components;
          state.counters["basis_count_per_component"] =
            components == 0 ? 0 : count / components;
          state.counters["nonzero_entries_per_input"] =
            Input == MatrixContractionInput::ComponentSparse ? (components == 0 ? 0 : 1)
                                                             : components;
          state.counters["nonzero_pair_outputs"] = components == 0
            ? 0
            : (Input == MatrixContractionInput::ComponentSparse
                  ? count * count / components
                  : count * count);
        }
        state.SetItemsProcessed(state.iterations() * (Batch ? count * count : 1));
      }

      /**
       * @brief Registers every bounded shape and the intended pair counts.
       * @param field Scalar-pair label.
       */
      static void registerCases(const char* field)
      {
        const auto name = std::string("SpatialMatrixContraction/") +
          (Input == MatrixContractionInput::ComponentDense       ? "ComponentDense/"
              : Input == MatrixContractionInput::ComponentSparse ? "ComponentSparse/"
              : Batch                                            ? "PairBatch/"
                                                                 : "Scalar/") +
          (KnownColumns  ? "KnownColumns/"
              : Original ? "Original/"
                         : "Public/") +
          field;
        auto* registration = benchmark::RegisterBenchmark(name.c_str(), &run);
        for (size_t rows = 0; rows <= 3; ++rows)
        {
          for (size_t columns = 0; columns <= 3; ++columns)
          {
            if constexpr (Batch)
            {
              if constexpr (Input != MatrixContractionInput::Dense)
              {
                for (const int64_t basisCount : {6, 27, 64})
                {
                  const int64_t components = static_cast<int64_t>(rows * columns);
                  registration->Args(
                    {static_cast<int64_t>(rows), static_cast<int64_t>(columns),
                      std::max(int64_t(1), components) * basisCount});
                }
              }
              else
              {
                for (const int64_t count : {6, 12, 30})
                {
                  registration->Args(
                    {static_cast<int64_t>(rows), static_cast<int64_t>(columns), count});
                }
              }
            }
            else
            {
              registration->Args(
                {static_cast<int64_t>(rows), static_cast<int64_t>(columns), 1});
            }
          }
        }
        registration->UseRealTime();
      }
  };

  /**
   * @brief Registers public and original protocols for a scalar pair.
   * @tparam Left Left scalar field.
   * @tparam Right Right scalar field.
   * @param field Scalar-pair label.
   */
  template <class Left, class Right>
  void registerMatrixContractions(const char* field)
  {
    SpatialMatrixContraction<Left, Right, false, false>::registerCases(field);
    SpatialMatrixContraction<Left, Right, true, false>::registerCases(field);
    SpatialMatrixContraction<Left, Right, false, true>::registerCases(field);
    SpatialMatrixContraction<Left, Right, true, true>::registerCases(field);
    SpatialMatrixContraction<Left, Right, false, false, true>::registerCases(field);
    SpatialMatrixContraction<Left, Right, false, true, true>::registerCases(field);
    SpatialMatrixContraction<Left, Right, false, true, false,
      MatrixContractionInput::ComponentDense>::registerCases(field);
    SpatialMatrixContraction<Left, Right, true, true, false,
      MatrixContractionInput::ComponentDense>::registerCases(field);
    SpatialMatrixContraction<Left, Right, false, true, false,
      MatrixContractionInput::ComponentSparse>::registerCases(field);
    SpatialMatrixContraction<Left, Right, true, true, false,
      MatrixContractionInput::ComponentSparse>::registerCases(field);
  }
}

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Tests::Benchmarks;
  registerMatrixContractions<Real, Real>("RealReal");
  registerMatrixContractions<Real, Complex>("RealComplex");
  registerMatrixContractions<Complex, Real>("ComplexReal");
  registerMatrixContractions<Complex, Complex>("ComplexComplex");
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  benchmark::AddCustomContext(
    "scope", "Conjugating matrix product; weighted pair batches include zeroing");
  const auto cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || matrixContractionFailed ? 1 : 0;
}
