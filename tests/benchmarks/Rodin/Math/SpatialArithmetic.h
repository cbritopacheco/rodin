/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Matched spatial arithmetic workloads and independent numerical guards.
 */
#ifndef RODIN_RODIN_MATH_SPATIALARITHMETIC_H
#define RODIN_RODIN_MATH_SPATIALARITHMETIC_H

#include <benchmark/benchmark.h>
#include <limits>
#include <string>
#include <vector>
#include <utility>
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialVector.h"

namespace Rodin::Tests::Benchmarks
{
  /// @brief Storage and extent knowledge available to the timed operation.
  enum class SpatialStorage
  {
    Spatial,
    Bounded,
    Dynamic,
    Fixed
  };

  /// @brief Public value operations with distinct algebraic conventions.
  enum class SpatialOperation
  {
    Zero,
    Copy,
    Add,
    Subtract,
    Negate,
    ScaleLeft,
    ScaleRight,
    Divide,
    AddAssign,
    SubtractAssign,
    ScaleAssign,
    DivideAssign,
    MathSum,
    MathSubtract,
    MathNegate,
    MathScale,
    MathDivide,
    MathProduct,
    ScalarDot,
    MathDot,
    RightDot,
    BilinearProduct,
    SquaredNorm,
    Norm,
    Conjugate,
    MathConjugate,
    Transpose,
    StableNorm,
    BlueNorm,
    LpNorm1,
    LpNorm2,
    LpNorm3,
    Normalized,
    Normalize,
    Cross,
    Adjoint,
    Trace,
    Determinant,
    Inverse,
    Solve,
    PseudoInverse,
    Row,
    Column,
    MatrixVector,
    VectorMatrix,
    MatrixProduct
  };

  /**
   * @brief Compares owned results over identical retained input arrays.
   *
   * ## Architecture
   *
   * Every case prepares nonsymmetric dyadic inputs and an independent reference
   * outside timing, checks a perturbed reference is rejected, measures a complete
   * sweep, and checks all outputs again. Shapes remain runtime values for Spatial,
   * Bounded and Dynamic; Fixed grants the compiler additional extent knowledge.
   * Spatial padding is NaN-poisoned. Reference contractions, cross products,
   * transposes and matrix products use scalar loops with explicit conjugation.
   * Solves/inverses use Eigen full-pivot LU as an untimed independent oracle.
   *
   * Each timed operation constructs an owned result and assigns it to retained
   * output. Dynamic result allocations are included; input/output-array allocation
   * and oracle evaluation are excluded. Assignment scopes include copying the
   * left value before modifying it. Zero means zero construction, not an
   * uninitialized Eigen constructor. Memory barriers surround a complete sweep.
   * Counts one and 32 distinguish isolated calls from amortized loop/barrier costs.
   * Destination variants instead evaluate safe Eigen expressions directly into
   * prepared disjoint output through noalias(). They isolate temporary ownership
   * from arithmetic and are labeled separately from owned-result measurements.
   *
   * MathDot and RightDot compute right-conjugating inner products. BilinearProduct
   * computes an explicitly unconjugated entrywise product. Eigen RightDot uses
   * rhs.dot(lhs), since Eigen's dot conjugates its left operand. Matrix products
   * and Cross are bilinear. Complex Eigen cross is conjugated back to match the
   * spatial algebraic cross. Pseudoinverse comparators use the same full-rank
   * normal-equation algorithm; this is not a rank-revealing solver comparison.
   *
   * @tparam Scalar Real or complex scalar field.
   * @tparam Matrix Whether the measured value is a matrix rather than a vector.
   * @tparam Storage Owning representation used inside timing.
   * @tparam Rows Fixed row count, or Dynamic for runtime extents.
   * @tparam Cols Fixed column count, or Dynamic for runtime matrix extents.
   */
  template <class Scalar, bool Matrix, SpatialStorage Storage, int Rows = Eigen::Dynamic,
    int Cols = Eigen::Dynamic>
  class SpatialArithmetic
  {
    using Spatial = std::conditional_t<Matrix, Math::SpatialMatrix<Scalar>,
      Math::SpatialVector<Scalar>>;
    using EigenValue =
      Eigen::Matrix<Scalar, Storage == SpatialStorage::Fixed ? Rows : Eigen::Dynamic,
        Matrix ? (Storage == SpatialStorage::Fixed ? Cols : Eigen::Dynamic) : 1,
        Matrix && Storage == SpatialStorage::Fixed && Rows == 1 && Cols != 1
          ? Eigen::RowMajor
          : Eigen::ColMajor,
        Storage == SpatialStorage::Dynamic
          ? Eigen::Dynamic
          : (Storage == SpatialStorage::Fixed ? Rows : 3),
        Matrix ? (Storage == SpatialStorage::Dynamic
                     ? Eigen::Dynamic
                     : (Storage == SpatialStorage::Fixed ? Cols : 3))
               : 1>;
    using Value =
      std::conditional_t<Storage == SpatialStorage::Spatial, Spatial, EigenValue>;
    using Reference = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using ReferenceVector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Op = SpatialOperation;

    using Vector = std::conditional_t<Storage == SpatialStorage::Spatial,
      Math::SpatialVector<Scalar>,
      Eigen::Matrix<Scalar,
        Storage == SpatialStorage::Fixed ? (Matrix ? Cols : Rows) : Eigen::Dynamic, 1,
        Eigen::ColMajor,
        Storage == SpatialStorage::Dynamic
          ? Eigen::Dynamic
          : (Storage == SpatialStorage::Fixed ? (Matrix ? Cols : Rows) : 3),
        1>>;
    using RowVector = std::conditional_t<Storage == SpatialStorage::Spatial,
      Math::SpatialVector<Scalar>,
      Eigen::Matrix<Scalar, Storage == SpatialStorage::Fixed ? Rows : Eigen::Dynamic, 1,
        Eigen::ColMajor,
        Storage == SpatialStorage::Dynamic
          ? Eigen::Dynamic
          : (Storage == SpatialStorage::Fixed ? Rows : 3),
        1>>;
    using Product = std::conditional_t<Storage == SpatialStorage::Spatial,
      Math::SpatialMatrix<Scalar>,
      Eigen::Matrix<Scalar, Storage == SpatialStorage::Fixed ? Cols : Eigen::Dynamic,
        Storage == SpatialStorage::Fixed ? Rows : Eigen::Dynamic,
        Storage == SpatialStorage::Fixed && Cols == 1 && Rows != 1 ? Eigen::RowMajor
                                                                   : Eigen::ColMajor,
        Storage == SpatialStorage::Dynamic
          ? Eigen::Dynamic
          : (Storage == SpatialStorage::Fixed ? Cols : 3),
        Storage == SpatialStorage::Dynamic
          ? Eigen::Dynamic
          : (Storage == SpatialStorage::Fixed ? Rows : 3)>>;

    struct Input
    {
      Value left;
      Value right;
      Scalar factor;
      Vector vector;
      RowVector rowVector;
      Product product;
    };

    public:
      /**
       * @brief Measures one shape, input count and operation.
       * @tparam Operation Mathematical operation selected before timing.
       * @tparam Destination Whether Eigen expressions target prepared disjoint output.
       * @param state Runtime shape, count and timing iterations.
       */
      template <Op Operation, bool Destination = false>
      static void run(benchmark::State& state)
      {
        const auto rows = static_cast<std::uint8_t>(state.range(0));
        const auto cols = static_cast<std::uint8_t>(Matrix ? state.range(1) : 1);
        const size_t count = state.range(Matrix ? 2 : 1);
        std::vector<Input> inputs(count);
        for (size_t n = 0; n < count; ++n)
        {
          auto& input = inputs[n];
          input.factor = coefficient(n, 3);
          input.left = make(rows, cols);
          input.right = make(rows, cols);
          input.vector.resize(Matrix ? cols : rows);
          input.rowVector.resize(rows);
          input.product.resize(cols, rows);
          if constexpr (Storage == SpatialStorage::Spatial)
          {
            input.left.getData().setConstant(
              Scalar(std::numeric_limits<Real>::quiet_NaN()));
            input.right.getData().setConstant(
              Scalar(std::numeric_limits<Real>::quiet_NaN()));
          }
          for (size_t i = 0; i < rows; ++i)
          {
            input.rowVector(i) = coefficient(n + i, 1);
            for (size_t j = 0; j < cols; ++j)
            {
              entry(input.left, i, j) =
                coefficient(n + 2 * i + 3 * j, 0) + (i == j ? Scalar(4) : Scalar(0));
              entry(input.right, i, j) = coefficient(n + 3 * i + j, 1);
              input.product(j, i) = entry(input.right, i, j);
            }
          }
          for (size_t i = 0; i < static_cast<size_t>(input.vector.size()); ++i)
            input.vector(i) = coefficient(n + i, 2);
        }
        using Result = decltype(own(evaluate<Operation, Destination>(inputs[0])));
        using Expected = decltype(reference<Operation>(inputs[0]));
        std::vector<Result> output(count);
        std::vector<Expected> expected;
        expected.reserve(count);
        for (const auto& input : inputs)
          expected.push_back(reference<Operation>(input));
        const auto compute = [&] {
          for (size_t i = 0; i < count; ++i)
          {
            if constexpr (Destination && Storage != SpatialStorage::Spatial)
              output[i].noalias() = evaluate<Operation, Destination>(inputs[i]);
            else
              output[i] = evaluate<Operation, Destination>(inputs[i]);
          }
        };
        const auto valid = [&] {
          for (size_t i = 0; i < count; ++i)
          {
            if (!close(output[i], expected[i]))
              return false;
          }
          return true;
        };
        compute();
        if (!valid())
        {
          failed = true;
          state.SkipWithError("Spatial arithmetic disagrees with independent oracle");
          return;
        }
        auto changed = expected[0];
        if constexpr (std::is_arithmetic_v<Expected> || std::is_same_v<Expected, Complex>)
          changed += Expected(1);
        else if (changed.size())
          changed(0, 0) += Scalar(1);
        if constexpr (std::is_arithmetic_v<Expected> || std::is_same_v<Expected, Complex>)
        {
          if (close(output[0], changed))
            failed = true;
        }
        else if (changed.size() && close(output[0], changed))
          failed = true;
        if (failed)
        {
          state.SkipWithError("Spatial arithmetic guard accepted a changed result");
          return;
        }
        for (auto iteration : state)
        {
          benchmark::DoNotOptimize(inputs.data());
          benchmark::ClobberMemory();
          compute();
          benchmark::DoNotOptimize(output.data());
          benchmark::ClobberMemory();
        }
        if (!valid())
        {
          failed = true;
          state.SkipWithError("Timed spatial arithmetic changed a result");
        }
        state.counters["operations"] = count;
        state.counters["input_value_bytes"] = sizeof(Value);
        state.counters["entries"] = rows * cols;
        state.SetItemsProcessed(state.iterations() * count);
      }

      /// @brief Whether any numerical guard failed in this executable.
      inline static bool failed = false;

    private:
      /**
       * @brief Constructs the representation's zero value.
       * @param rows Logical row count (vector size).
       * @param cols Logical column count (one for vectors).
       * @returns An owned zero value.
       */
      static Value make(std::uint8_t rows, std::uint8_t cols)
      {
        if constexpr (Storage == SpatialStorage::Spatial && !Matrix)
          return Value(rows);
        else
        {
          Value value(rows, cols);
          value.setZero();
          return value;
        }
      }

      /**
       * @brief Generates bounded dyadic real or complex entries.
       * @param index Input/component index.
       * @param offset Distinguishes independent operands.
       * @returns Nonzero deterministic coefficient.
       */
      static Scalar coefficient(size_t index, size_t offset)
      {
        const Real real = Real(1 + (index + offset) % 7) / 8;
        if constexpr (std::is_same_v<Scalar, Complex>)
          return Scalar(real, Real(1 + (2 * index + offset) % 5) / 16);
        else
          return real;
      }

      /**
       * @brief Accesses one active entry during untimed setup.
       * @param value Owned input.
       * @param i Row index.
       * @param j Column index (zero for vectors).
       * @returns Active entry reference.
       */
      static Scalar& entry(Value& value, size_t i, size_t j)
      {
        if constexpr (Storage == SpatialStorage::Spatial && !Matrix)
          return value(i);
        else
          return value(i, j);
      }

      /**
       * @brief Borrows the active Eigen view without copying its owner.
       * @param value Immutable input.
       * @returns Active view or Eigen owner reference.
       */
      static decltype(auto) view(const Value& value)
      {
        if constexpr (Storage == SpatialStorage::Spatial && Matrix)
          return value.getData().topLeftCorner(value.rows(), value.cols());
        else if constexpr (Storage == SpatialStorage::Spatial)
          return value.getData().head(value.size());
        else
          return (value);
      }

      /**
       * @brief Materializes Eigen expressions while preserving spatial values.
       * @tparam Destination Whether to retain an expression for direct assignment.
       * @tparam T Expression or owned spatial type.
       * @param value Operation result to materialize.
       * @returns An owned result.
       */
      template <bool Destination = false, class T>
      static decltype(auto) own(T&& value)
      {
        using Argument = std::remove_cvref_t<T>;
        if constexpr (Destination && std::is_same_v<Argument, Value> &&
          std::is_lvalue_reference_v<T>)
          return (value);
        else if constexpr (!Destination && requires { value.eval(); })
        {
          using Owned = std::decay_t<decltype(value.eval())>;
          return Owned(std::forward<T>(value));
        }
        else
          return Argument(std::forward<T>(value));
      }

      /**
       * @brief Executes the selected public operation or scalar-loop comparator.
       * @tparam Operation Timed operation.
       * @tparam Destination Whether to defer safe expressions into disjoint output.
       * @param input Prepared immutable operands.
       * @returns An owned scalar, vector or matrix result.
       */
      template <Op Operation, bool Destination = false>
      static decltype(auto) evaluate(const Input& input)
      {
        const auto& a = input.left;
        const auto& b = input.right;
        const Scalar factor = input.factor;
        if constexpr (Operation == Op::Zero)
        {
          if constexpr (Matrix)
            return make(a.rows(), a.cols());
          else
            return make(a.size(), 1);
        }
        else if constexpr (Operation == Op::Copy)
          return Value(a);
        else if constexpr (Operation == Op::Add)
          return own<Destination>(a + b);
        else if constexpr (Operation == Op::Subtract)
          return own<Destination>(a - b);
        else if constexpr (Operation == Op::Negate)
          return own<Destination>(-a);
        else if constexpr (Operation == Op::ScaleLeft)
          return own<Destination>(factor * a);
        else if constexpr (Operation == Op::ScaleRight)
          return own<Destination>(a * factor);
        else if constexpr (Operation == Op::Divide)
          return own<Destination>(a / factor);
        else if constexpr (Operation == Op::AddAssign ||
          Operation == Op::SubtractAssign || Operation == Op::ScaleAssign ||
          Operation == Op::DivideAssign || Operation == Op::Normalize)
        {
          Value result(a);
          if constexpr (Operation == Op::AddAssign)
            result += b;
          else if constexpr (Operation == Op::SubtractAssign)
            result -= b;
          else if constexpr (Operation == Op::ScaleAssign)
            result *= factor;
          else if constexpr (Operation == Op::DivideAssign)
            result /= factor;
          else
            result.normalize();
          return result;
        }
        else if constexpr (Operation == Op::MathSum)
          return own<Destination>(Math::sum(a, b));
        else if constexpr (Operation == Op::MathSubtract)
          return own<Destination>(Math::minus(a, b));
        else if constexpr (Operation == Op::MathNegate)
          return own<Destination>(Math::minus(a));
        else if constexpr (Operation == Op::MathScale)
          return own<Destination>(Math::mult(a, factor));
        else if constexpr (Operation == Op::MathDivide)
          return own<Destination>(Math::division(a, factor));
        else if constexpr (Operation == Op::ScalarDot)
        {
          const auto& left = view(a);
          const auto& right = view(b);
          Scalar result(0);
          for (Eigen::Index j = 0; j < left.cols(); ++j)
          {
            for (Eigen::Index i = 0; i < left.rows(); ++i)
              result += left(i, j) * Math::conj(right(i, j));
          }
          return result;
        }
        else if constexpr (Operation == Op::MathDot)
          return Math::dot(a, b);
        else if constexpr (Operation == Op::RightDot)
        {
          if constexpr (Storage == SpatialStorage::Spatial)
            return a.dot(b);
          else if constexpr (!Matrix)
            return b.dot(a);
          else
            return (a.array() * b.conjugate().array()).sum();
        }
        else if constexpr (Operation == Op::BilinearProduct)
          return (view(a).array() * view(b).array()).sum();
        else if constexpr (Operation == Op::SquaredNorm)
          return a.squaredNorm();
        else if constexpr (Operation == Op::Norm)
          return a.norm();
        else if constexpr (Operation == Op::Conjugate)
          return own<Destination>(a.conjugate());
        else if constexpr (Operation == Op::MathConjugate)
          return Math::conj(a);
        else if constexpr (Operation == Op::Transpose)
          return own<Destination>(a.transpose());
        else if constexpr (Operation == Op::Adjoint)
          return own<Destination>(a.adjoint());
        else if constexpr (Operation == Op::StableNorm)
          return a.stableNorm();
        else if constexpr (Operation == Op::BlueNorm)
          return a.blueNorm();
        else if constexpr (Operation == Op::LpNorm1)
          return a.template lpNorm<1>();
        else if constexpr (Operation == Op::LpNorm2)
          return a.template lpNorm<2>();
        else if constexpr (Operation == Op::LpNorm3)
          return a.template lpNorm<3>();
        else if constexpr (Operation == Op::Normalized)
          return a.normalized();
        else if constexpr (Operation == Op::Cross)
        {
          if constexpr (Storage == SpatialStorage::Spatial)
            return a.cross(b);
          else
          {
            return own<Destination>(
              a.template head<3>().cross(b.template head<3>()).conjugate());
          }
        }
        else if constexpr (Operation == Op::Trace)
          return a.trace();
        else if constexpr (Operation == Op::Determinant)
          return a.determinant();
        else if constexpr (Operation == Op::Inverse)
          return own<Destination>(a.inverse());
        else if constexpr (Operation == Op::Solve)
        {
          if constexpr (Storage == SpatialStorage::Spatial)
            return a.solve(input.vector);
          else
            return own<Destination>(a.partialPivLu().solve(input.vector));
        }
        else if constexpr (Operation == Op::PseudoInverse)
        {
          if constexpr (Storage == SpatialStorage::Spatial)
            return a.pseudoInverse();
          else
          {
            Value result(a.cols(), a.rows());
            if (a.rows() >= a.cols())
              result = (a.adjoint() * a).eval().inverse() * a.adjoint();
            else
              result = a.adjoint() * (a * a.adjoint()).eval().inverse();
            return result;
          }
        }
        else if constexpr (Operation == Op::Row)
        {
          if constexpr (Storage == SpatialStorage::Spatial)
            return a.row(0);
          else
            return own<Destination>(a.row(0).transpose());
        }
        else if constexpr (Operation == Op::Column)
          return own<Destination>(a.col(0));
        else if constexpr (Operation == Op::MatrixVector)
        {
          if constexpr (Storage == SpatialStorage::Spatial)
            return a * input.vector;
          else
            return own<Destination>(a * input.vector);
        }
        else if constexpr (Operation == Op::VectorMatrix)
        {
          if constexpr (Storage == SpatialStorage::Spatial)
            return input.rowVector * a;
          else
            return own<Destination>(input.rowVector.transpose() * a);
        }
        else
        {
          if constexpr (Storage == SpatialStorage::Spatial)
          {
            if constexpr (Operation == Op::MathProduct)
              return Math::mult(a, input.product);
            else
              return a * input.product;
          }
          else
          {
            if constexpr (Operation == Op::MathProduct)
              return own<Destination>(Math::mult(a, input.product));
            else
              return own<Destination>(a * input.product);
          }
        }
      }

      /**
       * @brief Copies active vector entries for the untimed oracle.
       * @tparam T Vector representation.
       * @param value Input vector.
       * @returns Owned reference vector.
       */
      template <class T>
      static ReferenceVector referenceVector(const T& value)
      {
        if constexpr (requires { value.getData(); })
          return value.getData().head(value.size());
        else
          return value;
      }

      /**
       * @brief Copies active matrix entries for the untimed oracle.
       * @tparam T Matrix representation.
       * @param value Input matrix.
       * @returns Owned reference matrix.
       */
      template <class T>
      static Reference referenceMatrix(const T& value)
      {
        if constexpr (requires { value.getData(); })
          return value.getData().topLeftCorner(value.rows(), value.cols());
        else
          return value;
      }

      /**
       * @brief Evaluates the independent untimed mathematical oracle.
       * @tparam Operation Operation whose contract is checked.
       * @param input Prepared immutable operands.
       * @returns Scalar or dense reference result.
       */
      template <Op Operation>
      static auto reference(const Input& input)
      {
        const Reference a = view(input.left), b = view(input.right);
        const Scalar factor = input.factor;
        if constexpr (Operation == Op::ScalarDot || Operation == Op::MathDot ||
          Operation == Op::RightDot || Operation == Op::BilinearProduct)
        {
          Scalar result(0);
          for (Eigen::Index i = 0; i < a.rows(); ++i)
          {
            for (Eigen::Index j = 0; j < a.cols(); ++j)
            {
              result += a(i, j) *
                (Operation == Op::BilinearProduct ? b(i, j) : Math::conj(b(i, j)));
            }
          }
          return result;
        }
        else if constexpr (Operation == Op::Norm || Operation == Op::SquaredNorm ||
          Operation == Op::StableNorm || Operation == Op::BlueNorm ||
          Operation == Op::LpNorm1 || Operation == Op::LpNorm2 ||
          Operation == Op::LpNorm3)
        {
          Real result = 0;
          constexpr int power = Operation == Op::LpNorm1 ? 1
            : Operation == Op::LpNorm3                   ? 3
                                                         : 2;
          for (Eigen::Index i = 0; i < a.rows(); ++i)
          {
            for (Eigen::Index j = 0; j < a.cols(); ++j)
              result += std::pow(std::abs(a(i, j)), power);
          }
          if constexpr (Operation == Op::SquaredNorm)
            return result;
          else
            return std::pow(result, Real(1) / power);
        }
        else if constexpr (Operation == Op::Trace)
        {
          Scalar result(0);
          for (Eigen::Index i = 0; i < a.rows(); ++i)
            result += a(i, i);
          return result;
        }
        else if constexpr (Operation == Op::Determinant)
        {
          if (!a.rows())
            return Scalar(1);
          return Scalar(a.fullPivLu().determinant());
        }
        else
        {
          Reference result(a.rows(), a.cols());
          if constexpr (Operation == Op::Inverse)
            result = a.fullPivLu().inverse();
          else if constexpr (Operation == Op::Solve)
            result = a.fullPivLu().solve(referenceVector(input.vector));
          else if constexpr (Operation == Op::PseudoInverse)
            result = a.completeOrthogonalDecomposition().pseudoInverse();
          else if constexpr (Operation == Op::MathProduct ||
            Operation == Op::MatrixProduct || Operation == Op::MatrixVector ||
            Operation == Op::VectorMatrix)
          {
            Reference left = a, right = referenceMatrix(input.product);
            if constexpr (Operation == Op::MatrixVector)
              right = referenceVector(input.vector);
            else if constexpr (Operation == Op::VectorMatrix)
            {
              left = referenceVector(input.rowVector).transpose();
              right = a;
            }
            result = Reference::Zero(left.rows(), right.cols());
            for (Eigen::Index i = 0; i < left.rows(); ++i)
            {
              for (Eigen::Index j = 0; j < right.cols(); ++j)
              {
                for (Eigen::Index k = 0; k < left.cols(); ++k)
                  result(i, j) += left(i, k) * right(k, j);
              }
            }
          }
          else if constexpr (Operation == Op::Cross)
          {
            for (Eigen::Index i = 0; i < 3; ++i)
            {
              result(i, 0) = a((i + 1) % 3, 0) * b((i + 2) % 3, 0) -
                a((i + 2) % 3, 0) * b((i + 1) % 3, 0);
            }
          }
          else if constexpr (Operation == Op::Transpose || Operation == Op::Adjoint)
          {
            result.resize(a.cols(), a.rows());
            for (Eigen::Index i = 0; i < a.rows(); ++i)
            {
              for (Eigen::Index j = 0; j < a.cols(); ++j)
                result(j, i) = Operation == Op::Adjoint ? Math::conj(a(i, j)) : a(i, j);
            }
          }
          else if constexpr (Operation == Op::Row)
          {
            result.resize(a.cols(), 1);
            for (Eigen::Index j = 0; j < a.cols(); ++j)
              result(j, 0) = a(0, j);
          }
          else if constexpr (Operation == Op::Column)
            result = a.col(0);
          else
          {
            Real norm = 0;
            if constexpr (Operation == Op::Normalized || Operation == Op::Normalize)
            {
              for (Eigen::Index i = 0; i < a.rows(); ++i)
              {
                for (Eigen::Index j = 0; j < a.cols(); ++j)
                  norm += std::norm(a(i, j));
              }
              norm = std::sqrt(norm);
            }
            for (Eigen::Index i = 0; i < a.rows(); ++i)
            {
              for (Eigen::Index j = 0; j < a.cols(); ++j)
              {
                if constexpr (Operation == Op::Zero)
                  result(i, j) = Scalar(0);
                else if constexpr (Operation == Op::MathSum || Operation == Op::Add ||
                  Operation == Op::AddAssign)
                  result(i, j) = a(i, j) + b(i, j);
                else if constexpr (Operation == Op::MathSubtract ||
                  Operation == Op::Subtract || Operation == Op::SubtractAssign)
                  result(i, j) = a(i, j) - b(i, j);
                else if constexpr (Operation == Op::Negate || Operation == Op::MathNegate)
                  result(i, j) = -a(i, j);
                else if constexpr (Operation == Op::MathScale ||
                  Operation == Op::ScaleLeft || Operation == Op::ScaleRight ||
                  Operation == Op::ScaleAssign)
                  result(i, j) = a(i, j) * factor;
                else if constexpr (Operation == Op::MathDivide ||
                  Operation == Op::Divide || Operation == Op::DivideAssign)
                  result(i, j) = a(i, j) / factor;
                else if constexpr (Operation == Op::Normalized ||
                  Operation == Op::Normalize)
                  result(i, j) = a(i, j) / norm;
                else if constexpr (Operation == Op::Conjugate ||
                  Operation == Op::MathConjugate)
                  result(i, j) = Math::conj(a(i, j));
                else
                  result(i, j) = a(i, j);
              }
            }
          }
          return result;
        }
      }

      /**
       * @brief Checks dimensions, finite values and roundoff-scaled agreement.
       *
       * The relative budget of 1024 machine epsilons is a fixture policy for
       * independently ordered reductions and factorizations on bounded,
       * diagonally strengthened inputs of dimension at most three. It is
       * neither a general conditioning bound nor a performance threshold.
       *
       * @tparam Actual Measured result type.
       * @tparam Expected Independent oracle result type.
       * @param actual Computed result.
       * @param expected Reference result.
       * @returns Whether both results agree within the small-operation budget.
       */
      template <class Actual, class Expected>
      static bool close(const Actual& actual, const Expected& expected)
      {
        constexpr Real tolerance = 1024 * std::numeric_limits<Real>::epsilon();
        if constexpr (std::is_arithmetic_v<Expected> || std::is_same_v<Expected, Complex>)
          return std::abs(actual - expected) <= tolerance * (1 + std::abs(expected));
        else
        {
          const auto data = [&] {
            if constexpr (requires { actual.getData(); })
            {
              if constexpr (requires { actual.rows(); })
              {
                return Reference(
                  actual.getData().topLeftCorner(actual.rows(), actual.cols()));
              }
              else
                return Reference(actual.getData().head(actual.size()));
            }
            else
              return Reference(actual);
          }();
          return data.rows() == expected.rows() && data.cols() == expected.cols() &&
            (data - expected).norm() <= tolerance * (1 + expected.norm());
        }
      }
  };
}

#endif
