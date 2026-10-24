/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Matrix.h
 * @brief Dense matrix classes and scalar aliases.
 *
 * This file provides Rodin matrix types with Eigen storage and expression support.
 * Matrices are stored in column-major order by default and support both dynamic
 * and fixed sizes.
 */
#ifndef RODIN_MATH_MATRIX_H
#define RODIN_MATH_MATRIX_H

#include <algorithm>
#include <initializer_list>
#include <type_traits>
#include <limits>
#include <utility>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/src/Core/ArrayBase.h>

#include "Rodin/Configure.h"

#include "Rodin/FormLanguage/Traits.h"

#include "Rodin/Math/Traits.h"
#include "Rodin/Types.h"

#include "ForwardDecls.h"
#include "Common.h"

namespace Eigen::internal
{
  /// @brief Retains Eigen map and reference support for dynamic Rodin matrices.
  template <class Scalar>
  struct traits<Rodin::Math::Matrix<Scalar>>
    : traits<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>
  {};

  /// @brief Retains Eigen map and reference support for fixed-size Rodin matrices.
  template <class Scalar, size_t Rows, size_t Cols>
  struct traits<Rodin::Math::FixedSizeMatrix<Scalar, Rows, Cols>>
    : traits<Eigen::Matrix<Scalar, Rows, Cols>>
  {};

  /// @brief Retains an already evaluated dynamic matrix without copying it.
  template <class Scalar>
  struct eval<Rodin::Math::Matrix<Scalar>, Dense>
  {
    /// @brief Reference to the existing matrix coefficients.
    using type = const Rodin::Math::Matrix<Scalar>&;
  };

  /// @brief Retains an already evaluated fixed-size matrix without copying it.
  template <class Scalar, size_t Rows, size_t Cols>
  struct eval<Rodin::Math::FixedSizeMatrix<Scalar, Rows, Cols>, Dense>
  {
    /// @brief Reference to the existing matrix coefficients.
    using type = const Rodin::Math::FixedSizeMatrix<Scalar, Rows, Cols>&;
  };

  /// @brief Evaluates dynamic Rodin matrices through their Eigen storage.
  template <class Scalar>
  struct evaluator<Rodin::Math::Matrix<Scalar>>
    : evaluator<Eigen::PlainObjectBase<Rodin::Math::Matrix<Scalar>>>
  {
    /// @brief Matrix whose coefficients are evaluated.
    using XprType = Rodin::Math::Matrix<Scalar>;

    /// @brief Native Eigen storage evaluator.
    using Parent = evaluator<Eigen::PlainObjectBase<XprType>>;

    /// @brief Inherits evaluation of the Eigen storage.
    using Parent::Parent;
  };

  /// @brief Evaluates fixed-size Rodin matrices through their Eigen storage.
  template <class Scalar, size_t Rows, size_t Cols>
  struct evaluator<Rodin::Math::FixedSizeMatrix<Scalar, Rows, Cols>>
    : evaluator<Eigen::PlainObjectBase<Rodin::Math::FixedSizeMatrix<Scalar, Rows, Cols>>>
  {
    /// @brief Matrix whose coefficients are evaluated.
    using XprType = Rodin::Math::FixedSizeMatrix<Scalar, Rows, Cols>;

    /// @brief Native Eigen storage evaluator.
    using Parent = evaluator<Eigen::PlainObjectBase<XprType>>;

    /// @brief Inherits evaluation of the Eigen storage.
    using Parent::Parent;
  };
}

namespace Rodin::Math
{
  /**
   * @brief Dynamic-size dense matrix type.
   *
   * A dense matrix with both dimensions determined at runtime. The matrix
   * is stored in column-major order (Fortran style) and supports standard
   * linear algebra operations without additional storage.
   * Dot products conjugate the right operand.
   *
   * @tparam ScalarType The element type (e.g., Real, Complex)
   */
  template <class ScalarType>
  class Matrix final : public Eigen::PlainObjectBase<Matrix<ScalarType>>
  {
    public:
      /// @brief Eigen storage and expression base.
      using Parent = Eigen::PlainObjectBase<Matrix>;

      /// @brief Base alias used by Eigen's storage operations.
      using Base = Parent;

      /// @brief Scalar coefficient type.
      using Scalar = ScalarType;

      /// @brief Reference retained when nested in an Eigen expression.
      using Nested = const Matrix&;

      /// @brief Assigns Eigen expressions into the owned storage.
      using Parent::operator=;

      /// @brief Constructs a matrix with Eigen's default initialization.
      Matrix()
      {
        if constexpr (requires { Parent::_check_template_params(); })
        {
          Parent::_check_template_params();
        }
#if defined(EIGEN_INITIALIZE_MATRICES_BY_ZERO)
        this->setZero();
#elif defined(EIGEN_INITIALIZE_MATRICES_BY_NAN)
        this->setConstant(std::numeric_limits<Scalar>::quiet_NaN());
#endif
      }

      /**
       * @brief Copies the matrix and its owned coefficients.
       * @param other Matrix to copy.
       */
      Matrix(const Matrix& other) = default;

      /**
       * @brief Moves the matrix and its owned coefficients.
       * @param other Matrix to move.
       */
      Matrix(Matrix&& other) noexcept = default;

      /**
       * @brief Constructs a matrix from dimensions, coefficients or an Eigen object.
       * @tparam T Argument type understood by Eigen's matrix constructor.
       * @param value Dimension, coefficient, coefficient array or Eigen object.
       */
      template <class T>
      explicit Matrix(const T& value)
      {
        if constexpr (requires { Parent::_check_template_params(); })
        {
          Parent::_check_template_params();
        }
        Parent::template _init1<T>(value);
      }

      /**
       * @brief Constructs a matrix from dimensions or two vector coefficients.
       * @tparam T0 First argument type.
       * @tparam T1 Second argument type.
       * @param first Row count or first coefficient of a two-component vector.
       * @param second Column count or second coefficient of a two-component vector.
       */
      template <class T0, class T1>
      Matrix(const T0& first, const T1& second)
      {
        if constexpr (requires { Parent::_check_template_params(); })
        {
          Parent::_check_template_params();
        }
        Parent::template _init2<T0, T1>(first, second);
      }

      /**
       * @brief Copies the coefficients of an Eigen dense expression.
       * @tparam Derived Eigen expression type.
       * @param other Expression to evaluate.
       */
      template <class Derived>
      Matrix(const Eigen::DenseBase<Derived>& other)
        : Parent(other)
      {}

      /**
       * @brief Evaluates an Eigen operator, such as a permutation or solve.
       * @tparam Derived Eigen operator type.
       * @param other Operator to evaluate.
       */
      template <class Derived>
      Matrix(const Eigen::EigenBase<Derived>& other)
        : Parent(other)
      {}

      /**
       * @brief Evaluates an Eigen return-by-value operation.
       * @tparam Derived Eigen operation type.
       * @param other Operation to evaluate.
       */
      template <class Derived>
      Matrix(const Eigen::ReturnByValue<Derived>& other)
        : Parent(other)
      {}

      /**
       * @brief Constructs a matrix from coefficients grouped by row.
       * @param rows Coefficient rows.
       */
      explicit Matrix(const std::initializer_list<std::initializer_list<Scalar>>& rows)
        : Parent(rows)
      {}

      /**
       * @brief Takes ownership of a native Eigen matrix without reallocating.
       * @param other Native matrix to move.
       */
      Matrix(Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>&& other) noexcept
        : Matrix()
      {
        Eigen::internal::matrix_swap_impl<Matrix,
          Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>, true>::run(*this, other);
      }

      /**
       * @brief Exchanges owned storage with a native Eigen matrix.
       * @param other Native matrix to move.
       * @returns This matrix.
       */
      Matrix& operator=(
        Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>&& other) noexcept
      {
        Eigen::internal::matrix_swap_impl<Matrix,
          Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>, true>::run(*this, other);
        return *this;
      }

      /**
       * @brief Copies the matrix coefficients.
       * @param other Matrix to copy.
       * @returns This matrix.
       */
      Matrix& operator=(const Matrix& other)
      {
        return Parent::_set(other);
      }

      /**
       * @brief Moves the matrix coefficients.
       * @param other Matrix to move.
       * @returns This matrix.
       */
      Matrix& operator=(Matrix&& other) noexcept = default;

      /**
       * @brief Assigns a dense expression, resizing when necessary.
       * @tparam Derived Eigen expression type.
       * @param other Expression to evaluate.
       * @returns This matrix.
       */
      template <class Derived>
      Matrix& operator=(const Eigen::DenseBase<Derived>& other)
      {
        return Parent::_set(other);
      }

      /**
       * @brief Assigns an Eigen operator, retaining the Rodin matrix return type.
       * @tparam Derived Eigen operator type.
       * @param other Operator to evaluate.
       * @returns This matrix.
       */
      template <class Derived>
      Matrix& operator=(const Eigen::EigenBase<Derived>& other)
      {
        Parent::operator=(other);
        return *this;
      }

      /**
       * @brief Assigns an Eigen return-by-value operation.
       * @tparam Derived Eigen operation type.
       * @param other Operation to evaluate.
       * @returns This matrix.
       */
      template <class Derived>
      Matrix& operator=(const Eigen::ReturnByValue<Derived>& other)
      {
        Parent::operator=(other);
        return *this;
      }

      /**
       * @brief Resizes while preserving the common top-left block.
       * @param rows New row count.
       * @param cols New column count.
       */
      void conservativeResize(Eigen::Index rows, Eigen::Index cols)
      {
        if (this->rows() == rows && this->cols() == cols)
        {
          return;
        }
        if (std::is_trivially_copyable_v<Scalar> &&
          (Parent::IsRowMajor ? this->cols() == cols : this->rows() == rows))
        {
#if EIGEN_VERSION_AT_LEAST(3, 4, 90)
#ifndef EIGEN_NO_DEBUG
          Eigen::internal::check_rows_cols_for_overflow<Parent::MaxSizeAtCompileTime,
            Parent::MaxRowsAtCompileTime, Parent::MaxColsAtCompileTime>::run(rows, cols);
#endif
#else
          Eigen::internal::check_rows_cols_for_overflow<
            Parent::MaxSizeAtCompileTime>::run(rows, cols);
#endif
          this->m_storage.conservativeResize(rows * cols, rows, cols);
        }
        else
        {
          Matrix result(rows, cols);
          const auto commonRows = std::min(rows, this->rows());
          const auto commonCols = std::min(cols, this->cols());
          result.topLeftCorner(commonRows, commonCols) =
            this->topLeftCorner(commonRows, commonCols);
          this->swap(result);
        }
      }

      /**
       * @brief Changes the row count while preserving existing coefficients.
       * @param rows New row count.
       * @param unchanged Retains the column count.
       */
      void conservativeResize(Eigen::Index rows, Eigen::NoChange_t unchanged)
      {
        (void)unchanged;
        conservativeResize(rows, this->cols());
      }

      /**
       * @brief Changes the column count while preserving existing coefficients.
       * @param unchanged Retains the row count.
       * @param cols New column count.
       */
      void conservativeResize(Eigen::NoChange_t unchanged, Eigen::Index cols)
      {
        (void)unchanged;
        conservativeResize(this->rows(), cols);
      }

      /**
       * @brief Resizes to an expression, retaining old values and copying new entries.
       * @tparam Derived Eigen dense expression type.
       * @param other Expression supplying the new shape and appended coefficients.
       */
      template <class Derived>
      void conservativeResizeLike(const Eigen::DenseBase<Derived>& other)
      {
        if (this->rows() == other.rows() && this->cols() == other.cols())
        {
          return;
        }
        if (std::is_trivially_copyable_v<Scalar> &&
          (Parent::IsRowMajor ? this->cols() == other.cols()
                              : this->rows() == other.rows()))
        {
#if EIGEN_VERSION_AT_LEAST(3, 4, 90)
#ifndef EIGEN_NO_DEBUG
          Eigen::internal::check_rows_cols_for_overflow<Parent::MaxSizeAtCompileTime,
            Parent::MaxRowsAtCompileTime, Parent::MaxColsAtCompileTime>::run(other.rows(),
            other.cols());
#endif
#else
          Eigen::internal::check_rows_cols_for_overflow<
            Parent::MaxSizeAtCompileTime>::run(other.rows(), other.cols());
#endif
          const auto addedRows = other.rows() - this->rows();
          const auto addedCols = other.cols() - this->cols();
          this->m_storage.conservativeResize(other.size(), other.rows(), other.cols());
          if (addedRows > 0)
          {
            this->bottomRows(addedRows) = other.bottomRows(addedRows);
          }
          else if (addedCols > 0)
          {
            this->rightCols(addedCols) = other.rightCols(addedCols);
          }
        }
        else
        {
          Matrix result(other);
          const auto commonRows = std::min(result.rows(), this->rows());
          const auto commonCols = std::min(result.cols(), this->cols());
          result.topLeftCorner(commonRows, commonCols) =
            this->topLeftCorner(commonRows, commonCols);
          this->swap(result);
        }
      }

      /**
       * @brief Stride between consecutive coefficients in the inner dimension.
       * @returns One, since inner coefficients are contiguous.
       */
      constexpr Eigen::Index innerStride() const noexcept
      {
        return 1;
      }

      /**
       * @brief Stride between consecutive columns, or rows for row-major storage.
       * @returns The number of coefficients in the inner dimension.
       */
      constexpr Eigen::Index outerStride() const noexcept
      {
        return this->innerSize();
      }

      /**
       * @brief Contracts with the conjugate of the right operand.
       *
       * Computes @f$\sum_{i,j} A_{ij}\overline{B_{ij}}@f$, consistently with
       * @ref Rodin::Math::dot. Eigen expressions and views retain their Eigen
       * types; use @ref Rodin::Math::dot to contract those with this convention.
       *
       * @tparam Derived Eigen type of the right operand.
       * @param other Matrix of the same shape.
       * @returns The right-conjugating inner product.
       */
      template <class Derived>
      auto dot(const Eigen::MatrixBase<Derived>& other) const
      {
        return Math::dot(*this, other);
      }
  };

  /**
   * @brief Real-valued dense matrix.
   *
   * Convenience alias for Matrix<Real>.
   */
  using RealMatrix = Matrix<Real>;

  /**
   * @brief Complex-valued dense matrix.
   *
   * Convenience alias for Matrix<Complex>, commonly used in eigenvalue
   * problems, frequency-domain analysis, and other applications requiring
   * complex arithmetic.
   */
  using ComplexMatrix = Matrix<Complex>;

  /**
   * @brief Fixed-size matrix type.
   *
   * A compile-time fixed-size matrix. Both dimensions are known at compile time,
   * enabling optimizations and stack allocation without additional storage.
   * Dot products conjugate the right operand.
   *
   * Example uses:
   * - 2×2 rotation matrices
   * - 3×3 transformation matrices
   * - Small dense blocks in larger systems
   *
   * @tparam ScalarType The element type
   * @tparam Rows The number of rows (must be known at compile time)
   * @tparam Cols The number of columns (must be known at compile time)
   */
  template <class ScalarType, size_t Rows, size_t Cols>
  class FixedSizeMatrix final
    : public Eigen::PlainObjectBase<FixedSizeMatrix<ScalarType, Rows, Cols>>
  {
    public:
      /// @brief Eigen storage and expression base.
      using Parent = Eigen::PlainObjectBase<FixedSizeMatrix>;

      /// @brief Base alias used by Eigen's storage operations.
      using Base = Parent;

      /// @brief Scalar coefficient type.
      using Scalar = ScalarType;

      /// @brief Reference retained when nested in an Eigen expression.
      using Nested = const FixedSizeMatrix&;

      /// @brief Assigns Eigen expressions into the owned storage.
      using Parent::operator=;

      /// @brief Constructs a matrix with Eigen's default initialization.
      FixedSizeMatrix()
      {
        if constexpr (requires { Parent::_check_template_params(); })
        {
          Parent::_check_template_params();
        }
#if defined(EIGEN_INITIALIZE_MATRICES_BY_ZERO)
        this->setZero();
#elif defined(EIGEN_INITIALIZE_MATRICES_BY_NAN)
        this->setConstant(std::numeric_limits<Scalar>::quiet_NaN());
#endif
      }

      /**
       * @brief Copies the matrix and its owned coefficients.
       * @param other Matrix to copy.
       */
      FixedSizeMatrix(const FixedSizeMatrix& other) = default;

      /**
       * @brief Moves the matrix and its owned coefficients.
       * @param other Matrix to move.
       */
      FixedSizeMatrix(FixedSizeMatrix&& other) noexcept = default;

      /**
       * @brief Constructs a matrix from dimensions, coefficients or an Eigen object.
       * @tparam T Argument type understood by Eigen's matrix constructor.
       * @param value Dimension, coefficient, coefficient array or Eigen object.
       */
      template <class T>
      explicit FixedSizeMatrix(const T& value)
      {
        if constexpr (requires { Parent::_check_template_params(); })
        {
          Parent::_check_template_params();
        }
        Parent::template _init1<T>(value);
      }

      /**
       * @brief Constructs a matrix from dimensions or two vector coefficients.
       * @tparam T0 First argument type.
       * @tparam T1 Second argument type.
       * @param first Row count or first coefficient of a two-component vector.
       * @param second Column count or second coefficient of a two-component vector.
       */
      template <class T0, class T1>
      FixedSizeMatrix(const T0& first, const T1& second)
      {
        if constexpr (requires { Parent::_check_template_params(); })
        {
          Parent::_check_template_params();
        }
        Parent::template _init2<T0, T1>(first, second);
      }

      /**
       * @brief Copies the coefficients of an Eigen dense expression.
       * @tparam Derived Eigen expression type.
       * @param other Expression to evaluate.
       */
      template <class Derived>
      FixedSizeMatrix(const Eigen::DenseBase<Derived>& other)
        : Parent(other)
      {}

      /**
       * @brief Evaluates an Eigen operator, such as a permutation or solve.
       * @tparam Derived Eigen operator type.
       * @param other Operator to evaluate.
       */
      template <class Derived>
      FixedSizeMatrix(const Eigen::EigenBase<Derived>& other)
        : Parent(other)
      {}

      /**
       * @brief Evaluates an Eigen return-by-value operation.
       * @tparam Derived Eigen operation type.
       * @param other Operation to evaluate.
       */
      template <class Derived>
      FixedSizeMatrix(const Eigen::ReturnByValue<Derived>& other)
        : Parent(other)
      {}

      /**
       * @brief Constructs a matrix from coefficients grouped by row.
       * @param rows Coefficient rows.
       */
      explicit FixedSizeMatrix(
        const std::initializer_list<std::initializer_list<Scalar>>& rows)
        : Parent(rows)
      {}

      /**
       * @brief Constructs a three-component row or column vector.
       * @param first First coefficient.
       * @param second Second coefficient.
       * @param third Third coefficient.
       */
      FixedSizeMatrix(const Scalar& first, const Scalar& second, const Scalar& third)
      {
        static_assert((Rows == 1 || Cols == 1) && Rows * Cols == 3);
        this->coeffRef(0) = first;
        this->coeffRef(1) = second;
        this->coeffRef(2) = third;
      }

      /**
       * @brief Constructs a fixed-size vector from four or more coefficients.
       * @tparam Args Remaining coefficient types.
       * @param first First coefficient.
       * @param second Second coefficient.
       * @param third Third coefficient.
       * @param fourth Fourth coefficient.
       * @param args Remaining coefficients.
       */
      template <class... Args>
      FixedSizeMatrix(const Scalar& first, const Scalar& second, const Scalar& third,
        const Scalar& fourth, const Args&... args)
      {
        static_assert((Rows == 1 || Cols == 1) && Rows * Cols == sizeof...(Args) + 4);
        this->coeffRef(0) = first;
        this->coeffRef(1) = second;
        this->coeffRef(2) = third;
        this->coeffRef(3) = fourth;
        Eigen::Index index = 4;
        ((this->coeffRef(index++) = args), ...);
      }

      /**
       * @brief Copies the matrix coefficients.
       * @param other Matrix to copy.
       * @returns This matrix.
       */
      FixedSizeMatrix& operator=(const FixedSizeMatrix& other)
      {
        return Parent::_set(other);
      }

      /**
       * @brief Moves the matrix coefficients.
       * @param other Matrix to move.
       * @returns This matrix.
       */
      FixedSizeMatrix& operator=(FixedSizeMatrix&& other) noexcept = default;

      /**
       * @brief Assigns a dense expression, resizing when necessary.
       * @tparam Derived Eigen expression type.
       * @param other Expression to evaluate.
       * @returns This matrix.
       */
      template <class Derived>
      FixedSizeMatrix& operator=(const Eigen::DenseBase<Derived>& other)
      {
        return Parent::_set(other);
      }

      /**
       * @brief Assigns an Eigen operator, retaining the Rodin matrix return type.
       * @tparam Derived Eigen operator type.
       * @param other Operator to evaluate.
       * @returns This matrix.
       */
      template <class Derived>
      FixedSizeMatrix& operator=(const Eigen::EigenBase<Derived>& other)
      {
        Parent::operator=(other);
        return *this;
      }

      /**
       * @brief Assigns an Eigen return-by-value operation.
       * @tparam Derived Eigen operation type.
       * @param other Operation to evaluate.
       * @returns This matrix.
       */
      template <class Derived>
      FixedSizeMatrix& operator=(const Eigen::ReturnByValue<Derived>& other)
      {
        Parent::operator=(other);
        return *this;
      }

      /**
       * @brief Assigns an Eigen rotation through its matrix representation.
       * @tparam Derived Eigen rotation type.
       * @param rotation Rotation to evaluate.
       * @returns This matrix.
       */
      template <class Derived>
      FixedSizeMatrix& operator=(const Eigen::RotationBase<Derived, Cols>& rotation)
      {
        return *this = rotation.toRotationMatrix();
      }

      /**
       * @brief Stride between consecutive coefficients in the inner dimension.
       * @returns One, since inner coefficients are contiguous.
       */
      constexpr Eigen::Index innerStride() const noexcept
      {
        return 1;
      }

      /**
       * @brief Stride between consecutive columns, or rows for row-major storage.
       * @returns The number of coefficients in the inner dimension.
       */
      constexpr Eigen::Index outerStride() const noexcept
      {
        return this->innerSize();
      }

      /**
       * @brief Contracts with the conjugate of the right operand.
       *
       * Computes @f$\sum_{i,j} A_{ij}\overline{B_{ij}}@f$, consistently with
       * @ref Rodin::Math::dot. Eigen expressions and views retain their Eigen
       * types; use @ref Rodin::Math::dot to contract those with this convention.
       *
       * @tparam Derived Eigen type of the right operand.
       * @param other Matrix of the same shape.
       * @returns The right-conjugating inner product.
       */
      template <class Derived>
      auto dot(const Eigen::MatrixBase<Derived>& other) const
      {
        return Math::dot(*this, other);
      }
  };
}

namespace Rodin::FormLanguage
{
  /// @brief Scalar type of a native Eigen dense matrix, including fixed and bounded types.
  template <class Number, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
  struct Traits<Eigen::Matrix<Number, Rows, Cols, Options, MaxRows, MaxCols>>
  {
    /// @brief Scalar value type.
    using ScalarType = Number;
  };

  /// @brief Type traits for a dense Math::Matrix: exposes the scalar type.
  template <class Number>
  struct Traits<Math::Matrix<Number>>
  {
    /// @brief Scalar value type.
    using ScalarType = Number;
  };

  /// @brief Scalar type of a fixed-size Rodin dense matrix.
  template <class Number, size_t Rows, size_t Cols>
  struct Traits<Math::FixedSizeMatrix<Number, Rows, Cols>>
  {
    /// @brief Scalar value type.
    using ScalarType = Number;
  };
}

#endif
