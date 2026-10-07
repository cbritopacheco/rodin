/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file MatrixFunction.h
 * @brief Matrix-valued functions for variational formulations.
 *
 * This file defines MatrixFunctionBase and MatrixFunction for representing
 * functions mapping points to matrices: @f$ A: \Omega \to \mathbb{R}^{m \times n} @f$.
 * These are used for tensors, stress/strain fields, and coefficient matrices.
 */
#ifndef RODIN_VARIATIONAL_MATRIXFUNCTION_H
#define RODIN_VARIATIONAL_MATRIXFUNCTION_H

#include <set>
#include <optional>

#include "Rodin/Alert.h"
#include "Rodin/Math/SpatialMatrix.h"

#include "ForwardDecls.h"
#include "Function.h"

namespace Rodin::FormLanguage
{
  /**
   * @brief Type traits for @c MatrixFunctionBase: exposes the scalar type and the
   * derived type.
   */
  template <class Scalar, class Derived>
  struct Traits<Variational::MatrixFunctionBase<Scalar, Derived>>
  {
      /// @brief Scalar value type.
      using ScalarType = Scalar;
      /// @brief Derived CRTP function type.
      using DerivedType = Derived;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup MatrixFunctionSpecializations MatrixFunction Template Specializations
   * @brief Template specializations of the MatrixFunction class.
   * @see <a href="_matrix_function_8h.html">MatrixFunction</a>
   */

  /**
   * @brief Base class for matrix-valued functions.
   *
   * MatrixFunctionBase extends FunctionBase to represent matrix-valued functions:
   * @f[
   *    A: \Omega \to \mathbb{R}^{m \times n}
   * @f]
   * where @f$ m @f$ is the number of rows and @f$ n @f$ is the number of columns.
   *
   * These functions are used in finite element analysis for:
   * - **Material properties**: Diffusion tensors, conductivity matrices
   * - **Stress/strain**: Stress tensor @f$ \sigma(x) @f$, strain tensor @f$ \varepsilon(x) @f$
   * - **Gradients of vector fields**: @f$ \nabla \mathbf{u} @f$
   * - **Jacobians**: Transformation Jacobian matrices
   *
   * @tparam Scalar The scalar entry type (typically Real or Complex)
   * @tparam Derived The derived class following CRTP pattern
   *
   * ## Component Access
   * Matrix entries can be accessed via:
   * - `A(i, j)` for the entry at row i, column j
   * - `A(i)` for the i-th row (returning a vector function)
   *
   * @see <a href="class_rodin_1_1_variational_1_1_function_base.html">FunctionBase</a>
   * @see <a href="_vector_function_8h.html">VectorFunction</a>
   * @see <a href="_transpose_8h.html">Transpose</a>
   */
  template <class Scalar, class Derived>
  class MatrixFunctionBase : public FunctionBase<MatrixFunctionBase<Scalar, Derived>>
  {
    public:
      /// @brief Type of scalar entries
      using ScalarType = Scalar;

      /// @brief Parent class type
      using Parent = FunctionBase<MatrixFunctionBase<ScalarType, Derived>>;

      // Import traceOf methods from parent.
      using Parent::traceOf;

      // Import operator() from parent.
      using Parent::operator();

      /// @brief Default constructor
      MatrixFunctionBase() = default;

      /**
       * @brief Copy constructor
       * @param[in] other Matrix function to copy from
       */
      MatrixFunctionBase(const MatrixFunctionBase& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor
       * @param[in] other Matrix function to move from
       */
      MatrixFunctionBase(MatrixFunctionBase&& other)
        : Parent(std::move(other))
      {}

      /// @brief Virtual destructor
      virtual ~MatrixFunctionBase() = default;

      /**
       * @brief Evaluates the matrix function at a point.
       *
       * CRTP method delegating to derived class implementation.
       *
       * @param[in] p Point at which to evaluate
       * @returns Matrix value at the point
       */
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
      }

      /**
       * @brief Evaluates the expression at an integration point.
       * @param ip Integration point at which the expression is evaluated.
       * @returns Value of the expression at the supplied evaluation point.
       */
      constexpr
      auto getValue(const IntegrationPoint& ip) const
      {
        if constexpr (requires (const Derived& f, const IntegrationPoint& q) { f.getValue(q); })
          return static_cast<const Derived&>(*this).getValue(ip);
        else
          return static_cast<const Derived&>(*this).getValue(ip.getPoint());
      }

      /**
       * @brief Sets the trace domain for the function.
       *
       * @tparam Args Variadic template for trace domain specification
       * @param[in] args Arguments specifying the trace domain
       * @returns Reference to derived object (for method chaining)
       */
      template <class ... Args>
      constexpr
      Derived& traceOf(const Args& ... args)
      {
        return static_cast<Derived&>(*this).traceOf(args...);
      }

      /**
       * @brief Gets the number of rows in the matrix.
       * @returns Number of rows
       */
      constexpr
      size_t getRows() const
      {
        return static_cast<const Derived&>(*this).getRows();
      }

      /**
       * @brief Gets the number of columns in the matrix.
       * @returns Number of columns
       */
      constexpr
      size_t getColumns() const
      {
        return static_cast<const Derived&>(*this).getColumns();
      }

      /**
       * @brief Returns the polynomial order used on a mesh entity.
       * @param polytope Mesh entity used by this operation.
       * @returns Polynomial order on the entity, or an empty optional when no order is available.
       */
      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& polytope) const noexcept
      {
        return static_cast<const Derived&>(*this).getOrder(polytope);
      }

      /**
       * @brief Creates a polymorphic copy of the function.
       * @returns Pointer to newly allocated copy
       */
      virtual MatrixFunctionBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }
  };

  /**
   * @brief Matrix-valued constant function.
   * @ingroup MatrixFunctionSpecializations
   */
  template <class Scalar>
  class MatrixFunction<Math::Matrix<Scalar>> final
    : public MatrixFunctionBase<Scalar, MatrixFunction<Math::Matrix<Scalar>>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Scalar;

      /// @brief Matrix (operator) type of the linear system.
      using MatrixType = Math::Matrix<ScalarType>;

      /// @brief Parent class type.
      using Parent = MatrixFunctionBase<Scalar, MatrixFunction<MatrixType>>;

      using Parent::traceOf;

      /**
       * @brief Constructs the MatrixFunction from the given arguments.
       * @param matrix Matrix operand.
       */
      MatrixFunction(const MatrixType& matrix)
        : m_matrix(matrix)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      MatrixFunction(const MatrixFunction& other)
        : Parent(other),
          m_matrix(other.m_matrix)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      MatrixFunction(MatrixFunction&& other)
        : Parent(std::move(other)),
          m_matrix(std::move(other.m_matrix))
      {}

      /**
       * @brief Evaluates the expression at a geometric point.
       * @returns Value of the expression at the supplied evaluation point.
       */
      constexpr
      MatrixType getValue(const Geometry::Point&) const
      {
        return m_matrix;
      }

      /**
       * @brief Gets the number of rows.
       * @returns The number of rows.
       */
      constexpr
      size_t getRows() const
      {
        return m_matrix.rows();
      }

      /**
       * @brief Gets the number of columns in the matrix
       * @returns Number of columns
       */
      constexpr
      size_t getColumns() const
      {
        return m_matrix.cols();
      }

      /**
       * @brief Returns the polynomial order used on a mesh entity.
       * @returns Polynomial order on the entity, or an empty optional when no order is available.
       */
      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& ) const noexcept
      {
        return 0;
      }

      MatrixFunction* copy() const noexcept override
      {
        return new MatrixFunction(*this);
      }

    private:
      const MatrixType m_matrix;
  };

  /// @brief Deduction guide for @c MatrixFunction.
  template <class Scalar>
  MatrixFunction(const Math::Matrix<Scalar>&)
    -> MatrixFunction<Math::Matrix<Scalar>>;
}

namespace Rodin::Variational
{
  /// @brief Matrix-range finite element or expression specialization.
  template <class Scalar>
  class MatrixFunction<Math::SpatialMatrix<Scalar>> final
    : public MatrixFunctionBase<Scalar, MatrixFunction<Math::SpatialMatrix<Scalar>>>
  {
    public:
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = Scalar;

      /// @brief Fixed-capacity matrix value type.
      using MatrixType = Math::SpatialMatrix<ScalarType>;

      /// @brief CRTP or finite element base class.
      using Parent = MatrixFunctionBase<Scalar, MatrixFunction<MatrixType>>;

      using Parent::traceOf;

      /**
       * @brief Constructs a constant or callable matrix coefficient, or copies its value.
       * @param matrix Matrix operand.
       */
      MatrixFunction(const MatrixType& matrix)
        : m_matrix(matrix)
      {}

      /**
       * @brief Constructs a constant or callable matrix coefficient, or copies its value.
       * @param other Object to copy from.
       */
      MatrixFunction(const MatrixFunction& other)
        : Parent(other),
          m_matrix(other.m_matrix)
      {}

      /**
       * @brief Constructs a constant or callable matrix coefficient, or copies its value.
       * @param other Object to move from.
       */
      MatrixFunction(MatrixFunction&& other)
        : Parent(std::move(other)),
          m_matrix(std::move(other.m_matrix))
      {}

      /**
       * @brief Evaluates the expression at the supplied physical or integration point.
       * @returns Value of the expression at the supplied evaluation point.
       */
      constexpr MatrixType getValue(const Geometry::Point&) const
      {
        return m_matrix;
      }

      /**
       * @brief Returns the number of matrix rows.
       * @returns The number of matrix rows.
       */
      constexpr size_t getRows() const
      {
        return m_matrix.rows();
      }

      /**
       * @brief Gets the number of columns in the matrix
       * @returns Number of columns
       */
      constexpr size_t getColumns() const
      {
        return m_matrix.cols();
      }

      /**
       * @brief Returns the polynomial order when it is known.
       * @returns Polynomial order on the entity, or an empty optional when no order is available.
       */
      constexpr Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return 0;
      }

      MatrixFunction* copy() const noexcept override
      {
        return new MatrixFunction(*this);
      }

    private:
      const MatrixType m_matrix;
  };

  /// @brief Deduces the matrix space or coefficient type from constructor arguments.
  template <class Scalar>
  MatrixFunction(
    const Math::SpatialMatrix<Scalar>&) -> MatrixFunction<Math::SpatialMatrix<Scalar>>;
}

namespace Rodin::Variational
{
  /**
   * @ingroup MatrixFunctionSpecializations
   * @brief Callable matrix coefficient with explicit row and column extents.
   */
  template <class F>
  class MatrixFunction final
    : public MatrixFunctionBase<typename FormLanguage::Traits<std::invoke_result_t<F,
                                  const Geometry::Point&>>::ScalarType,
        MatrixFunction<F>>
  {
    public:
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename FormLanguage::Traits<
        std::invoke_result_t<F, const Geometry::Point&>>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialMatrix<ScalarType>;
      /// @brief CRTP or finite element base class.
      using Parent = MatrixFunctionBase<ScalarType, MatrixFunction>;
      /**
       * @brief Constructs a constant or callable matrix coefficient, or copies its value.
       * @param rows Number of rows.
       * @param function Function to evaluate.
       * @param columns Number of matrix columns.
       */
      MatrixFunction(size_t rows, size_t columns, F function)
        : m_rows(rows),
          m_columns(columns),
          m_function(std::move(function))
      {
        if (rows == 0 || columns == 0 || rows > RODIN_MAXIMAL_SPACE_DIMENSION ||
          columns > RODIN_MAXIMAL_SPACE_DIMENSION)
          Alert::Exception() << "Invalid spatial matrix coefficient dimensions."
                             << Alert::Raise;
      }
      /**
       * @brief Constructs a constant or callable matrix coefficient, or copies its value.
       * @param other Object to copy from.
       */
      MatrixFunction(const MatrixFunction& other)
        : Parent(other),
          m_rows(other.m_rows),
          m_columns(other.m_columns),
          m_function(other.m_function),
          m_order(other.m_order)
      {}
      /**
       * @brief Constructs a constant or callable matrix coefficient, or copies its value.
       * @param other Object to move from.
       */
      MatrixFunction(MatrixFunction&& other)
        : Parent(std::move(other)),
          m_rows(other.m_rows),
          m_columns(other.m_columns),
          m_function(std::move(other.m_function)),
          m_order(other.m_order)
      {}
      /**
       * @brief Evaluates the expression at the supplied physical or integration point.
       * @param point Point at which the operation is evaluated.
       * @returns Value of the expression at the supplied evaluation point.
       */
      RangeType getValue(const Geometry::Point& point) const
      {
        const auto result = m_function(point);
        if (result.rows() != m_rows || result.cols() != m_columns)
          Alert::Exception() << "Matrix coefficient returned incompatible dimensions."
                             << Alert::Raise;
        return RangeType(result);
      }
      /**
       * @brief Evaluates the expression at the supplied physical or integration point.
       * @param point Point at which the operation is evaluated.
       * @returns Value of the expression at the supplied evaluation point.
       */
      RangeType getValue(const IntegrationPoint& point) const
      {
        return getValue(point.getPoint());
      }
      /**
       * @brief Returns the number of matrix rows.
       * @returns The number of matrix rows.
       */
      size_t getRows() const
      {
        return m_rows;
      }
      /**
       * @brief Returns the number of matrix columns.
       * @returns The number of matrix columns.
       */
      size_t getColumns() const
      {
        return m_columns;
      }
      /**
       * @brief Declares polynomial order for coefficient quadrature selection.
       * @param order Polynomial order.
       * @returns Reference to this object after the operation.
       */
      MatrixFunction& setOrder(size_t order)
      {
        m_order = order;
        return *this;
      }
      /**
       * @brief Returns the polynomial order when it is known.
       * @returns Polynomial order on the entity, or an empty optional when no order is available.
       */
      Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return m_order;
      }
      MatrixFunction* copy() const noexcept override
      {
        return new MatrixFunction(*this);
      }

    private:
      size_t m_rows, m_columns;
      F m_function;
      Optional<size_t> m_order;
  };
  /// @brief Deduces the matrix space or coefficient type from constructor arguments.
  template <class F>
  MatrixFunction(size_t, size_t, F) -> MatrixFunction<F>;
}

#endif
