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

#include "ForwardDecls.h"
#include "Function.h"

namespace Rodin::FormLanguage
{
  template <class Scalar, class Derived>
  struct Traits<Variational::MatrixFunctionBase<Scalar, Derived>>
  {
    using ScalarType = Scalar;
    using DerivedType = Derived;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup MatrixFunctionSpecializations MatrixFunction Template Specializations
   * @brief Template specializations of the MatrixFunction class.
   * @see MatrixFunction
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
   * @see FunctionBase, VectorFunction, Transpose
   */
  template <class Scalar, class Derived>
  class MatrixFunctionBase : public FunctionBase<MatrixFunctionBase<Scalar, Derived>>
  {
    public:
      /// @brief Type of scalar entries
      using ScalarType = Scalar;

      /// @brief Parent class type
      using Parent = FunctionBase<MatrixFunctionBase<ScalarType, Derived>>;

      /// @brief Import traceOf methods from parent
      using Parent::traceOf;

      /// @brief Import operator() from parent
      using Parent::operator();

      /// @brief Default constructor
      MatrixFunctionBase() = default;

      /// @brief Copy constructor
      /// @param[in] other Matrix function to copy from
      MatrixFunctionBase(const MatrixFunctionBase& other)
        : Parent(other)
      {}

      /// @brief Move constructor
      /// @param[in] other Matrix function to move from
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
      decltype(auto) getValue(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
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
   * @ingroup MatrixFunctionSpecializations
   */
  template <class Scalar>
  class MatrixFunction<Math::Matrix<Scalar>> final
    : public MatrixFunctionBase<Scalar, MatrixFunction<Math::Matrix<Scalar>>>
  {
    public:
      using ScalarType = Scalar;

      using MatrixType = Math::Matrix<ScalarType>;

      using Parent = MatrixFunctionBase<Scalar, MatrixFunction<MatrixType>>;

      using Parent::traceOf;

      MatrixFunction(const MatrixType& matrix)
        : m_matrix(matrix)
      {}

      MatrixFunction(const MatrixFunction& other)
        : Parent(other),
          m_matrix(other.m_matrix)
      {}

      MatrixFunction(MatrixFunction&& other)
        : Parent(std::move(other)),
          m_matrix(std::move(other.m_matrix))
      {}

      constexpr
      const auto& getValue(const Geometry::Point&) const
      {
        return m_matrix;
      }

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

  template <class Scalar>
  MatrixFunction(const Math::Matrix<Scalar>&)
    -> MatrixFunction<Math::Matrix<Scalar>>;
}

#endif
