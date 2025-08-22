/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
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

  template <class Scalar, class Derived>
  class MatrixFunctionBase : public FunctionBase<MatrixFunctionBase<Scalar, Derived>>
  {
    public:
      using ScalarType = Scalar;

      using Parent = FunctionBase<MatrixFunctionBase<ScalarType, Derived>>;

      using Parent::traceOf;

      using Parent::operator();

      MatrixFunctionBase() = default;

      MatrixFunctionBase(const MatrixFunctionBase& other)
        : Parent(other)
      {}

      MatrixFunctionBase(MatrixFunctionBase&& other)
        : Parent(std::move(other))
      {}

      virtual ~MatrixFunctionBase() = default;

      constexpr
      decltype(auto) getValue(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
      }

      template <class ... Args>
      constexpr
      Derived& traceOf(const Args& ... args)
      {
        return static_cast<Derived&>(*this).traceOf(args...);
      }

      /**
       * @brief Gets the number of rows in the matrix
       * @returns Number of rows
       */
      constexpr
      size_t getRows() const
      {
        return static_cast<const Derived&>(*this).getRows();
      }

      /**
       * @brief Gets the number of columns in the matrix
       * @returns Number of columns
       */
      constexpr
      size_t getColumns() const
      {
        return static_cast<const Derived&>(*this).getColumns();
      }

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
      decltype(auto) getValue(const Geometry::Point&) const
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
