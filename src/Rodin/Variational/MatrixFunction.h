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
#include "RangeShape.h"

namespace Rodin::Variational
{
  /**
   * @defgroup MatrixFunctionSpecializations MatrixFunction Template Specializations
   * @brief Template specializations of the MatrixFunction class.
   * @see MatrixFunction
   */

  template <class Derived>
  class MatrixFunctionBase : public FunctionBase<MatrixFunctionBase<Derived>>
  {
    public:
      using Parent = FunctionBase<MatrixFunctionBase<Derived>>;

      MatrixFunctionBase() = default;

      MatrixFunctionBase(const MatrixFunctionBase& other)
        : Parent(other)
      {}

      MatrixFunctionBase(MatrixFunctionBase&& other)
        : Parent(std::move(other))
      {}

      virtual ~MatrixFunctionBase() = default;

      inline
      const Derived& getDerived() const
      {
        return static_cast<const Derived&>(*this);
      }

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
      }

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return { getRows(), getColumns() };
      }

      /**
       * @brief Gets the number of rows in the matrix
       * @returns Number of rows
       */
      inline
      constexpr
      size_t getRows() const
      {
        return static_cast<const Derived&>(*this).getRows();
      }

      /**
       * @brief Gets the number of columns in the matrix
       * @returns Number of columns
       */
      inline
      constexpr
      size_t getColumns() const
      {
        return static_cast<const Derived&>(*this).getColumns();
      }

      inline
      constexpr
      MatrixFunctionBase& traceOf(Geometry::Attribute attr)
      {
        Parent::traceOf(attr);
        return *this;
      }

      inline
      constexpr
      MatrixFunctionBase& traceOf(const FlatSet<Geometry::Attribute>& attrs)
      {
        Parent::traceOf(attrs);
        return *this;
      }

      virtual inline MatrixFunctionBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }
  };

  /**
   * @ingroup MatrixFunctionSpecializations
   */
  template <>
  class MatrixFunction<Math::Matrix<Scalar>> final
    : public MatrixFunctionBase<MatrixFunction<Math::Matrix<Scalar>>>
  {
    public:
      using Parent = MatrixFunctionBase<MatrixFunction<Math::Matrix<Scalar>>>;

      MatrixFunction(std::reference_wrapper<const Math::Matrix<Scalar>> matrix)
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

      inline
      const Math::Matrix<Scalar>& getValue(const Geometry::Point&) const
      {
        return m_matrix.get();
      }

      inline
      constexpr
      MatrixFunction& traceOf(Geometry::Attribute)
      {
        return *this;
      }

      inline
      constexpr
      MatrixFunction& traceOf(const FlatSet<Geometry::Attribute>& attr)
      {
        return *this;
      }

      inline
      constexpr
      size_t getRows() const
      {
        return m_matrix.get().rows();
      }

      /**
       * @brief Gets the number of columns in the matrix
       * @returns Number of columns
       */
      inline
      constexpr
      size_t getColumns() const
      {
        return m_matrix.get().cols();
      }

      inline MatrixFunction* copy() const noexcept override
      {
        return new MatrixFunction(*this);
      }

    private:
      std::reference_wrapper<const Math::Matrix<Scalar>> m_matrix;
  };

  MatrixFunction(std::reference_wrapper<const Math::Matrix<Scalar>>)
    -> MatrixFunction<Math::Matrix<Scalar>>;
}

#endif
