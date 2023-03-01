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
#include <mfem.hpp>

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

      constexpr
      MatrixFunctionBase() = default;

      constexpr
      MatrixFunctionBase(const MatrixFunctionBase& other)
        : Parent(other)
      {}

      constexpr
      MatrixFunctionBase(MatrixFunctionBase&& other)
        : Parent(std::move(other))
      {}

      virtual ~MatrixFunctionBase() = default;

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
      MatrixFunctionBase& traceOf(const std::set<Geometry::Attribute>& attrs)
      {
        Parent::traceOf(attrs);
        return *this;
      }

      virtual MatrixFunctionBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }
  };
}

#endif
