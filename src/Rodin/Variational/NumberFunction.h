/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_NUMBERFUNCTION_H
#define RODIN_VARIATIONAL_NUMBERFUNCTION_H

#include <map>
#include <set>
#include <memory>
#include <optional>
#include <type_traits>

#include "Rodin/Geometry/Polytope.h"

#include "ForwardDecls.h"

#include "Function.h"
#include "RangeShape.h"

namespace Rodin::Variational
{
  /**
   * @defgroup NumberFunctionSpecializations NumberFunction Template Specializations
   * @brief Template specializations of the NumberFunction class.
   * @see NumberFunction
   */

  template <class Number, class Derived>
  class NumberFunctionBase : public FunctionBase<NumberFunctionBase<Number, Derived>>
  {
    public:
      using NumberType = Number;

      using RangeType = NumberType;

      using Parent = FunctionBase<NumberFunctionBase<NumberType, Derived>>;

      using Parent::traceOf;

      NumberFunctionBase() = default;

      NumberFunctionBase(const NumberFunctionBase& other)
        : Parent(other)
      {}

      NumberFunctionBase(NumberFunctionBase&& other)
        : Parent(std::move(other))
      {}

      virtual ~NumberFunctionBase() = default;

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
      void getValue(Math::Vector<Number>&, const Geometry::Point&) const = delete;

      inline
      constexpr
      void getValue(Math::Matrix<Number>&, const Geometry::Point&) const = delete;

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return { 1, 1 };
      }

      virtual NumberFunctionBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }
  };
}

#endif

