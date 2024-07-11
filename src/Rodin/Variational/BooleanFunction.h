/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_BOOLEANFUNCTION_H
#define RODIN_VARIATIONAL_BOOLEANFUNCTION_H

#include "ForwardDecls.h"
#include "Function.h"
#include "RangeShape.h"

namespace Rodin::Variational
{
  /**
   * @defgroup BooleanFunctionSpecializations BooleanFunction Template Specializations
   * @brief Template specializations of the BooleanFunction class.
   * @see BooleanFunction
   */

  template <class Derived>
  class BooleanFunctionBase
    : public FunctionBase<BooleanFunctionBase<Derived>>
  {
    public:
      using Parent = FunctionBase<BooleanFunctionBase<Derived>>;

      BooleanFunctionBase() = default;

      BooleanFunctionBase(const BooleanFunctionBase& other)
        : Parent(other)
      {}

      BooleanFunctionBase(BooleanFunctionBase&& other)
        : Parent(std::move(other))
      {}

      virtual ~BooleanFunctionBase() = default;

      inline
      constexpr
      BooleanFunctionBase& traceOf(Geometry::Attribute attr)
      {
        Parent::traceOf(attr);
        return *this;
      }

      inline
      constexpr
      BooleanFunctionBase& traceOf(const FlatSet<Geometry::Attribute>& attrs)
      {
        Parent::traceOf(attrs);
        return *this;
      }

      /**
       * @note CRTP function to be overriden in Derived class.
       */
      inline
      const Derived& getDerived() const
      {
        return static_cast<const Derived&>(*this);
      }

      /**
       * @note CRTP function to be overriden in Derived class.
       */
      inline
      constexpr
      Boolean getValue(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
      }

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return { 1, 1 };
      }

      virtual BooleanFunctionBase* copy() const noexcept override = 0;
  };

  /**
   * @ingroup BooleanFunctionSpecializations
   */
  template <>
  class BooleanFunction<Boolean> final
    : public BooleanFunctionBase<BooleanFunction<Boolean>>
  {
    public:
      using Parent = BooleanFunctionBase<BooleanFunction<Boolean>>;

      BooleanFunction(Boolean v)
        : m_v(v)
      {}

      BooleanFunction(const BooleanFunction& other)
        : Parent(other),
          m_v(other.m_v)
      {}

      BooleanFunction(BooleanFunction&& other)
        : Parent(std::move(other)),
          m_v(other.m_v)
      {}

      inline
      constexpr
      Boolean getValue(const Geometry::Point&) const
      {
        return m_v;
      }

      inline BooleanFunction* copy() const noexcept final override
      {
        return new BooleanFunction(*this);
      }

    private:
      const Boolean m_v;
  };

  BooleanFunction(Boolean) -> BooleanFunction<Boolean>;
}

#endif

