/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_ABS_H
#define RODIN_VARIATIONAL_ABS_H

#include <cmath>
#include "ForwardDecls.h"
#include "Function.h"
#include "RealFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup AbsSpecializations Abs Template Specializations
   * @brief Template specializations of the Abs class.
   * @see Abs
   */

  /**
   * @ingroup AbsSpecializations
   */
  template <class NestedDerived>
  class Abs<FunctionBase<NestedDerived>>
    : public RealFunctionBase<Abs<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using Parent = RealFunctionBase<Abs<OperandType>>;

      using OperandRangeType = typename FormLanguage::Traits<OperandType>::RangeType;

      Abs(const OperandType& v)
        : m_v(v.copy())
      {}

      Abs(const Abs& other)
        : Parent(other),
          m_v(other.m_v->copy())
      {}

      Abs(Abs&& other)
        : Parent(std::move(other)),
          m_v(std::move(other.m_v))
      {}

      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return Math::abs(getOperand().getValue(p));
      }

      const OperandType& getOperand() const
      {
        assert(m_v);
        return *m_v;
      }

      Abs* copy() const noexcept override
      {
        return new Abs(*this);
      }

    private:
      std::unique_ptr<OperandType> m_v;
  };

  template <class NestedDerived>
  Abs(const FunctionBase<NestedDerived>&) -> Abs<FunctionBase<NestedDerived>>;

  template <class NestedDerived>
  constexpr auto
  abs(const FunctionBase<NestedDerived>& op)
  {
    return Abs(op);
  }
}

#endif


