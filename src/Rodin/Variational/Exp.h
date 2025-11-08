/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_EXP_H
#define RODIN_VARIATIONAL_EXP_H

#include "Rodin/Math/Common.h"

#include "ForwardDecls.h"
#include "Function.h"
#include "RealFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup ExpSpecializations Exp Template Specializations
   * @brief Template specializations of the Exp class.
   * @see Exp
   */

  /**
   * @ingroup ExpSpecializations
   */
  template <class NestedDerived>
  class Exp<FunctionBase<NestedDerived>>
    : public RealFunctionBase<Exp<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using Parent = RealFunctionBase<Exp<OperandType>>;

      using OperandRangeType = typename FormLanguage::Traits<OperandType>::RangeType;

      Exp(const OperandType& v)
        : m_v(v.copy())
      {}

      Exp(const Exp& other)
        : Parent(other),
          m_v(other.m_v->copy())
      {}

      Exp(Exp&& other)
        : Parent(std::move(other)),
          m_v(std::move(other.m_v))
      {}

      template <class ... Args>
      constexpr
      Exp& traceOf(const Args& ... args)
      {
        m_v->traceOf(args...);
        return *this;
      }

      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return Math::exp(this->getOperand().getValue(p));
      }

      const OperandType& getOperand() const
      {
        assert(m_v);
        return *m_v;
      }

      Exp* copy() const noexcept override
      {
        return new Exp(*this);
      }

    private:
      std::unique_ptr<OperandType> m_v;
  };

  template <class NestedDerived>
  Exp(const FunctionBase<NestedDerived>&) -> Exp<FunctionBase<NestedDerived>>;

  template <class NestedDerived>
  constexpr auto
  exp(const FunctionBase<NestedDerived>& op)
  {
    return Exp(op);
  }
}

#endif



