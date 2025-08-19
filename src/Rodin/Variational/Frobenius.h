/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FROBENIUS_H
#define RODIN_VARIATIONAL_FROBENIUS_H

#include "Rodin/Math/Common.h"

#include "ForwardDecls.h"
#include "RealFunction.h"
#include "Function.h"

namespace Rodin::Variational
{
  /**
   * @defgroup FrobeniusSpecializations Frobenius Template Specializations
   * @brief Template specializations of the Frobenius class.
   * @see Frobenius
   */

  /**
   * @ingroup FrobeniusSpecializations
   */
  template <class NestedDerived>
  class Frobenius<FunctionBase<NestedDerived>>
    : public RealFunctionBase<Frobenius<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using OperandRangeType = typename FormLanguage::Traits<OperandType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<OperandType>::ScalarType;

      using Parent = RealFunctionBase<Frobenius<OperandType>>;

      Frobenius(const OperandType& v)
        : m_v(v.copy())
      {}

      Frobenius(const Frobenius& other)
        : Parent(other),
          m_v(other.m_v->copy())
      {}

      Frobenius(Frobenius&& other)
        : Parent(std::move(other)),
          m_v(std::move(other.m_v))
      {}

      decltype(auto) getValue(const Geometry::Point& p) const
      {
        if constexpr (std::is_same_v<OperandRangeType, Real>)
        {
          return Math::abs(this->getOperand().getValue(p));
        }
        else
        {
          static thread_local OperandRangeType s_v;
          s_v = this->getOperand().getValue(p);
          return s_v.norm();
        }
      }

      const OperandType& getOperand() const
      {
        assert(m_v);
        return *m_v;
      }

      Frobenius* copy() const noexcept override
      {
        return new Frobenius(*this);
      }

    private:
      std::unique_ptr<OperandType> m_v;
  };

  template <class NestedDerived>
  Frobenius(const FunctionBase<NestedDerived>&) -> Frobenius<FunctionBase<NestedDerived>>;
}

#endif

