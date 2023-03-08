/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FROBENIUS_H
#define RODIN_VARIATIONAL_FROBENIUS_H

#include "ForwardDecls.h"
#include "Function.h"
#include "ScalarFunction.h"

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
    : public ScalarFunctionBase<Frobenius<FunctionBase<NestedDerived>>>
  {
    public:
      using Operand = FunctionBase<NestedDerived>;
      using Parent = ScalarFunctionBase<Frobenius<Operand>>;

      Frobenius(const Operand& v)
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

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        using OperandRange = typename FormLanguage::Traits<Operand>::RangeType;
        if constexpr (std::is_same_v<OperandRange, Scalar>)
        {
          return std::abs(getOperand().getValue(p));
        }
        else
        {
          return getOperand().getValue(p).norm();
        }
      }

      const Operand& getOperand() const
      {
        assert(m_v);
        return *m_v;
      }

      inline Frobenius* copy() const noexcept override
      {
        return new Frobenius(*this);
      }

    private:
      std::unique_ptr<Operand> m_v;
  };

  template <class NestedDerived>
  Frobenius(const FunctionBase<NestedDerived>&) -> Frobenius<FunctionBase<NestedDerived>>;
}

#endif

