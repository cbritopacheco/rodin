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
      using Operand = FunctionBase<FunctionBase<NestedDerived>>;
      using Parent = ScalarFunctionBase<Frobenius<FunctionBase<NestedDerived>>>;

      constexpr
      Frobenius(const Operand& v)
        : m_v(v)
      {}

      constexpr
      Frobenius(const Frobenius& other)
        : Parent(other),
          m_v(other.m_v)
      {}

      constexpr
      Frobenius(Frobenius&& other)
        : Parent(std::move(other)),
          m_v(std::move(other.m_v))
      {}

      inline
      constexpr
      Frobenius& traceOf(Geometry::Attribute attrs)
      {
        m_v.traceOf(attrs);
        return *this;
      }

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        using OperandRange = FormLanguage::RangeOf<typename FormLanguage::Traits<Operand>::ResultType>;
        if constexpr (std::is_same_v<OperandRange, Scalar>)
        {
          return std::abs(m_v.getValue(p));
        }
        else
        {
          return m_v.getValue(p).norm();
        }
      }

    private:
      Operand m_v;
  };

  template <class NestedDerived>
  Frobenius(const FunctionBase<NestedDerived>&) -> Frobenius<FunctionBase<NestedDerived>>;
}

#endif

