/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_Cos_H
#define RODIN_VARIATIONAL_Cos_H

#include "Rodin/Math.h"
#include "ForwardDecls.h"
#include "Function.h"
#include "ScalarFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup CosSpecializations Cos Template Specializations
   * @brief Template specializations of the Cos class.
   * @see Cos
   */

  /**
   * @ingroup CosSpecializations
   */
  template <class NestedDerived>
  class Cos<FunctionBase<NestedDerived>>
    : public ScalarFunctionBase<Cos<FunctionBase<NestedDerived>>>
  {
    public:
      using Operand = FunctionBase<NestedDerived>;
      using Parent = ScalarFunctionBase<Cos<FunctionBase<NestedDerived>>>;

      Cos(const Operand& v)
        : m_v(v)
      {}

      Cos(const Cos& other)
        : Parent(other),
          m_v(other.m_v)
      {}

      Cos(Cos&& other)
        : Parent(std::move(other)),
          m_v(std::move(other.m_v))
      {}

      inline
      constexpr
      Cos& traceOf(Geometry::Attribute attrs)
      {
        m_v.traceOf(attrs);
        return *this;
      }

      inline
      auto getValue(const Geometry::Point& p) const
      {
        return std::cos(Scalar(m_v.getValue(p)));
      }

    private:
      Operand m_v;
  };

  template <class NestedDerived>
  Cos(const FunctionBase<NestedDerived>&) -> Cos<FunctionBase<NestedDerived>>;
}

#endif
