/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_TANGENT_H
#define RODIN_VARIATIONAL_TANGENT_H

#include "Rodin/Math.h"
#include "ForwardDecls.h"
#include "Function.h"
#include "ScalarFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup TanSpecializations Tan Template Specializations
   * @brief Template specializations of the Tan class.
   * @see Tan
   */

  /**
   * @ingroup TanSpecializations
   */
  template <class NestedDerived>
  class Tan<FunctionBase<NestedDerived>>
    : public ScalarFunctionBase<Tan<FunctionBase<NestedDerived>>>
  {
    public:
      using Operand = FunctionBase<NestedDerived>;
      using Parent = ScalarFunctionBase<Tan<FunctionBase<NestedDerived>>>;

      constexpr
      Tan(const Operand& v)
        : m_v(v)
      {}

      constexpr
      Tan(const Tan& other)
        : Parent(other),
          m_v(other.m_v)
      {}

      constexpr
      Tan(Tan&& other)
        : Parent(std::move(other)),
          m_v(std::move(other.m_v))
      {}

      inline
      constexpr
      Tan& traceOf(Geometry::Attribute attrs)
      {
        m_v.traceOf(attrs);
        return *this;
      }

      inline
      Scalar getValue(const Geometry::Point& p) const
      {
        return std::tan(Scalar(m_v.getValue(p)));
      }

      inline
      Tan* copy() const noexcept
      override
      {
        return new Tan(*this);
      }

    private:
      Operand m_v;
  };

  template <class NestedDerived>
  Tan(const FunctionBase<NestedDerived>&) -> Tan<FunctionBase<NestedDerived>>;
}

#endif

