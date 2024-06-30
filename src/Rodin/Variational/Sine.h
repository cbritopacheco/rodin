/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_SIN_H
#define RODIN_VARIATIONAL_SIN_H

#include "Rodin/Math.h"
#include "ForwardDecls.h"
#include "Function.h"
#include "ScalarFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup SinSpecializations Sin Template Specializations
   * @brief Template specializations of the Sin class.
   * @see Sin
   */

  /**
   * @ingroup SinSpecializations
   */
  template <class NestedDerived>
  class Sin<FunctionBase<NestedDerived>> final
    : public ScalarFunctionBase<Sin<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using Parent = ScalarFunctionBase<Sin<FunctionBase<NestedDerived>>>;

      Sin(const OperandType& v)
        : m_v(v.copy())
      {}

      Sin(const Sin& other)
        : Parent(other),
          m_v(other.m_v->copy())
      {}

      Sin(Sin&& other)
        : Parent(std::move(other)),
          m_v(std::move(other.m_v))
      {}

      inline
      constexpr
      Sin& traceOf(Geometry::Attribute attrs)
      {
        m_v.traceOf(attrs);
        return *this;
      }

      inline
      auto getValue(const Geometry::Point& p) const
      {
        return std::sin(static_cast<Scalar>(getOperand().getValue(p)));
      }

      inline
      const OperandType& getOperand() const
      {
        assert(m_v);
        return *m_v;
      }

      inline Sin* copy() const noexcept override
      {
        return new Sin(*this);
      }

    private:
      std::unique_ptr<OperandType> m_v;
  };

  template <class NestedDerived>
  Sin(const FunctionBase<NestedDerived>&) -> Sin<FunctionBase<NestedDerived>>;

  /**
   * @brief Helper function to construct objects of type Sin.
   */
  template <class NestedDerived>
  auto sin(const FunctionBase<NestedDerived>& f)
  {
    return Sin(f);
  }
}

#endif

