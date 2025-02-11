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
#include "RealFunction.h"

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
    : public RealFunctionBase<Sin<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using Parent = RealFunctionBase<Sin<FunctionBase<NestedDerived>>>;

      Sin(const OperandType& v)
        : m_operand(v.copy())
      {}

      Sin(const Sin& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      Sin(Sin&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}

      constexpr
      Sin& traceOf(Geometry::Attribute attr)
      {
        m_operand->traceOf(attr);
        return *this;
      }

      constexpr
      Sin& traceOf(const FlatSet<Geometry::Attribute>& attrs)
      {
        m_operand->traceOf(attrs);
        return *this;
      }

      Real getValue(const Geometry::Point& p) const
      {
        return Math::sin(getOperand().getValue(p));
      }

      const OperandType& getOperand() const
      {
        assert(m_operand);
        return *m_operand;
      }

      Sin* copy() const noexcept override
      {
        return new Sin(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
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

