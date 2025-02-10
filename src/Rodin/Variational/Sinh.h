/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_SINH_H
#define RODIN_VARIATIONAL_SINH_H

#include "Rodin/Math.h"
#include "ForwardDecls.h"
#include "Function.h"
#include "RealFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup SinhSpecializations Sinh Template Specializations
   * @brief Template specializations of the Sinh class.
   * @see Sinh
   */

  /**
   * @ingroup SinhSpecializations
   */
  template <class NestedDerived>
  class Sinh<FunctionBase<NestedDerived>> final
    : public RealFunctionBase<Sinh<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using Parent = RealFunctionBase<Sinh<FunctionBase<NestedDerived>>>;

      Sinh(const OperandType& v)
        : m_operand(v.copy())
      {}

      Sinh(const Sinh& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      Sinh(Sinh&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}

      constexpr
      Sinh& traceOf(Geometry::Attribute attr)
      {
        m_operand->traceOf(attr);
        return *this;
      }

      constexpr
      Sinh& traceOf(const FlatSet<Geometry::Attribute>& attrs)
      {
        m_operand->traceOf(attrs);
        return *this;
      }

      Real getValue(const Geometry::Point& p) const
      {
        return Math::sinh(getOperand().getValue(p));
      }

      const OperandType& getOperand() const
      {
        assert(m_operand);
        return *m_operand;
      }

      Sinh* copy() const noexcept override
      {
        return new Sinh(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
  };

  template <class NestedDerived>
  Sinh(const FunctionBase<NestedDerived>&) -> Sinh<FunctionBase<NestedDerived>>;

  /**
   * @brief Helper function to construct objects of type Sinh.
   */
  template <class NestedDerived>
  auto sinh(const FunctionBase<NestedDerived>& f)
  {
    return Sinh(f);
  }
}

#endif
