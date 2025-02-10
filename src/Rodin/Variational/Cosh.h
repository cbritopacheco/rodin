/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_COSH_H
#define RODIN_VARIATIONAL_COSH_H

#include "Rodin/Math.h"
#include "ForwardDecls.h"
#include "Function.h"
#include "RealFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup CoshSpecializations Cosh Template Specializations
   * @brief Template specializations of the Cosh class.
   * @see Cosh
   */

  /**
   * @ingroup CoshSpecializations
   */
  template <class NestedDerived>
  class Cosh<FunctionBase<NestedDerived>> final
    : public RealFunctionBase<Cosh<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using Parent = RealFunctionBase<Cosh<FunctionBase<NestedDerived>>>;

      Cosh(const OperandType& v)
        : m_operand(v.copy())
      {}

      Cosh(const Cosh& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      Cosh(Cosh&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}

      constexpr
      Cosh& traceOf(Geometry::Attribute attr)
      {
        m_operand->traceOf(attr);
        return *this;
      }

      constexpr
      Cosh& traceOf(const FlatSet<Geometry::Attribute>& attrs)
      {
        m_operand->traceOf(attrs);
        return *this;
      }

      Real getValue(const Geometry::Point& p) const
      {
        return Math::cosh(getOperand().getValue(p));
      }

      const OperandType& getOperand() const
      {
        assert(m_operand);
        return *m_operand;
      }

      Cosh* copy() const noexcept override
      {
        return new Cosh(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
  };

  template <class NestedDerived>
  Cosh(const FunctionBase<NestedDerived>&) -> Cosh<FunctionBase<NestedDerived>>;

  /**
   * @brief Helper function to construct objects of type Cosh.
   */
  template <class NestedDerived>
  auto cosh(const FunctionBase<NestedDerived>& f)
  {
    return Cosh(f);
  }
}

#endif
