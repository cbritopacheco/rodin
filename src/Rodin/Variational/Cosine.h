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
#include "RealFunction.h"

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
  class Cos<FunctionBase<NestedDerived>> final
    : public RealFunctionBase<Cos<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using Parent = RealFunctionBase<Cos<FunctionBase<NestedDerived>>>;

      Cos(const OperandType& v)
        : m_operand(v.copy())
      {}

      Cos(const Cos& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      Cos(Cos&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}

      inline
      constexpr
      Cos& traceOf(Geometry::Attribute attr)
      {
        m_operand->traceOf(attr);
        return *this;
      }

      inline
      constexpr
      Cos& traceOf(const FlatSet<Geometry::Attribute>& attrs)
      {
        m_operand->traceOf(attrs);
        return *this;
      }

      inline
      Real getValue(const Geometry::Point& p) const
      {
        return Math::cos(getOperand().getValue(p));
      }

      inline
      const OperandType& getOperand() const
      {
        assert(m_operand);
        return *m_operand;
      }

      inline Cos* copy() const noexcept override
      {
        return new Cos(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
  };

  template <class NestedDerived>
  Cos(const FunctionBase<NestedDerived>&) -> Cos<FunctionBase<NestedDerived>>;

  /**
   * @brief Helper function to construct objects of type Cos.
   */
  template <class NestedDerived>
  auto cos(const FunctionBase<NestedDerived>& f)
  {
    return Cos(f);
  }
}

#endif
