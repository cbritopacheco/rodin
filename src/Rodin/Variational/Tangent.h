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
#include "RealFunction.h"

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
    : public RealFunctionBase<Tan<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using ScalarType = Real;

      using Parent = RealFunctionBase<Tan<FunctionBase<NestedDerived>>>;

      constexpr
      Tan(const OperandType& v)
        : m_operand(v.copy())
      {}

      constexpr
      Tan(const Tan& other)
        : Parent(other),
          m_operand(other.m_v->copy())
      {}

      constexpr
      Tan(Tan&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_v))
      {}

      inline
      constexpr
      Tan& traceOf(Geometry::Attribute attr)
      {
        m_operand->traceOf(attr);
        return *this;
      }

      inline
      constexpr
      Tan& traceOf(const FlatSet<Geometry::Attribute>& attrs)
      {
        m_operand->traceOf(attrs);
        return *this;
      }

      inline
      ScalarType getValue(const Geometry::Point& p) const
      {
        return Math::tan(Real(getOperand().getValue(p)));
      }

      inline
      const OperandType& getOperand() const
      {
        assert(m_operand);
        return *m_operand;
      }

      inline
      Tan* copy() const noexcept
      override
      {
        return new Tan(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
  };

  template <class NestedDerived>
  Tan(const FunctionBase<NestedDerived>&) -> Tan<FunctionBase<NestedDerived>>;
}

#endif

