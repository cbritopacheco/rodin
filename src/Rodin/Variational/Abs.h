/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_ABS_H
#define RODIN_VARIATIONAL_ABS_H

#include <cmath>
#include "ForwardDecls.h"
#include "Function.h"
#include "ScalarFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup AbsSpecializations Abs Template Specializations
   * @brief Template specializations of the Abs class.
   * @see Abs
   */

  /**
   * @ingroup AbsSpecializations
   */
  template <class NestedDerived>
  class Abs<FunctionBase<NestedDerived>>
    : public ScalarFunctionBase<Abs<FunctionBase<NestedDerived>>>
  {
    public:
      using Operand = FunctionBase<NestedDerived>;
      using Parent = ScalarFunctionBase<Abs<Operand>>;

      using OperandRange = typename FormLanguage::Traits<Operand>::RangeType;
      static_assert(std::is_same_v<OperandRange, Scalar>);

      Abs(const Operand& v)
        : m_v(v.copy())
      {}

      Abs(const Abs& other)
        : Parent(other),
          m_v(other.m_v->copy())
      {}

      Abs(Abs&& other)
        : Parent(std::move(other)),
          m_v(std::move(other.m_v))
      {}

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return std::abs(getOperand().getValue(p));
      }

      const Operand& getOperand() const
      {
        assert(m_v);
        return *m_v;
      }

      inline Abs* copy() const noexcept override
      {
        return new Abs(*this);
      }

    private:
      std::unique_ptr<Operand> m_v;
  };

  template <class NestedDerived>
  Abs(const FunctionBase<NestedDerived>&) -> Abs<FunctionBase<NestedDerived>>;
}

namespace Rodin::Math
{
  template <class NestedDerived>
  inline
  constexpr auto
  abs(const Rodin::Variational::FunctionBase<NestedDerived>& op)
  {
    return Variational::Abs<Rodin::Variational::FunctionBase<NestedDerived>>(op);
  }
}

#endif


