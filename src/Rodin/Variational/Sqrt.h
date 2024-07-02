/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_SQRT_H
#define RODIN_VARIATIONAL_SQRT_H

#include "Rodin/Math.h"
#include "ForwardDecls.h"
#include "Function.h"
#include "RealFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup SqrtSpecializations Sqrt Template Specializations
   * @brief Template specializations of the Sqrt class.
   * @see Sqrt
   */

  /**
   * @ingroup SqrtSpecializations
   */
  template <class NestedDerived>
  class Sqrt<FunctionBase<NestedDerived>> final
    : public RealFunctionBase<Sqrt<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using Parent = RealFunctionBase<Sqrt<FunctionBase<NestedDerived>>>;

      Sqrt(const OperandType& v)
        : m_v(v.copy())
      {}

      Sqrt(const Sqrt& other)
        : Parent(other),
          m_v(other.m_v->copy())
      {}

      Sqrt(Sqrt&& other)
        : Parent(std::move(other)),
          m_v(std::move(other.m_v))
      {}

      inline
      constexpr
      Sqrt& traceOf(Geometry::Attribute attrs)
      {
        m_v.traceOf(attrs);
        return *this;
      }

      inline
      auto getValue(const Geometry::Point& p) const
      {
        return std::sqrt(static_cast<Real>(getOperand().getValue(p)));
      }

      inline
      const OperandType& getOperand() const
      {
        assert(m_v);
        return *m_v;
      }

      inline Sqrt* copy() const noexcept override
      {
        return new Sqrt(*this);
      }

    private:
      std::unique_ptr<OperandType> m_v;
  };

  template <class NestedDerived>
  Sqrt(const FunctionBase<NestedDerived>&) -> Sqrt<FunctionBase<NestedDerived>>;

  /**
   * @brief Helper function to construct objects of type Sqrt.
   */
  template <class NestedDerived>
  auto sqrt(const FunctionBase<NestedDerived>& f)
  {
    return Sqrt(f);
  }
}

#endif


