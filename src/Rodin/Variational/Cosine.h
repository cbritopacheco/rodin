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
  class Cos<FunctionBase<NestedDerived>> final
    : public ScalarFunctionBase<Cos<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using Parent = ScalarFunctionBase<Cos<FunctionBase<NestedDerived>>>;

      Cos(const OperandType& v)
        : m_v(v.copy())
      {}

      Cos(const Cos& other)
        : Parent(other),
          m_v(other.m_v->copy())
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
        return std::cos(static_cast<Scalar>(getOperand().getValue(p)));
      }

      inline
      const OperandType& getOperand() const
      {
        assert(m_v);
        return *m_v;
      }

      inline Cos* copy() const noexcept override
      {
        return new Cos(*this);
      }

    private:
      std::unique_ptr<OperandType> m_v;
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
