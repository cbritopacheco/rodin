/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_SQRT_H
#define RODIN_VARIATIONAL_SQRT_H

/**
 * @file
 * @brief Square root operations for functions.
 */

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
   * @brief Square root operation for functions.
   *
   * This class represents the square root operation applied to a scalar function:
   * @f[
   *    \text{Sqrt}(f)(x) = \sqrt{f(x)}
   * @f]
   *
   * The input function must be real-valued and non-negative for real results.
   * The result is always a real-valued scalar function.
   *
   * @tparam NestedDerived Type of the operand function
   *
   * @note For vector or matrix functions, consider using the Frobenius norm first.
   *
   * @see RealFunctionBase, Frobenius
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

      template <class ... Args>
      constexpr
      Sqrt& traceOf(const Args& ... args)
      {
        m_v->traceOf(args...);
        return *this;
      }

      Real getValue(const Geometry::Point& p) const
      {
        return Math::sqrt(this->getOperand().getValue(p));
      }

      const OperandType& getOperand() const
      {
        assert(m_v);
        return *m_v;
      }

      Sqrt* copy() const noexcept override
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


