/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file LEQ.h
 * @brief Less-than-or-equal comparison operator for functions.
 *
 * This file defines the LEQ class, which implements the less-than-or-equal
 * comparison operation between functions. This operator creates a boolean
 * function indicating where one function does not exceed another.
 *
 * ## Mathematical Foundation
 * For functions @f$ f, g : \Omega \to \mathbb{R} @f$, the LEQ operator:
 * @f[
 *   (f \leq g)(x) = \begin{cases} 1 & \text{if } f(x) \leq g(x) \\ 0 & \text{otherwise} \end{cases}
 * @f]
 *
 * ## Applications
 * - Bounded region indicators: @f$ \mathbb{1}_{u \leq M} @f$
 * - Feasibility checking for optimization problems
 * - Level set definitions: @f$ \{x : \phi(x) \leq 0\} @f$
 * - Truncation and clipping operations
 *
 * ## Usage Example
 * ```cpp
 * auto bounded_region = (u <= 1.0);  // Boolean: true where u ≤ 1
 * ```
 */
#ifndef RODIN_VARIATIONAL_LEQ_H
#define RODIN_VARIATIONAL_LEQ_H

#include "ForwardDecls.h"
#include "BooleanFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup LEQSpecializations LEQ Template Specializations
   * @brief Template specializations of the LEQ class.
   * @see LEQ
   */

  /**
   * @ingroup LEQSpecializations
   */
  template <class LHSDerived, class RHSDerived>
  class LEQ<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>> final
    : public BooleanFunctionBase<LEQ<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using LHSType = FunctionBase<LHSDerived>;

      using RHSType = FunctionBase<RHSDerived>;

      using Parent = BooleanFunctionBase<LEQ<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>;

      LEQ(const LHSType& lhs, const RHSType& rhs)
        : m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {}

      LEQ(const LEQ& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()),
          m_rhs(other.m_rhs->copy())
      {}

      LEQ(LEQ&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)),
          m_rhs(std::move(other.m_rhs))
      {}

      constexpr
      Boolean getValue(const Geometry::Point& p) const
      {
        return getLHS().getValue(p) <= getRHS().getValue(p);
      }

      const auto& getLHS() const
      {
        assert(m_lhs);
        return *m_lhs;
      }

      const auto& getRHS() const
      {
        assert(m_rhs);
        return *m_rhs;
      }

      LEQ* copy() const noexcept override
      {
        return new LEQ(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      std::unique_ptr<RHSType> m_rhs;
  };

  template <class LHSDerived, class RHSDerived>
  LEQ(const FunctionBase<LHSDerived>&, const FunctionBase<RHSDerived>&)
    -> LEQ<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>;

  template <class LHSDerived, class RHSDerived>
  constexpr
  auto
  operator<=(const FunctionBase<LHSDerived>& lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return LEQ(lhs, rhs);
  }

  template <class Number, class RHSDerived,
           typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  constexpr
  auto
  operator<=(Number lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return LEQ(RealFunction(lhs), rhs);
  }

  template <class LHSDerived, class Number,
           typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  constexpr
  auto
  operator<=(const FunctionBase<LHSDerived>& lhs, Number rhs)
  {
    return LEQ(lhs, RealFunction(rhs));
  }
}

#endif

