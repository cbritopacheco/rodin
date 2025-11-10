/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file GEQ.h
 * @brief Greater-than-or-equal comparison operator for functions.
 *
 * This file defines the GEQ class, which implements the greater-than-or-equal
 * comparison operation between functions. This operator creates a boolean
 * function indicating where one function is at least as large as another.
 *
 * ## Mathematical Foundation
 * For functions @f$ f, g : \Omega \to \mathbb{R} @f$, the GEQ operator:
 * @f[
 *   (f \geq g)(x) = \begin{cases} 1 & \text{if } f(x) \geq g(x) \\ 0 & \text{otherwise} \end{cases}
 * @f]
 *
 * ## Applications
 * - Non-negative region indicators: @f$ \mathbb{1}_{u \geq 0} @f$
 * - Constraint satisfaction checking
 * - Admissible set definitions
 * - Complementarity conditions
 *
 * ## Usage Example
 * ```cpp
 * auto u = /* some function */;
 * auto non_negative = (u >= 0.0);  // Boolean: true where u is non-negative
 * ```
 */
#ifndef RODIN_VARIATIONAL_GEQ_H
#define RODIN_VARIATIONAL_GEQ_H

#include "ForwardDecls.h"
#include "BooleanFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup GEQSpecializations GEQ Template Specializations
   * @brief Template specializations of the GEQ class.
   * @see GEQ
   */

  /**
   * @ingroup GEQSpecializations
   */
  template <class LHSDerived, class RHSDerived>
  class GEQ<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>> final
    : public BooleanFunctionBase<GEQ<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using LHSType = FunctionBase<LHSDerived>;

      using RHSType = FunctionBase<RHSDerived>;

      using Parent = BooleanFunctionBase<GEQ<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>;

      GEQ(const LHSType& lhs, const RHSType& rhs)
        : m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {}

      GEQ(const GEQ& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      GEQ(GEQ&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)),
          m_rhs(std::move(other.m_rhs))
      {}

      constexpr
      Boolean getValue(const Geometry::Point& p) const
      {
        return getLHS().getValue(p) >= getRHS().getValue(p);
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

      GEQ* copy() const noexcept final override
      {
        return new GEQ(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      std::unique_ptr<RHSType> m_rhs;
  };

  /**
   * @brief CTAD for GEQ.
   */
  template <class LHSDerived, class RHSDerived>
  GEQ(const FunctionBase<LHSDerived>&, const FunctionBase<RHSDerived>&)
    -> GEQ<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>;

  template <class LHSDerived, class RHSDerived>
  constexpr
  auto
  operator>=(const FunctionBase<LHSDerived>& lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return GEQ(lhs, rhs);
  }

  template <class Number, class RHSDerived,
           typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  constexpr
  auto
  operator>=(Number lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return GEQ(RealFunction(lhs), rhs);
  }

  template <class LHSDerived, class Number,
           typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  constexpr
  auto
  operator>=(const FunctionBase<LHSDerived>& lhs, Number rhs)
  {
    return GEQ(lhs, RealFunction(rhs));
  }
}

#endif

