/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file NEQ.h
 * @brief Inequality comparison operator for functions.
 *
 * This file defines the NEQ class, which implements the not-equal comparison
 * operation between functions. This operator tests whether two function values
 * differ at each point in the domain.
 *
 * ## Mathematical Foundation
 * For functions @f$ f, g : \Omega \to \mathbb{R} @f$, the inequality operator:
 * @f[
 *   (f \neq g)(x) = \begin{cases} 1 & \text{if } f(x) \neq g(x) \\ 0 & \text{otherwise} \end{cases}
 * @f]
 *
 * ## Numerical Considerations
 * For floating-point comparisons, consider using approximate equality tests
 * rather than strict inequality, due to numerical precision limitations.
 *
 * ## Applications
 * - Complement of indicator functions
 * - Non-zero region detection
 * - Logical complement of equality conditions
 *
 * ## Usage Example
 * ```cpp
 * auto is_nonzero = (u != 0.0);  // Boolean function, true where u is non-zero
 * ```
 */
#ifndef RODIN_VARIATIONAL_NEQ_H
#define RODIN_VARIATIONAL_NEQ_H

#include "ForwardDecls.h"
#include "BooleanFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup NEQSpecializations NEQ Template Specializations
   * @brief Template specializations of the NEQ class.
   * @see NEQ
   */

  /**
   * @ingroup NEQSpecializations
   * @brief Logical NEQ operator between two instances of FunctionBase
   */
  template <class LHSDerived, class RHSDerived>
  class NEQ<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>> final
    : public BooleanFunctionBase<NEQ<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using LHSType = FunctionBase<LHSDerived>;

      using RHSType = FunctionBase<RHSDerived>;

      using Parent = BooleanFunctionBase<NEQ<LHSType, RHSType>>;

      NEQ(const LHSType& lhs, const RHSType& rhs)
        : m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {}

      NEQ(const NEQ& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      NEQ(NEQ&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)),
          m_rhs(std::move(other.m_rhs))
      {}

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

      constexpr
      Boolean getValue(const Geometry::Point& p) const
      {
        return getLHS().getValue(p) != getRHS().getValue(p);
      }

      NEQ* copy() const noexcept final override
      {
        return new NEQ(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      std::unique_ptr<RHSType> m_rhs;
  };

  template <class LHSDerived, class RHSDerived>
  NEQ(const FunctionBase<LHSDerived>&, const FunctionBase<RHSDerived>&)
    -> NEQ<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>;

  template <class LHSDerived, class RHSDerived>
  constexpr
  auto
  operator!=(const FunctionBase<LHSDerived>& lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return NEQ(lhs, rhs);
  }

  template <class RHSDerived>
  constexpr
  auto
  operator!=(Boolean lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return NEQ(BooleanFunction(lhs), rhs);
  }

  template <class LHSDerived>
  constexpr
  auto
  operator!=(const FunctionBase<LHSDerived>& lhs, Boolean rhs)
  {
    return NEQ(lhs, BooleanFunction(rhs));
  }

  template <class Number, class RHSDerived,
           typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  constexpr
  auto
  operator!=(Number lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return NEQ(RealFunction(lhs), rhs);
  }

  template <class LHSDerived, class Number,
           typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  constexpr
  auto
  operator!=(const FunctionBase<LHSDerived>& lhs, Number rhs)
  {
    return NEQ(lhs, RealFunction(rhs));
  }
}

#endif
