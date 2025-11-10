/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file LT.h
 * @brief Less-than comparison operator for functions.
 *
 * This file defines the LT class, which implements the less-than comparison
 * operation between functions. This operator creates a boolean function
 * indicating where one function is less than another.
 *
 * ## Mathematical Foundation
 * For functions @f$ f, g : \Omega \to \mathbb{R} @f$, the less-than operator:
 * @f[
 *   (f < g)(x) = \begin{cases} 1 & \text{if } f(x) < g(x) \\ 0 & \text{otherwise} \end{cases}
 * @f]
 *
 * ## Applications
 * - Negative region indicators: @f$ \mathbb{1}_{u < 0} @f$
 * - Subdomain identification
 * - Constraint violation detection
 * - Thresholding operations
 *
 * ## Usage Example
 * ```cpp
 * auto u = /* some function */;
 * auto negative_region = (u < 0.0);  // Boolean: true where u is negative
 * ```
 */
#ifndef RODIN_VARIATIONAL_LT_H
#define RODIN_VARIATIONAL_LT_H

#include "ForwardDecls.h"
#include "BooleanFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup LTSpecializations LT Template Specializations
   * @brief Template specializations of the LT class.
   * @see LT
   */

  /**
   * @ingroup LTSpecializations
   */
  template <class LHSDerived, class RHSDerived>
  class LT<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>> final
    : public BooleanFunctionBase<LT<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using LHSType = FunctionBase<LHSDerived>;

      using RHSType = FunctionBase<RHSDerived>;

      using Parent = BooleanFunctionBase<LT<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>;

      LT(const LHSType& lhs, const RHSType& rhs)
        : m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {}

      LT(const LT& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()),
          m_rhs(other.m_rhs->copy())
      {}

      LT(LT&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)),
          m_rhs(std::move(other.m_rhs))
      {}

      constexpr
      Boolean getValue(const Geometry::Point& p) const
      {
        return getLHS().getValue(p) < getRHS().getValue(p);
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

      LT* copy() const noexcept override
      {
        return new LT(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      std::unique_ptr<RHSType> m_rhs;
  };

  /**
   * @brief CTAD for LT.
   */
  template <class LHSDerived, class RHSDerived>
  LT(const FunctionBase<LHSDerived>&, const FunctionBase<RHSDerived>&)
    -> LT<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>;

  template <class LHSDerived, class RHSDerived>
  constexpr
  auto
  operator<(const FunctionBase<LHSDerived>& lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return LT(lhs, rhs);
  }

  template <class Number, class RHSDerived,
           typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  constexpr
  auto
  operator<(Number lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return LT(RealFunction(lhs), rhs);
  }

  template <class LHSDerived, class Number,
           typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  constexpr
  auto
  operator<(const FunctionBase<LHSDerived>& lhs, Number rhs)
  {
    return LT(lhs, RealFunction(rhs));
  }
}

#endif

