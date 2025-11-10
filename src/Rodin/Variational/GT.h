/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file GT.h
 * @brief Greater-than comparison operator for functions.
 *
 * This file defines the GT class, which implements the greater-than comparison
 * operation between functions. This operator creates a boolean function that
 * indicates where one function exceeds another.
 *
 * ## Mathematical Foundation
 * For functions @f$ f, g : \Omega \to \mathbb{R} @f$, the greater-than operator:
 * @f[
 *   (f > g)(x) = \begin{cases} 1 & \text{if } f(x) > g(x) \\ 0 & \text{otherwise} \end{cases}
 * @f]
 *
 * ## Applications
 * - Indicator functions for regions: @f$ \mathbb{1}_{u > 0} @f$
 * - Contact problems: identify contact/separation zones
 * - Adaptive mesh refinement: flag elements where error exceeds tolerance
 * - Phase field methods: distinguish phases
 *
 * ## Usage Example
 * ```cpp
 * auto u = /* some function */;
 * auto positive_region = (u > 0.0);  // Boolean: true where u is positive
 * 
 * // Adaptive refinement indicator
 * auto refine_flag = (error_estimate > tolerance);
 * ```
 */
#ifndef RODIN_VARIATIONAL_GT_H
#define RODIN_VARIATIONAL_GT_H

#include "ForwardDecls.h"
#include "BooleanFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup GTSpecializations GT Template Specializations
   * @brief Template specializations of the GT class.
   * @see GT
   */

  /**
   * @ingroup GTSpecializations
   */
  template <class LHSDerived, class RHSDerived>
  class GT<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>> final
    : public BooleanFunctionBase<GT<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using LHSType = FunctionBase<LHSDerived>;

      using RHSType = FunctionBase<RHSDerived>;

      using Parent = BooleanFunctionBase<GT<LHSType, RHSType>>;

      GT(const LHSType& lhs, const RHSType& rhs)
        : m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {}

      GT(const GT& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()),
          m_rhs(other.m_rhs->copy())
      {}

      GT(GT&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)),
          m_rhs(std::move(other.m_rhs))
      {}

      constexpr
      Boolean getValue(const Geometry::Point& p) const
      {
        return getLHS().getValue(p) > getRHS().getValue(p);
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

      GT* copy() const noexcept final override
      {
        return new GT(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      std::unique_ptr<RHSType> m_rhs;
  };

  /**
   * @brief CTAD for GT.
   */
  template <class LHSDerived, class RHSDerived>
  GT(const FunctionBase<LHSDerived>&, const FunctionBase<RHSDerived>&)
    -> GT<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>;

  template <class LHSDerived, class RHSDerived>
  constexpr
  auto
  operator>(const FunctionBase<LHSDerived>& lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return GT(lhs, rhs);
  }

  template <class Number, class RHSDerived,
           typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  constexpr
  auto
  operator>(Number lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return GT(RealFunction(lhs), rhs);
  }

  template <class LHSDerived, class Number,
           typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  constexpr
  auto
  operator>(const FunctionBase<LHSDerived>& lhs, Number rhs)
  {
    return GT(lhs, RealFunction(rhs));
  }
}

#endif

