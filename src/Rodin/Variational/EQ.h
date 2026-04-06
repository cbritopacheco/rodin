/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file EQ.h
 * @brief Equality comparison operator for functions.
 *
 * This file defines the EQ class, which implements the equality comparison
 * operation between functions. This operator tests whether two function values
 * are equal at each point in the domain.
 *
 * ## Mathematical Foundation
 * For functions @f$ f, g : \Omega \to \mathbb{R} @f$, the equality operator:
 * @f[
 *   (f = g)(x) = \begin{cases} 1 & \text{if } f(x) = g(x) \\ 0 & \text{otherwise} \end{cases}
 * @f]
 *
 * ## Numerical Considerations
 * For floating-point comparisons, consider using approximate equality with
 * a tolerance due to numerical precision limitations.
 *
 * ## Applications
 * - Indicator functions for specific values
 * - Level set definitions
 * - Region identification based on function values
 *
 * ## Usage Example
 * ```cpp
 * auto is_zero = (u == 0.0);  // Boolean function, true where u equals zero
 * ```
 */
#ifndef RODIN_VARIATIONAL_EQ_H
#define RODIN_VARIATIONAL_EQ_H

#include "ForwardDecls.h"
#include "BooleanFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup EQSpecializations EQ Template Specializations
   * @brief Template specializations of the EQ class.
   * @see EQ
   */

  /**
   * @ingroup EQSpecializations
   * @brief Logical EQ operator between two instances of FunctionBase
   */
  template <class LHSDerived, class RHSDerived>
  class EQ<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>> final
    : public BooleanFunctionBase<EQ<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using LHSType = FunctionBase<LHSDerived>;

      using RHSType = FunctionBase<RHSDerived>;

      using Parent = BooleanFunctionBase<EQ<LHSType, RHSType>>;

      EQ(const LHSType& lhs, const RHSType& rhs)
        : m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {}

      EQ(const EQ& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      EQ(EQ&& other)
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
      auto getValue(const Geometry::Point& p) const
      {
        return this->getLHS().getValue(p) == this->getRHS().getValue(p);
      }

      EQ* copy() const noexcept final override
      {
        return new EQ(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      std::unique_ptr<RHSType> m_rhs;
  };

  template <class LHSDerived, class RHSDerived>
  EQ(const FunctionBase<LHSDerived>&, const FunctionBase<RHSDerived>&)
    -> EQ<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>;

  template <class LHSDerived, class RHSDerived>
  constexpr
  auto
  operator==(const FunctionBase<LHSDerived>& lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return EQ(lhs, rhs);
  }

  template <class RHSDerived>
  constexpr
  auto
  operator==(Boolean lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return EQ(BooleanFunction(lhs), rhs);
  }

  template <class LHSDerived>
  constexpr
  auto
  operator==(const FunctionBase<LHSDerived>& lhs, Boolean rhs)
  {
    return EQ(lhs, BooleanFunction(rhs));
  }

  template <class Number, class RHSDerived,
           typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  constexpr
  auto
  operator==(Number lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return EQ(RealFunction(lhs), rhs);
  }

  template <class LHSDerived, class Number,
           typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  constexpr
  auto
  operator==(const FunctionBase<LHSDerived>& lhs, Number rhs)
  {
    return EQ(lhs, RealFunction(rhs));
  }
}

#endif


