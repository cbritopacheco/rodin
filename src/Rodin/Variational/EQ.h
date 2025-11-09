/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_EQ_H
#define RODIN_VARIATIONAL_EQ_H

/**
 * @file
 * @brief Equality comparison operator for functions.
 *
 * This file defines the equality operator @f$ == @f$ for comparing functions.
 * The result is a Boolean-valued function that evaluates to true where the
 * operands are equal and false otherwise.
 */

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
   * @brief Equality comparison operator between two functions.
   *
   * Represents the equality test @f$ f == g @f$ which evaluates to a Boolean
   * function:
   * @f[
   *   (f == g)(x) = \begin{cases}
   *     \text{true} & \text{if } f(x) = g(x) \\
   *     \text{false} & \text{otherwise}
   *   \end{cases}
   * @f]
   *
   * This is useful for:
   * - Defining indicator functions
   * - Implementing conditional logic in variational forms
   * - Creating characteristic functions of regions
   *
   * @tparam LHSDerived Type of left-hand side function
   * @tparam RHSDerived Type of right-hand side function
   *
   * @see GT, LT, GEQ, LEQ, AND, OR
   */
  template <class LHSDerived, class RHSDerived>
  class EQ<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>> final
    : public FunctionBase<EQ<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using LHSType = FunctionBase<LHSDerived>;

      using RHSType = FunctionBase<RHSDerived>;

      using Parent = FunctionBase<EQ<LHSType, RHSType>>;

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
}

#endif


