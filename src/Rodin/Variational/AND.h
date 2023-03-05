/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_AND_H
#define RODIN_VARIATIONAL_AND_H

#include "ForwardDecls.h"
#include "GridFunction.h"
#include "BooleanFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup ANDSpecializations AND Template Specializations
   * @brief Template specializations of the AND class.
   * @see AND
   */

  /**
   * @ingroup ANDSpecializations
   * @brief Logical AND operator between two instances of BooleanFunctionBase
   */
  template <class LHSDerived, class RHSDerived>
  class AND<BooleanFunctionBase<LHSDerived>, BooleanFunctionBase<RHSDerived>> final
    : public BooleanFunctionBase<AND<BooleanFunctionBase<LHSDerived>, BooleanFunctionBase<RHSDerived>>>
  {
    public:
      using LHS = BooleanFunctionBase<LHSDerived>;
      using RHS = BooleanFunctionBase<RHSDerived>;
      using Parent = BooleanFunctionBase<OR<LHS, RHS>>;

      AND(const LHS& lhs, const RHS& rhs)
        : m_lhs(lhs), m_rhs(rhs)
      {}

      AND(const AND& other)
        : Parent(other),
          m_lhs(other.m_lhs),
          m_rhs(other.m_rhs)
      {}

      AND(AND&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)),
          m_rhs(std::move(other.m_rhs))
      {}

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return Boolean(m_lhs.getValue(p)) && Boolean(m_rhs.getValue(p));
      }

    private:
      LHS m_lhs;
      RHS m_rhs;
  };

  template <class LHSDerived, class RHSDerived>
  AND(const BooleanFunctionBase<LHSDerived>&, const BooleanFunctionBase<RHSDerived>&)
    -> AND<BooleanFunctionBase<LHSDerived>, BooleanFunctionBase<RHSDerived>>;

  template <class LHSDerived, class RHSDerived>
  inline
  constexpr
  auto
  operator||(const BooleanFunctionBase<LHSDerived>& lhs, const BooleanFunctionBase<RHSDerived>& rhs)
  {
    return AND(lhs, rhs);
  }

  template <class RHSDerived>
  inline
  constexpr
  auto
  operator||(Boolean lhs, const BooleanFunctionBase<RHSDerived>& rhs)
  {
    return AND(BooleanFunction(lhs), rhs);
  }

  template <class LHSDerived>
  inline
  constexpr
  auto
  operator||(const BooleanFunctionBase<LHSDerived>& lhs, Boolean rhs)
  {
    return AND(lhs, BooleanFunction(rhs));
  }
}

#endif


