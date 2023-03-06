/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_OR_H
#define RODIN_VARIATIONAL_OR_H

#include "ForwardDecls.h"
#include "GridFunction.h"
#include "BooleanFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup ORSpecializations OR Template Specializations
   * @brief Template specializations of the OR class.
   * @see OR
   */

  template <class LHSDerived, class RHSDerived>
  class OR<BooleanFunctionBase<LHSDerived>, BooleanFunctionBase<RHSDerived>> final
    : public BooleanFunctionBase<OR<BooleanFunctionBase<LHSDerived>, BooleanFunctionBase<RHSDerived>>>
  {
    public:
      using LHS = BooleanFunctionBase<LHSDerived>;
      using RHS = BooleanFunctionBase<RHSDerived>;
      using Parent = BooleanFunctionBase<OR<LHS, RHS>>;

      constexpr
      OR(const LHS& lhs, const RHS& rhs)
        : m_lhs(lhs), m_rhs(rhs)
      {}

      constexpr
      OR(const OR& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()),
          m_rhs(other.m_rhs->copy())
      {}

      constexpr
      OR(OR&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)),
          m_rhs(std::move(other.m_rhs))
      {}

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return Boolean(m_lhs.getValue(p)) || Boolean(m_rhs.getValue(p));
      }

      inline
      OR* copy() const noexcept
      override
      {
        return new OR(*this);
      }

    private:
      LHS m_lhs;
      RHS m_rhs;
  };

  template <class LHSDerived, class RHSDerived>
  OR(const BooleanFunctionBase<LHSDerived>&, const BooleanFunctionBase<RHSDerived>&)
    -> OR<BooleanFunctionBase<LHSDerived>, BooleanFunctionBase<RHSDerived>>;

  template <class LHSDerived, class RHSDerived>
  inline
  constexpr
  auto
  operator&&(const BooleanFunctionBase<LHSDerived>& lhs, const BooleanFunctionBase<RHSDerived>& rhs)
  {
    return OR(lhs, rhs);
  }

  template <class RHSDerived>
  inline
  constexpr
  auto
  operator&&(Boolean lhs, const BooleanFunctionBase<RHSDerived>& rhs)
  {
    return OR(BooleanFunction(lhs), rhs);
  }

  template <class LHSDerived>
  inline
  constexpr
  auto
  operator&&(const BooleanFunctionBase<LHSDerived>& lhs, Boolean rhs)
  {
    return OR(lhs, BooleanFunction(rhs));
  }
}

#endif



