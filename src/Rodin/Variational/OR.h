/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_OR_H
#define RODIN_VARIATIONAL_OR_H

#include "ForwardDecls.h"
#include "BooleanFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup ORSpecializations OR Template Specializations
   * @brief Template specializations of the OR class.
   * @see OR
   */

  /**
   * @ingroup ORSpecializations
   * @brief Logical OR operator between two instances of BooleanFunctionBase
   */
  template <class LHSDerived, class RHSDerived>
  class OR<BooleanFunctionBase<LHSDerived>, BooleanFunctionBase<RHSDerived>> final
    : public BooleanFunctionBase<OR<BooleanFunctionBase<LHSDerived>, BooleanFunctionBase<RHSDerived>>>
  {
    public:
      using LHSType = BooleanFunctionBase<LHSDerived>;

      using RHSType = BooleanFunctionBase<RHSDerived>;

      using Parent = BooleanFunctionBase<OR<LHSType, RHSType>>;

      constexpr
      OR(const LHSType& lhs, const RHSType& rhs)
        : m_lhs(lhs.copy()), m_rhs(rhs.copy())
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
        return getLHS().getValue(p) || getRHS().getValue(p);
      }

      OR* copy() const noexcept final override
      {
        return new OR(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      std::unique_ptr<RHSType> m_rhs;
  };

  template <class LHSDerived, class RHSDerived>
  OR(const BooleanFunctionBase<LHSDerived>&, const BooleanFunctionBase<RHSDerived>&)
    -> OR<BooleanFunctionBase<LHSDerived>, BooleanFunctionBase<RHSDerived>>;

  template <class LHSDerived, class RHSDerived>
  constexpr
  auto
  operator||(const BooleanFunctionBase<LHSDerived>& lhs, const BooleanFunctionBase<RHSDerived>& rhs)
  {
    return OR(lhs, rhs);
  }

  template <class RHSDerived>
  constexpr
  auto
  operator||(Boolean lhs, const BooleanFunctionBase<RHSDerived>& rhs)
  {
    return OR(BooleanFunction(lhs), rhs);
  }

  template <class LHSDerived>
  constexpr
  auto
  operator||(const BooleanFunctionBase<LHSDerived>& lhs, Boolean rhs)
  {
    return OR(lhs, BooleanFunction(rhs));
  }
}

#endif



