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
      using LHSType = BooleanFunctionBase<LHSDerived>;

      using RHSType = BooleanFunctionBase<RHSDerived>;

      using Parent = BooleanFunctionBase<AND<LHSType, RHSType>>;

      AND(const LHSType& lhs, const RHSType& rhs)
        : m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {}

      AND(const AND& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      AND(AND&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)),
          m_rhs(std::move(other.m_rhs))
      {}

      inline
      const auto& getLHS() const
      {
        assert(m_lhs);
        return *m_lhs;
      }

      inline
      const auto& getRHS() const
      {
        assert(m_rhs);
        return *m_rhs;
      }

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return getLHS().getValue(p) && getRHS().getValue(p);
      }

      inline AND* copy() const noexcept final override
      {
        return new AND(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      std::unique_ptr<RHSType> m_rhs;
  };

  template <class LHSDerived, class RHSDerived>
  AND(const BooleanFunctionBase<LHSDerived>&, const BooleanFunctionBase<RHSDerived>&)
    -> AND<BooleanFunctionBase<LHSDerived>, BooleanFunctionBase<RHSDerived>>;

  template <class LHSDerived, class RHSDerived>
  inline
  constexpr
  auto
  operator&&(const BooleanFunctionBase<LHSDerived>& lhs, const BooleanFunctionBase<RHSDerived>& rhs)
  {
    return AND(lhs, rhs);
  }

  template <class RHSDerived>
  inline
  constexpr
  auto
  operator&&(Boolean lhs, const BooleanFunctionBase<RHSDerived>& rhs)
  {
    return AND(BooleanFunction(lhs), rhs);
  }

  template <class LHSDerived>
  inline
  constexpr
  auto
  operator&&(const BooleanFunctionBase<LHSDerived>& lhs, Boolean rhs)
  {
    return AND(lhs, BooleanFunction(rhs));
  }
}

#endif


