/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_EQ_H
#define RODIN_VARIATIONAL_EQ_H

#include "ForwardDecls.h"
#include "GridFunction.h"
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
    : public FunctionBase<EQ<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using LHS = FunctionBase<LHSDerived>;
      using RHS = FunctionBase<RHSDerived>;
      using Parent = FunctionBase<EQ<LHS, RHS>>;

      EQ(const LHS& lhs, const RHS& rhs)
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
        return getLHS().getValue(p) == getRHS().getValue(p);
      }

      inline EQ* copy() const noexcept final override
      {
        return new EQ(*this);
      }

    private:
      std::unique_ptr<LHS> m_lhs;
      std::unique_ptr<RHS> m_rhs;
  };

  template <class LHSDerived, class RHSDerived>
  EQ(const FunctionBase<LHSDerived>&, const FunctionBase<RHSDerived>&)
    -> EQ<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>;

  template <class LHSDerived, class RHSDerived>
  inline
  constexpr
  auto
  operator==(const FunctionBase<LHSDerived>& lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return EQ(lhs, rhs);
  }

  template <class RHSDerived>
  inline
  constexpr
  auto
  operator==(Boolean lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return EQ(BooleanFunction(lhs), rhs);
  }

  template <class LHSDerived>
  inline
  constexpr
  auto
  operator==(const FunctionBase<LHSDerived>& lhs, Boolean rhs)
  {
    return EQ(lhs, BooleanFunction(rhs));
  }
}

#endif


