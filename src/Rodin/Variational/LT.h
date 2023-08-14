/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_LT_H
#define RODIN_VARIATIONAL_LT_H

#include "ForwardDecls.h"
#include "GridFunction.h"
#include "BooleanFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup LTSpecializations LT Template Specializations
   * @brief Template specializations of the LT class.
   * @see LT
   */

  template <class LHSDerived, class RHSDerived>
  class LT<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>> final
    : public BooleanFunctionBase<LT<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using Parent = BooleanFunctionBase<LT<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>;
      using LHS = FunctionBase<LHSDerived>;
      using RHS = FunctionBase<RHSDerived>;

      LT(const LHS& lhs, const RHS& rhs)
        : m_lhs(lhs), m_rhs(rhs)
      {}

      LT(const LT& other)
        : Parent(other),
          m_lhs(other.m_lhs),
          m_rhs(other.m_rhs)
      {}

      LT(LT&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)),
          m_rhs(std::move(other.m_rhs))
      {}

      inline
      constexpr
      Boolean getValue(const Geometry::Point& p) const
      {
        return Scalar(m_lhs.getValue(p)) < Scalar(m_rhs.getValue(p));
      }

      inline
      LT* copy() const noexcept
      override
      {
        return new LT(*this);
      }

    private:
      LHS m_lhs;
      RHS m_rhs;
  };

  template <class LHSDerived, class RHSDerived>
  LT(const FunctionBase<LHSDerived>&, const FunctionBase<RHSDerived>&)
    -> LT<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>;

  template <class LHSDerived, class RHSDerived>
  inline
  constexpr
  auto
  operator<(const FunctionBase<LHSDerived>& lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return LT(lhs, rhs);
  }

  template <class Number, class RHSDerived,
           typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  inline
  constexpr
  auto
  operator<(Number lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return LT(ScalarFunction(lhs), rhs);
  }

  template <class LHSDerived, class Number,
           typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  inline
  constexpr
  auto
  operator<(const FunctionBase<LHSDerived>& lhs, Number rhs)
  {
    return LT(lhs, ScalarFunction(rhs));
  }
}

#endif

