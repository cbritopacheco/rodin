/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_LEQ_H
#define RODIN_VARIATIONAL_LEQ_H

#include "ForwardDecls.h"
#include "GridFunction.h"
#include "BooleanFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup LEQSpecializations LEQ Template Specializations
   * @brief Template specializations of the LEQ class.
   * @see LEQ
   */

  /**
   * @ingroup LEQSpecializations
   */
  template <class LHSDerived, class RHSDerived>
  class LEQ<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>> final
    : public BooleanFunctionBase<LEQ<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using Parent = BooleanFunctionBase<LEQ<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>;
      using LHS = FunctionBase<LHSDerived>;
      using RHS = FunctionBase<RHSDerived>;

      LEQ(const LHS& lhs, const RHS& rhs)
        : m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {}

      LEQ(const LEQ& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()),
          m_rhs(other.m_rhs->copy())
      {}

      LEQ(LEQ&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)),
          m_rhs(std::move(other.m_rhs))
      {}

      inline
      constexpr
      Boolean getValue(const Geometry::Point& p) const
      {
        return getLHS().getValue(p) <= getRHS().getValue(p);
      }

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

      inline LEQ* copy() const noexcept override
      {
        return new LEQ(*this);
      }

    private:
      std::unique_ptr<LHS> m_lhs;
      std::unique_ptr<RHS> m_rhs;
  };

  template <class LHSDerived, class RHSDerived>
  LEQ(const FunctionBase<LHSDerived>&, const FunctionBase<RHSDerived>&)
    -> LEQ<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>;

  template <class LHSDerived, class RHSDerived>
  inline
  constexpr
  auto
  operator<=(const FunctionBase<LHSDerived>& lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return LEQ(lhs, rhs);
  }

  template <class Number, class RHSDerived,
           typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  inline
  constexpr
  auto
  operator<=(Number lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return LEQ(ScalarFunction(lhs), rhs);
  }

  template <class LHSDerived, class Number,
           typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  inline
  constexpr
  auto
  operator<=(const FunctionBase<LHSDerived>& lhs, Number rhs)
  {
    return LEQ(lhs, ScalarFunction(rhs));
  }
}

#endif

