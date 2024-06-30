/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GT_H
#define RODIN_VARIATIONAL_GT_H

#include "ForwardDecls.h"
#include "GridFunction.h"
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

      inline
      constexpr
      Boolean getValue(const Geometry::Point& p) const
      {
        return getLHS().getValue(p) > getRHS().getValue(p);
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

      inline GT* copy() const noexcept final override
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
  inline
  constexpr
  auto
  operator>(const FunctionBase<LHSDerived>& lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return GT(lhs, rhs);
  }

  template <class Number, class RHSDerived,
           typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  inline
  constexpr
  auto
  operator>(Number lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return GT(ScalarFunction(lhs), rhs);
  }

  template <class LHSDerived, class Number,
           typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  inline
  constexpr
  auto
  operator>(const FunctionBase<LHSDerived>& lhs, Number rhs)
  {
    return GT(lhs, ScalarFunction(rhs));
  }
}

#endif

