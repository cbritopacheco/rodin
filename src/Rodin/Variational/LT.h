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

  /**
   * @ingroup LTSpecializations
   */
  template <class LHSDerived, class RHSDerived>
  class LT<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>> final
    : public BooleanFunctionBase<LT<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using LHSType = FunctionBase<LHSDerived>;

      using RHSType = FunctionBase<RHSDerived>;

      using Parent = BooleanFunctionBase<LT<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>;

      LT(const LHSType& lhs, const RHSType& rhs)
        : m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {}

      LT(const LT& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()),
          m_rhs(other.m_rhs->copy())
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
        return getLHS().getValue(p) < getRHS().getValue(p);
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

      inline LT* copy() const noexcept override
      {
        return new LT(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      std::unique_ptr<RHSType> m_rhs;
  };

  /**
   * @brief CTAD for LT.
   */
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
    return LT(RealFunction(lhs), rhs);
  }

  template <class LHSDerived, class Number,
           typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  inline
  constexpr
  auto
  operator<(const FunctionBase<LHSDerived>& lhs, Number rhs)
  {
    return LT(lhs, RealFunction(rhs));
  }
}

#endif

