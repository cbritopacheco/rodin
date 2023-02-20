/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DIVISION_H
#define RODIN_VARIATIONAL_DIVISION_H

#include "ForwardDecls.h"
#include "Function.h"

namespace Rodin::Variational
{
  /**
   * @defgroup DivisionSpecializations Division Template Specializations
   * @brief Template specializations of the Division class.
   * @see Division
   */

  /**
   * @ingroup DivSpecializations
   * @brief Division of a FunctionBase by a FunctionBase.
   */
  template <class LHSDerived, class RHSDerived>
  class Division<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>
    : public FunctionBase<Division<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using Parent = FunctionBase<Division<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>;
      using LHS = FunctionBase<LHSDerived>;
      using RHS = FunctionBase<RHSDerived>;

      constexpr
      Division(const FunctionBase<LHSDerived>& lhs, const FunctionBase<RHSDerived>& rhs)
        : m_lhs(lhs), m_rhs(rhs)
      {}

      constexpr
      Division(const Division& other)
        : Parent(other),
          m_lhs(other.m_lhs), m_rhs(other.m_rhs)
      {}

      constexpr
      Division(Division&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return m_lhs.getRangeShape();
      }

      inline
      constexpr
      Division& traceOf(Geometry::Attribute attr)
      {
        m_lhs.traceOf(attr);
        m_rhs.traceOf(attr);
        return *this;
      }

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        using RHSRange = FormLanguage::RangeOf<typename FormLanguage::Traits<RHS>::ResultType>;
        static_assert(FormLanguage::IsScalarRange<RHSRange>::Value);
        return m_lhs.getValue(p) / m_rhs.getValue(p);
      }

    private:
      FunctionBase<LHSDerived> m_lhs;
      FunctionBase<RHSDerived> m_rhs;
  };
  template <class LHSDerived, class RHSDerived>
  Division(const FunctionBase<LHSDerived>&, const FunctionBase<RHSDerived>&)
    -> Division<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>;

  template <class LHSDerived, class RHSDerived>
  inline
  auto
  operator/(const FunctionBase<LHSDerived>& lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return Division(lhs, rhs);
  }

  template <class LHSDerived, class Number,
    typename = std::enable_if_t<
      std::is_arithmetic_v<Number>, Division<LHSDerived, ScalarFunction<Number>>>>
  inline
  auto
  operator/(const FunctionBase<LHSDerived>& lhs, Number rhs)
  {
    return Division(lhs, ScalarFunction(rhs));
  }

  template <class Number, class RHSDerived,
    typename = std::enable_if_t<
      std::is_arithmetic_v<Number>, Division<RHSDerived, ScalarFunction<Number>>>>
  inline
  auto
  operator/(Number lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return Division(ScalarFunction(lhs), rhs);
  }
}
#endif
