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
    using Parent = FunctionBase<Division<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>;
    static_assert(
        std::is_convertible_v<FormLanguage::Traits<FunctionBase<RHSDerived>>::ResultType, Scalar>);
    public:
      Division(const FunctionBase<LHSDerived>& lhs, const FunctionBase<RHSDerived>& rhs)
        : m_lhs(lhs), m_rhs(rhs)
      {}

      Division(const Division& other);

      Division(Division&& other);

      constexpr
      RangeShape getRangeShape() const
      {
        return m_lhs.getRangeShape();
      }

      Division& traceOf(Geometry::Attribute attr) override;

      auto getValue(const Geometry::Point& p) const override
      {
        auto v = m_lhs->getValue(p);
        v /= m_rhs->getValue(p).scalar();
        return v;
      }

      Division* copy() const noexcept override
      {
        return new Division(*this);
      }

    private:
      FunctionBase<LHSDerived> m_lhs;
      FunctionBase<RHSDerived> m_rhs;
  };
  Division(const FunctionBase&, const FunctionBase&)
    -> Division<FunctionBase, FunctionBase>;

  Division<FunctionBase, FunctionBase>
  operator/(const FunctionBase& lhs, const FunctionBase& rhs);

  template <class T>
  std::enable_if_t<std::is_arithmetic_v<T>,
    Division<FunctionBase, FunctionBase>>
  operator/(const FunctionBase& lhs, T rhs)
  {
    return Division(lhs, ScalarFunction(rhs));
  }

  template <class T>
  std::enable_if_t<std::is_arithmetic_v<T>,
    Division<FunctionBase, FunctionBase>>
  operator/(T lhs, const FunctionBase& rhs)
  {
    return Division(ScalarFunction(lhs), rhs);
  }
}
#endif
