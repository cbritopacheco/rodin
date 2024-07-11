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
   * @ingroup DivisionSpecializations
   * @brief Division of a FunctionBase by a FunctionBase.
   */
  template <class LHSDerived, class RHSDerived>
  class Division<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>
    : public FunctionBase<Division<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using LHSType = FunctionBase<LHSDerived>;

      using RHSType = FunctionBase<RHSDerived>;

      using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;

      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      using Parent =
        FunctionBase<Division<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>;

      Division(const FunctionBase<LHSDerived>& lhs, const FunctionBase<RHSDerived>& rhs)
        : m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {}


      Division(const Division& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      Division(Division&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      constexpr
      RangeShape getRangeShape() const
      {
        return getLHS().getRangeShape();
      }

      constexpr
      Division& traceOf(Geometry::Attribute attr)
      {
        getLHS().traceOf(attr);
        getRHS().traceOf(attr);
        return *this;
      }

      constexpr
      const LHSType& getLHS() const
      {
        assert(m_lhs);
        return *m_lhs;
      }

      constexpr
      const RHSType& getRHS() const
      {
        assert(m_rhs);
        return *m_rhs;
      }

      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return this->object(getLHS().getValue(p)) / this->object(getRHS().getValue(p));
      }

      template <class T>
      constexpr
      void getValue(T& res, const Geometry::Point& p) const
      {
        getLHS().getDerived().getValue(res, p);
        res /= getRHS().getValue(p);
      }

      inline
      Division* copy() const noexcept final override
      {
        return new Division(*this);
      }

    private:
      std::unique_ptr<FunctionBase<LHSDerived>> m_lhs;
      std::unique_ptr<FunctionBase<RHSDerived>> m_rhs;
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
      std::is_arithmetic_v<Number>, Division<LHSDerived, RealFunction<Number>>>>
  inline
  auto
  operator/(const FunctionBase<LHSDerived>& lhs, Number rhs)
  {
    return Division(lhs, RealFunction(rhs));
  }

  template <class Number, class RHSDerived,
    typename = std::enable_if_t<
      std::is_arithmetic_v<Number>, Division<RHSDerived, RealFunction<Number>>>>
  inline
  auto
  operator/(Number lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return Division(RealFunction(lhs), rhs);
  }
}
#endif
