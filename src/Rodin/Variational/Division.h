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
      using Parent = FunctionBase<Division<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>;
      using LHS = FunctionBase<LHSDerived>;
      using RHS = FunctionBase<RHSDerived>;
      using LHSRange = typename FormLanguage::Traits<LHS>::RangeType;
      using RHSRange = typename FormLanguage::Traits<RHS>::RangeType;

      static_assert(FormLanguage::IsScalarRange<RHSRange>::Value);

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

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return getLHS().getRangeShape();
      }

      inline
      constexpr
      Division& traceOf(Geometry::Attribute attr)
      {
        getLHS().traceOf(attr);
        getRHS().traceOf(attr);
        return *this;
      }

      inline
      constexpr
      const LHS& getLHS() const
      {
        assert(m_lhs);
        return *m_lhs;
      }

      inline
      constexpr
      const RHS& getRHS() const
      {
        assert(m_rhs);
        return *m_rhs;
      }

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return this->object(getLHS().getValue(p)) / this->object(getRHS().getValue(p));
      }

      inline
      constexpr
      void getValueByReference(Math::Vector& res, const Geometry::Point& p) const
      {
        static_assert(FormLanguage::IsVectorRange<LHSRange>::Value);
        getLHS().getValue(res, p);
        res /= getRHS().getValue(p);
      }

      inline
      constexpr
      void getValueByReference(Math::Matrix& res, const Geometry::Point& p) const
      {
        static_assert(FormLanguage::IsMatrixRange<LHSRange>::Value);
        getLHS().getValue(res, p);
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
