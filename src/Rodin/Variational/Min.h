/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_Min_H
#define RODIN_VARIATIONAL_Min_H

#include <cmath>

#include "ForwardDecls.h"
#include "Function.h"
#include "RealFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup MinSpecializations Min Template Specializations
   * @brief Template specializations of the Min class.
   * @see Min
   */

  /**
   * @ingroup MinSpecializations
   * @brief Represents the minimum between two arguments.
   */
  template <class LHSDerived, class RHSDerived>
  class Min<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>> final
    : public RealFunctionBase<Min<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using LHSType = FunctionBase<LHSDerived>;

      using RHSType = FunctionBase<RHSDerived>;

      using Parent = RealFunctionBase<Min<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>;

      Min(const LHSType& a, const RHSType& b)
        : m_lhs(a.copy()), m_rhs(b.copy())
      {}

      Min(const Min& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      Min(Min&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      inline
      constexpr
      Min& traceOf(Geometry::Attribute attrs)
      {
        m_lhs.traceOf(attrs);
        m_rhs.traceOf(attrs);
        return *this;
      }

      inline
      constexpr
      Real getValue(const Geometry::Point& p) const
      {
        const auto lhs = getLHS().getValue(p);
        const auto rhs = getRHS().getValue(p);
        if (lhs < rhs)
          return lhs;
        else
          return rhs;
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

      virtual Min* copy() const noexcept override
      {
        return new Min(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      std::unique_ptr<RHSType> m_rhs;
  };

  template <class LHSDerived, class RHSDerived>
  Min(const FunctionBase<LHSDerived>&, const FunctionBase<RHSDerived>&)
    -> Min<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>;

  template <class NestedDerived>
  class Min<FunctionBase<NestedDerived>, Real>
    : public RealFunctionBase<Min<FunctionBase<NestedDerived>, Real>>
  {
    public:
      using LHSType = FunctionBase<NestedDerived>;

      using RHSType = Real;

      using Parent = RealFunctionBase<Min<FunctionBase<NestedDerived>, RHSType>>;

      constexpr
      Min(const LHSType& a, const RHSType& b)
        : m_lhs(a.copy()), m_rhs(b)
      {}

      constexpr
      Min(const Min& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs)
      {}

      constexpr
      Min(Min&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      inline
      constexpr
      Min& traceOf(Geometry::Attribute attrs)
      {
        m_lhs.traceOf(attrs);
        return *this;
      }

      inline
      constexpr
      Real getValue(const Geometry::Point& p) const
      {
        const auto lhs = getLHS().getValue(p);
        const auto& rhs = getRHS();
        if (lhs < rhs)
          return lhs;
        else
          return rhs;
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
        return m_rhs;
      }

      virtual Min* copy() const noexcept override
      {
        return new Min(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      RHSType m_rhs;
  };

  template <class NestedDerived>
  Min(const FunctionBase<NestedDerived>&, Real) -> Min<FunctionBase<NestedDerived>, Real>;

  template <class NestedDerived>
  class Min<Real, FunctionBase<NestedDerived>>
    : public Min<FunctionBase<NestedDerived>, Real>
  {
    public:
      using LHSType = Real;

      using RHSType = FunctionBase<NestedDerived>;

      using Parent = Min<FunctionBase<NestedDerived>, Real>;

      constexpr
      Min(const LHSType& a, const RHSType& b)
        : Parent(b, a)
      {}

      constexpr
      Min(const Min& other)
        : Parent(other)
      {}

      constexpr
      Min(Min&& other)
        : Parent(std::move(other))
      {}

      virtual Min* copy() const noexcept override
      {
        return new Min(*this);
      }
  };

  template <class NestedDerived>
  Min(Real, const FunctionBase<NestedDerived>&) -> Min<Real, FunctionBase<NestedDerived>>;
}

#endif
