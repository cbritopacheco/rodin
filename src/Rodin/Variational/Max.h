/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_MAX_H
#define RODIN_VARIATIONAL_MAX_H

#include <cmath>

#include "ForwardDecls.h"
#include "Function.h"
#include "RealFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup MaxSpecializations Max Template Specializations
   * @brief Template specializations of the Max class.
   * @see Max
   */

  /**
   * @ingroup MaxSpecializations
   */
  template <class LHSDerived, class RHSDerived>
  class Max<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>> final
    : public FunctionBase<Max<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using LHSType = FunctionBase<LHSDerived>;

      using RHSType = FunctionBase<RHSDerived>;

      using Parent = FunctionBase<Max<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>;

      Max(const LHSType& a, const RHSType& b)
        : m_lhs(a.copy()), m_rhs(b.copy())
      {}

      Max(const Max& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      Max(Max&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      inline
      constexpr
      Max& traceOf(Geometry::Attribute attrs)
      {
        m_lhs.traceOf(attrs);
        m_rhs.traceOf(attrs);
        return *this;
      }

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        const auto lhs = getLHS().getValue(p);
        const auto rhs = getRHS().getValue(p);
        if (lhs < rhs)
          return rhs;
        else
          return lhs;
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

      virtual Max* copy() const noexcept override
      {
        return new Max(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      std::unique_ptr<RHSType> m_rhs;
  };

  template <class LHSDerived, class RHSDerived>
  Max(const FunctionBase<LHSDerived>&, const FunctionBase<RHSDerived>&)
    -> Max<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>;

  /**
   * @ingroup MaxSpecializations
   */
  template <class NestedDerived>
  class Max<FunctionBase<NestedDerived>, Real>
    : public RealFunctionBase<Max<FunctionBase<NestedDerived>, Real>>
  {
    public:
      using LHSType = FunctionBase<NestedDerived>;

      using RHSType = Real;

      using Parent = RealFunctionBase<Max<FunctionBase<NestedDerived>, RHSType>>;

      constexpr
      Max(const LHSType& a, const RHSType& b)
        : m_lhs(a.copy()), m_rhs(b)
      {}

      constexpr
      Max(const Max& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs)
      {}

      constexpr
      Max(Max&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      inline
      constexpr
      Max& traceOf(Geometry::Attribute attrs)
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
          return rhs;
        else
          return lhs;
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

      virtual Max* copy() const noexcept override
      {
        return new Max(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      RHSType m_rhs;
  };

  template <class NestedDerived>
  Max(const FunctionBase<NestedDerived>&, Real) -> Max<FunctionBase<NestedDerived>, Real>;

  template <class NestedDerived>
  class Max<Real, FunctionBase<NestedDerived>>
    : public Max<FunctionBase<NestedDerived>, Real>
  {
    public:
      using LHSType = Real;

      using RHSType = FunctionBase<NestedDerived>;

      using Parent = Max<FunctionBase<NestedDerived>, Real>;

      constexpr
      Max(const LHSType& a, const RHSType& b)
        : Parent(b, a)
      {}

      constexpr
      Max(const Max& other)
        : Parent(other)
      {}

      constexpr
      Max(Max&& other)
        : Parent(std::move(other))
      {}

      virtual Max* copy() const noexcept override
      {
        return new Max(*this);
      }
  };

  template <class NestedDerived>
  Max(Real, const FunctionBase<NestedDerived>&) -> Max<Real, FunctionBase<NestedDerived>>;
}

#endif
