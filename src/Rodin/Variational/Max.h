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
#include "ScalarFunction.h"

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
    : public ScalarFunctionBase<Max<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using LHS = FunctionBase<LHSDerived>;
      using RHS = FunctionBase<RHSDerived>;
      using Parent = ScalarFunctionBase<Max<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>;

      Max(const LHS& a, const RHS& b)
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
      Scalar getValue(const Geometry::Point& p) const
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
      std::unique_ptr<LHS> m_lhs;
      std::unique_ptr<RHS> m_rhs;
  };

  template <class LHSDerived, class RHSDerived>
  Max(const FunctionBase<LHSDerived>&, const FunctionBase<RHSDerived>&)
    -> Max<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>;

  /**
   * @ingroup MaxSpecializations
   */
  template <class NestedDerived>
  class Max<FunctionBase<NestedDerived>, Scalar>
    : public ScalarFunctionBase<Max<FunctionBase<NestedDerived>, Scalar>>
  {
    public:
      using LHS = FunctionBase<NestedDerived>;
      using RHS = Scalar;
      using Parent = ScalarFunctionBase<Max<FunctionBase<NestedDerived>, RHS>>;

      constexpr
      Max(const LHS& a, const RHS& b)
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
      Scalar getValue(const Geometry::Point& p) const
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
      std::unique_ptr<LHS> m_lhs;
      RHS m_rhs;
  };

  template <class NestedDerived>
  Max(const FunctionBase<NestedDerived>&, Scalar) -> Max<FunctionBase<NestedDerived>, Scalar>;

  template <class NestedDerived>
  class Max<Scalar, FunctionBase<NestedDerived>>
    : public Max<FunctionBase<NestedDerived>, Scalar>
  {
    public:
      using LHS = Scalar;
      using RHS = FunctionBase<NestedDerived>;
      using Parent = Max<FunctionBase<NestedDerived>, Scalar>;

      constexpr
      Max(const LHS& a, const RHS& b)
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
  Max(Scalar, const FunctionBase<NestedDerived>&) -> Max<Scalar, FunctionBase<NestedDerived>>;
}

#endif
