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
#include "ScalarFunction.h"

namespace Rodin::Variational
{
  /**
   * @brief Represents the Minimum between two arguments.
   */
  template <class LHSDerived, class RHSDerived>
  class Min<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>> final
    : public ScalarFunctionBase<Min<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using LHS = FunctionBase<LHSDerived>;
      using RHS = FunctionBase<RHSDerived>;
      using Parent = ScalarFunctionBase<Min<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>;

      constexpr
      Min(const LHS& a, const RHS& b)
        : m_lhs(a), m_rhs(b)
      {}

      constexpr
      Min(const Min& other)
        : Parent(other),
          m_lhs(other.m_lhs), m_rhs(other.m_rhs)
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
        m_rhs.traceOf(attrs);
        return *this;
      }

      inline
      constexpr
      Scalar getValue(const Geometry::Point& p) const
      {
        return std::min(Scalar(m_lhs.getValue(p)), Scalar(m_rhs.getValue(p).scalar()));
      }

    private:
      LHS m_lhs;
      RHS m_rhs;
  };

  template <class LHSDerived, class RHSDerived>
  Min(const FunctionBase<LHSDerived>&, const FunctionBase<RHSDerived>&)
    -> Min<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>;

  template <class NestedDerived>
  class Min<FunctionBase<NestedDerived>, Scalar>
    : public ScalarFunctionBase<Min<FunctionBase<NestedDerived>, Scalar>>
  {
    public:
      using LHS = FunctionBase<NestedDerived>;
      using RHS = Scalar;
      using Parent = ScalarFunctionBase<Min<FunctionBase<NestedDerived>, RHS>>;

      constexpr
      Min(const LHS& a, const RHS& b)
        : m_lhs(a), m_rhs(b)
      {}

      constexpr
      Min(const Min& other)
        : Parent(other),
          m_lhs(other.m_lhs), m_rhs(other.m_rhs)
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
      Scalar getValue(const Geometry::Point& p) const
      {
        return std::min(Scalar(m_lhs.getValue(p)), m_rhs);
      }

    private:
      LHS m_lhs;
      RHS m_rhs;
  };

  template <class NestedDerived>
  Min(const FunctionBase<NestedDerived>&, Scalar) -> Min<FunctionBase<NestedDerived>, Scalar>;

  template <class NestedDerived>
  class Min<Scalar, FunctionBase<NestedDerived>>
    : public Min<FunctionBase<NestedDerived>, Scalar>
  {
    public:
      using LHS = Scalar;
      using RHS = FunctionBase<NestedDerived>;
      using Parent = Min<FunctionBase<NestedDerived>, Scalar>;

      constexpr
      Min(const LHS& a, const RHS& b)
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
  };

  template <class NestedDerived>
  Min(Scalar, const FunctionBase<NestedDerived>&) -> Min<Scalar, FunctionBase<NestedDerived>>;
}

#endif
