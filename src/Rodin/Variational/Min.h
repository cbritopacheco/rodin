/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_Min_H
#define RODIN_VARIATIONAL_Min_H

/**
 * @file
 * @brief Minimum function operations.
 */

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
   * @brief Represents the minimum operation between two functions.
   *
   * This class represents the pointwise minimum of two functions:
   * @f[
   *    \text{Min}(f, g)(x) = \min(f(x), g(x))
   * @f]
   *
   * For scalar functions, this computes the minimum value at each point.
   * For vector or matrix functions, the operation is applied componentwise.
   *
   * Common applications include:
   * - Truncation and clamping operations
   * - Complementarity conditions
   * - Non-smooth optimization problems
   *
   * @tparam LHSDerived Type of the left operand function
   * @tparam RHSDerived Type of the right operand function
   *
   * @see Max, FunctionBase
   */
  template <class LHSDerived, class RHSDerived>
  class Min<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>> final
    : public FunctionBase<Min<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using LHSType = FunctionBase<LHSDerived>;

      using RHSType = FunctionBase<RHSDerived>;

      using Parent = FunctionBase<Min<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>;

      /**
       * @brief Constructs minimum operator for two functions.
       * @param a First function
       * @param b Second function
       */
      Min(const LHSType& a, const RHSType& b)
        : m_lhs(a.copy()), m_rhs(b.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param other Min operator to copy
       */
      Min(const Min& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      /**
       * @brief Move constructor.
       * @param other Min operator to move from
       */
      Min(Min&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      /**
       * @brief Restricts evaluation to specified mesh attributes.
       * @param attrs Mesh attributes for trace restriction
       * @returns Reference to this object
       */
      constexpr
      Min& traceOf(Geometry::Attribute attrs)
      {
        m_lhs.traceOf(attrs);
        m_rhs.traceOf(attrs);
        return *this;
      }

      /**
       * @brief Evaluates minimum at a point.
       * @param p Point at which to evaluate
       * @returns @f$ \min(f(p), g(p)) @f$
       */
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        const auto lhs = this->getLHS().getValue(p);
        const auto rhs = this->getRHS().getValue(p);
        if (lhs < rhs)
          return lhs;
        else
          return rhs;
      }

      Optional<size_t> getOrder(const Geometry::Polytope& polytope) const noexcept
      {
        const auto lo = getLHS().getOrder(polytope);
        const auto ro = getRHS().getOrder(polytope);
        if (lo && ro && *lo == 0 && *ro == 0)
          return size_t{0};
        return std::nullopt;
      }

      /**
       * @brief Gets the left operand function.
       * @returns Reference to first function
       */
      const auto& getLHS() const
      {
        assert(m_lhs);
        return *m_lhs;
      }

      /**
       * @brief Gets the right operand function.
       * @returns Reference to second function
       */
      const auto& getRHS() const
      {
        assert(m_rhs);
        return *m_rhs;
      }

      /**
       * @brief Polymorphic copy.
       * @returns Pointer to a copy of this object
       */
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

      using ScalarType = Real;

      using RHSType = ScalarType;

      using Parent = RealFunctionBase<Min<FunctionBase<NestedDerived>, RHSType>>;

      /**
       * @brief Constructs minimum operator for function and scalar.
       * @param a Function operand
       * @param b Scalar constant operand
       */
      constexpr
      Min(const LHSType& a, const RHSType& b)
        : m_lhs(a.copy()), m_rhs(b)
      {}

      /**
       * @brief Copy constructor.
       * @param other Min operator to copy
       */
      constexpr
      Min(const Min& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs)
      {}

      /**
       * @brief Move constructor.
       * @param other Min operator to move from
       */
      constexpr
      Min(Min&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      /**
       * @brief Restricts evaluation to specified mesh attributes.
       * @param attrs Mesh attributes for trace restriction
       * @returns Reference to this object
       */
      constexpr
      Min& traceOf(Geometry::Attribute attrs)
      {
        m_lhs.traceOf(attrs);
        return *this;
      }

      /**
       * @brief Evaluates minimum at a point.
       * @param p Point at which to evaluate
       * @returns @f$ \min(f(p), c) @f$ where @f$ c @f$ is the scalar constant
       */
      constexpr
      ScalarType getValue(const Geometry::Point& p) const
      {
        const auto lhs = getLHS().getValue(p);
        const auto& rhs = getRHS();
        if (lhs < rhs)
          return lhs;
        else
          return rhs;
      }

      /**
       * @brief Gets the function operand.
       * @returns Reference to the function
       */
      const auto& getLHS() const
      {
        assert(m_lhs);
        return *m_lhs;
      }

      /**
       * @brief Gets the scalar operand.
       * @returns Scalar constant value
       */
      const auto& getRHS() const
      {
        return m_rhs;
      }

      Optional<size_t> getOrder(const Geometry::Polytope& polytope) const noexcept
      {
        const auto lo = getLHS().getOrder(polytope);
        const auto ro = getRHS().getOrder(polytope);
        if (lo && ro && *lo == 0 && *ro == 0)
          return size_t{0};
        return std::nullopt;
      }

      /**
       * @brief Polymorphic copy.
       * @returns Pointer to a copy of this object
       */
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
