/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_MAX_H
#define RODIN_VARIATIONAL_MAX_H

/**
 * @file
 * @brief Maximum function operations.
 */

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
   * @brief Represents the maximum operation between two functions.
   *
   * This class represents the pointwise maximum of two functions:
   * @f[
   *    \text{Max}(f, g)(x) = \max(f(x), g(x))
   * @f]
   *
   * For scalar functions, this computes the maximum value at each point.
   * For vector or matrix functions, the operation is applied componentwise.
   *
   * Common applications include:
   * - Positive part functions: @f$ \max(f, 0) @f$
   * - Barrier methods and penalty functions
   * - Non-smooth optimization problems
   * - Contact and friction problems in mechanics
   *
   * @tparam LHSDerived Type of the left operand function
   * @tparam RHSDerived Type of the right operand function
   *
   * @see Min, FunctionBase
   */
  template <class LHSDerived, class RHSDerived>
  class Max<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>> final
    : public FunctionBase<Max<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using LHSType = FunctionBase<LHSDerived>;

      using RHSType = FunctionBase<RHSDerived>;

      using Parent = FunctionBase<Max<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>;

      /**
       * @brief Constructs maximum operator for two functions.
       * @param a First function
       * @param b Second function
       */
      Max(const LHSType& a, const RHSType& b)
        : m_lhs(a.copy()), m_rhs(b.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param other Max operator to copy
       */
      Max(const Max& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      /**
       * @brief Move constructor.
       * @param other Max operator to move from
       */
      Max(Max&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      /**
       * @brief Restricts evaluation to specified mesh attributes.
       * @param args Attribute arguments for trace restriction
       * @returns Reference to this object
       */
      template <class ... Args>
      constexpr
      Max& traceOf(const Args& ... args)
      {
        m_lhs->traceOf(args...);
        m_rhs->traceOf(args...);
        return *this;
      }

      /**
       * @brief Evaluates maximum at a point.
       * @param p Point at which to evaluate
       * @returns @f$ \max(f(p), g(p)) @f$
       */
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

      /**
       * @brief Constructs maximum operator for function and scalar.
       * @param a Function operand
       * @param b Scalar constant operand
       */
      constexpr
      Max(const LHSType& a, const RHSType& b)
        : m_lhs(a.copy()), m_rhs(b)
      {}

      /**
       * @brief Copy constructor.
       * @param other Max operator to copy
       */
      constexpr
      Max(const Max& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs)
      {}

      /**
       * @brief Move constructor.
       * @param other Max operator to move from
       */
      constexpr
      Max(Max&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      /**
       * @brief Restricts evaluation to specified mesh attributes.
       * @param args Attribute arguments for trace restriction
       * @returns Reference to this object
       */
      template <class ... Args>
      constexpr
      Max& traceOf(const Args& ... args)
      {
        m_lhs->traceOf(args...);
        return *this;
      }

      /**
       * @brief Evaluates maximum at a point.
       * @param p Point at which to evaluate
       * @returns @f$ \max(f(p), c) @f$ where @f$ c @f$ is the scalar constant
       */
      constexpr
      Real getValue(const Geometry::Point& p) const
      {
        const auto lhs = this->getLHS().getValue(p);
        const auto& rhs = this->getRHS();
        if (lhs < rhs)
          return rhs;
        else
          return lhs;
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

      /**
       * @brief Polymorphic copy.
       * @returns Pointer to a copy of this object
       */
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
