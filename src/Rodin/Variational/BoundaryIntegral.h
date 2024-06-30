/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_BOUNDARYINTEGRAL_H
#define RODIN_VARIATIONAL_BOUNDARYINTEGRAL_H

#include <cassert>
#include <set>
#include <utility>

#include "Rodin/FormLanguage/Base.h"

#include "ForwardDecls.h"
#include "QuadratureRule.h"

namespace Rodin::Variational
{
  /**
   * @defgroup BoundaryIntegralSpecializations BoundaryIntegral Template Specializations
   * @brief Template specializations of the BoundaryIntegral class.
   *
   * @see BoundaryIntegral
   */

  /**
   * @ingroup BoundaryIntegralSpecializations
   * @brief Integration of the dot product of a trial and test operators.
   *
   * Given two operators defined over trial and test spaces @f$ U_h
   * @f$ and @f$ V_h @f$,
   * @f[
   *   A : U_h \rightarrow \mathbb{R}^{p \times q}, \quad B : V_h \rightarrow \mathbb{R}^{p \times q},
   * @f]
   * this class represents the integral of their dot product:
   * @f[
   *   \int_{\mathcal{B}_h} A(u) : B(v) \ d\sigma(x) \ .
   * @f]
   */
  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  class BoundaryIntegral<Dot<
          ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
          ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>> final
    : public QuadratureRule<Dot<
          ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
          ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>
  {
    public:
      using LHSType = ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>;

      using RHSType = ShapeFunctionBase<RHSDerived, TestFES, TestSpace>;

      using IntegrandType = Dot<LHSType, RHSType>;

      using Parent = QuadratureRule<IntegrandType>;

      BoundaryIntegral(const LHSType& lhs, const RHSType& rhs)
        : BoundaryIntegral(Dot(lhs, rhs))
      {}

      BoundaryIntegral(const IntegrandType& prod)
        : Parent(prod)
      {}

      BoundaryIntegral(const BoundaryIntegral& other)
        : Parent(other)
      {}

      BoundaryIntegral(BoundaryIntegral&& other)
        : Parent(std::move(other))
      {}

      inline Integrator::Region getRegion() const override
      {
        return Integrator::Region::Boundary;
      }

      inline BoundaryIntegral* copy() const noexcept override
      {
        return new BoundaryIntegral(*this);
      }
  };

  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  BoundaryIntegral(const Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>&)
    -> BoundaryIntegral<Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>;

  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  BoundaryIntegral(const ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>&, const ShapeFunctionBase<RHSDerived, TestFES, TestSpace>&)
    -> BoundaryIntegral<Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>;

  /**
   * @ingroup BoundaryIntegralSpecializations
   * @brief Integration of a test operator.
   *
   * Given an operator defined over a test space @f$ V_h @f$
   * @f[
   *   A : V_h \rightarrow \mathbb{R},
   * @f]
   * this class will represent its integral
   * @f[
   *   \int_{\mathcal{B}_h} A(v) \ d\sigma(x) \ .
   * @f]
   */
  template <class NestedDerived, class FES>
  class BoundaryIntegral<ShapeFunctionBase<NestedDerived, FES, TestSpace>> final
    : public QuadratureRule<ShapeFunctionBase<NestedDerived, FES, TestSpace>>
  {
    public:
      using IntegrandType = ShapeFunctionBase<NestedDerived, FES, TestSpace>;

      using Parent = QuadratureRule<IntegrandType>;

      template <class LHSDerived, class RHSDerived>
      constexpr
      BoundaryIntegral(
          const FunctionBase<LHSDerived>& lhs,
          const ShapeFunctionBase<RHSDerived, FES, TestSpace>& rhs)
        : BoundaryIntegral(Dot(lhs, rhs))
      {}

      constexpr
      BoundaryIntegral(const IntegrandType& integrand)
        : Parent(integrand)
      {}

      constexpr
      BoundaryIntegral(const BoundaryIntegral& other)
        : Parent(other)
      {}

      constexpr
      BoundaryIntegral(BoundaryIntegral&& other)
        : Parent(std::move(other))
      {}

      inline Integrator::Region getRegion() const override
      {
        return Integrator::Region::Boundary;
      }

      inline BoundaryIntegral* copy() const noexcept override
      {
        return new BoundaryIntegral(*this);
      }
  };

  template <class NestedDerived, class FES>
  BoundaryIntegral(const ShapeFunctionBase<NestedDerived, FES, TestSpace>&)
    -> BoundaryIntegral<ShapeFunctionBase<NestedDerived, FES, TestSpace>>;

  template <class LHSDerived, class RHSDerived, class FES>
  BoundaryIntegral(
      const FunctionBase<LHSDerived>&,
      const ShapeFunctionBase<RHSDerived, FES, TestSpace>&)
    -> BoundaryIntegral<
        ShapeFunctionBase<Dot<
          FunctionBase<LHSDerived>,
          ShapeFunctionBase<RHSDerived, FES, TestSpace>>, FES, TestSpace>>;
}

#endif

