/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_INTERFACEINTEGRAL_H
#define RODIN_VARIATIONAL_INTERFACEINTEGRAL_H

#include <cassert>
#include <set>
#include <utility>

#include "Rodin/FormLanguage/Base.h"

#include "ForwardDecls.h"
#include "QuadratureRule.h"

namespace Rodin::Variational
{
  /**
   * @defgroup InterfaceIntegralSpecializations InterfaceIntegral Template Specializations
   * @brief Template specializations of the InterfaceIntegral class.
   *
   * @see InterfaceIntegral
   */

  /**
   * @ingroup InterfaceIntegralSpecializations
   * @brief Integration of the dot product of a trial and test operators.
   *
   * Given two operators defined over trial and test spaces @f$ U_h
   * @f$ and @f$ V_h @f$,
   * @f[
   *   A : U_h \rightarrow \mathbb{R}^{p \times q}, \quad B : V_h \rightarrow \mathbb{R}^{p \times q},
   * @f]
   * this class represents the integral of their dot product:
   * @f[
   *   \int_{\mathcal{F}_h} A(u) : B(v) \ d\sigma(x) \ .
   * @f]
   */
  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  class InterfaceIntegral<Dot<
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

      InterfaceIntegral(const LHSType& lhs, const RHSType& rhs)
        : InterfaceIntegral(Dot(lhs, rhs))
      {}

      InterfaceIntegral(const IntegrandType& prod)
        : Parent(prod)
      {}

      InterfaceIntegral(const InterfaceIntegral& other)
        : Parent(other)
      {}

      InterfaceIntegral(InterfaceIntegral&& other)
        : Parent(std::move(other))
      {}

      inline Integrator::Region getRegion() const override
      {
        return Integrator::Region::Interface;
      }

      inline InterfaceIntegral* copy() const noexcept override
      {
        return new InterfaceIntegral(*this);
      }
  };

  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  InterfaceIntegral(const Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>&)
    -> InterfaceIntegral<Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>;

  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  InterfaceIntegral(const ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>&, const ShapeFunctionBase<RHSDerived, TestFES, TestSpace>&)
    -> InterfaceIntegral<Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>;

  /**
   * @ingroup InterfaceIntegralSpecializations
   * @brief Integration of a test operator.
   *
   * Given an operator defined over a test space @f$ V_h @f$
   * @f[
   *   A : V_h \rightarrow \mathbb{R},
   * @f]
   * this class will represent its integral
   * @f[
   *   \int_{\mathcal{F}_h} A(v) \ d\sigma(x) \ .
   * @f]
   */
  template <class NestedDerived, class FES>
  class InterfaceIntegral<ShapeFunctionBase<NestedDerived, FES, TestSpace>> final
    : public QuadratureRule<ShapeFunctionBase<NestedDerived, FES, TestSpace>>
  {
    public:
      using IntegrandType = ShapeFunctionBase<NestedDerived, FES, TestSpace>;

      using Parent = QuadratureRule<IntegrandType>;

      template <class LHSDerived, class RHSDerived>
      constexpr
      InterfaceIntegral(
          const FunctionBase<LHSDerived>& lhs,
          const ShapeFunctionBase<RHSDerived, FES, TestSpace>& rhs)
        : InterfaceIntegral(Dot(lhs, rhs))
      {}

      constexpr
      InterfaceIntegral(const IntegrandType& integrand)
        : Parent(integrand)
      {}

      constexpr
      InterfaceIntegral(const InterfaceIntegral& other)
        : Parent(other)
      {}

      constexpr
      InterfaceIntegral(InterfaceIntegral&& other)
        : Parent(std::move(other))
      {}

      inline Integrator::Region getRegion() const override
      {
        return Integrator::Region::Interface;
      }

      inline InterfaceIntegral* copy() const noexcept override
      {
        return new InterfaceIntegral(*this);
      }
  };

  template <class NestedDerived, class FES>
  InterfaceIntegral(const ShapeFunctionBase<NestedDerived, FES, TestSpace>&)
    -> InterfaceIntegral<ShapeFunctionBase<NestedDerived, FES, TestSpace>>;

  template <class LHSDerived, class RHSDerived, class FES>
  InterfaceIntegral(
      const FunctionBase<LHSDerived>&,
      const ShapeFunctionBase<RHSDerived, FES, TestSpace>&)
    -> InterfaceIntegral<
        ShapeFunctionBase<Dot<
          FunctionBase<LHSDerived>,
          ShapeFunctionBase<RHSDerived, FES, TestSpace>>, FES, TestSpace>>;
}

#endif

