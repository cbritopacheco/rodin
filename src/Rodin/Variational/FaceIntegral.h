/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FACEINTEGRAL_H
#define RODIN_VARIATIONAL_FACEINTEGRAL_H

#include <cassert>
#include <set>
#include <utility>

#include "Rodin/FormLanguage/Base.h"

#include "ForwardDecls.h"
#include "QuadratureRule.h"

namespace Rodin::Variational
{
  /**
   * @defgroup FaceIntegralSpecializations FaceIntegral Template Specializations
   * @brief Template specializations of the FaceIntegral class.
   *
   * @see FaceIntegral
   */

  /**
   * @ingroup FaceIntegralSpecializations
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
  class FaceIntegral<Dot<
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

      FaceIntegral(const LHSType& lhs, const RHSType& rhs)
        : FaceIntegral(Dot(lhs, rhs))
      {}

      FaceIntegral(const IntegrandType& prod)
        : Parent(prod)
      {}

      FaceIntegral(const FaceIntegral& other)
        : Parent(other)
      {}

      FaceIntegral(FaceIntegral&& other)
        : Parent(std::move(other))
      {}

      inline Integrator::Region getRegion() const override
      {
        return Integrator::Region::Faces;
      }

      inline FaceIntegral* copy() const noexcept override
      {
        return new FaceIntegral(*this);
      }

    private:
      std::unique_ptr<IntegrandType> m_prod;
  };

  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  FaceIntegral(const Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>&)
    -> FaceIntegral<Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>;

  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  FaceIntegral(const ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>&, const ShapeFunctionBase<RHSDerived, TestFES, TestSpace>&)
    -> FaceIntegral<Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>;

  /**
   * @ingroup FaceIntegralSpecializations
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
  class FaceIntegral<ShapeFunctionBase<NestedDerived, FES, TestSpace>> final
    : public QuadratureRule<ShapeFunctionBase<NestedDerived, FES, TestSpace>>
  {
    public:
      using IntegrandType = ShapeFunctionBase<NestedDerived, FES, TestSpace>;

      using Parent = QuadratureRule<IntegrandType>;

      template <class LHSDerived, class RHSDerived>
      FaceIntegral(
          const FunctionBase<LHSDerived>& lhs,
          const ShapeFunctionBase<RHSDerived, FES, TestSpace>& rhs)
        : FaceIntegral(Dot(lhs, rhs))
      {}

      FaceIntegral(const IntegrandType& integrand)
        : Parent(integrand)
      {}

      FaceIntegral(const FaceIntegral& other)
        : Parent(other)
      {}

      FaceIntegral(FaceIntegral&& other)
        : Parent(std::move(other))
      {}

      inline Integrator::Region getRegion() const override
      {
        return Integrator::Region::Faces;
      }

      inline FaceIntegral* copy() const noexcept override
      {
        return new FaceIntegral(*this);
      }

    private:
      std::unique_ptr<IntegrandType> m_integrand;
  };

  template <class NestedDerived, class FES>
  FaceIntegral(const ShapeFunctionBase<NestedDerived, FES, TestSpace>&)
    -> FaceIntegral<ShapeFunctionBase<NestedDerived, FES, TestSpace>>;

  template <class LHSDerived, class RHSDerived, class FES>
  FaceIntegral(const FunctionBase<LHSDerived>&, const ShapeFunctionBase<RHSDerived, FES, TestSpace>&)
    -> FaceIntegral<ShapeFunctionBase<Dot<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, FES, TestSpace>>, FES, TestSpace>>;
}

#endif
