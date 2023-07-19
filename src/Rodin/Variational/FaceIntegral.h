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

  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  class FaceIntegral<Dot<
          ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
          ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>> final
    : public QuadratureRule<Dot<
          ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
          ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>
  {
    public:
      using LHS = ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>;
      using RHS = ShapeFunctionBase<RHSDerived, TestFES, TestSpace>;
      using Integrand = Dot<LHS, RHS>;
      using Parent = QuadratureRule<Integrand>;

      FaceIntegral(const LHS& lhs, const RHS& rhs)
        : FaceIntegral(Dot(lhs, rhs))
      {}

      FaceIntegral(const Integrand& prod)
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
      std::unique_ptr<Integrand> m_prod;
  };

  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  FaceIntegral(const Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>&)
    -> FaceIntegral<Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>;

  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  FaceIntegral(const ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>&, const ShapeFunctionBase<RHSDerived, TestFES, TestSpace>&)
    -> FaceIntegral<Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>;

  template <class NestedDerived, class FES>
  class FaceIntegral<ShapeFunctionBase<NestedDerived, FES, TestSpace>> final
    : public QuadratureRule<ShapeFunctionBase<NestedDerived, FES, TestSpace>>
  {
    public:
      using Integrand = ShapeFunctionBase<NestedDerived, FES, TestSpace>;
      using Parent = QuadratureRule<Integrand>;

      template <class LHSDerived, class RHSDerived>
      FaceIntegral(const FunctionBase<LHSDerived>& lhs, const ShapeFunctionBase<RHSDerived, FES, TestSpace>& rhs)
        : FaceIntegral(Dot(lhs, rhs))
      {}

      FaceIntegral(const Integrand& integrand)
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
      std::unique_ptr<Integrand> m_integrand;
  };

  template <class NestedDerived, class FES>
  FaceIntegral(const ShapeFunctionBase<NestedDerived, FES, TestSpace>&)
    -> FaceIntegral<ShapeFunctionBase<NestedDerived, FES, TestSpace>>;

  template <class LHSDerived, class RHSDerived, class FES>
  FaceIntegral(const FunctionBase<LHSDerived>&, const ShapeFunctionBase<RHSDerived, FES, TestSpace>&)
    -> FaceIntegral<ShapeFunctionBase<Dot<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, FES, TestSpace>>, FES, TestSpace>>;
}

#endif
