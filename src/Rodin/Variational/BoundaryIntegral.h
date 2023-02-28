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
#include <mfem.hpp>

#include "Rodin/FormLanguage/Base.h"


#include "GaussianQuadrature.h"

namespace Rodin::Variational
{
  /**
   * @defgroup BoundaryIntegralSpecializations BoundaryIntegral Template Specializations
   * @brief Template specializations of the BoundaryIntegral class.
   *
   * @see BoundaryIntegral
   */

  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  class BoundaryIntegral<Dot<
          ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
          ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>> final
    : public GaussianQuadrature<Dot<
          ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
          ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>
  {
    public:
      using LHS = ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>;
      using RHS = ShapeFunctionBase<RHSDerived, TestFES, TestSpace>;
      using Integrand = Dot<LHS, RHS>;
      using Parent = GaussianQuadrature<Integrand>;

      constexpr
      BoundaryIntegral(const LHS& lhs, const RHS& rhs)
        : BoundaryIntegral(Dot(lhs, rhs))
      {}

      constexpr
      BoundaryIntegral(const Integrand& prod)
        : Parent(prod)
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

    private:
      std::unique_ptr<Integrand> m_prod;
  };

  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  BoundaryIntegral(const Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>&)
    -> BoundaryIntegral<Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>;

  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  BoundaryIntegral(const ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>&, const ShapeFunctionBase<RHSDerived, TestFES, TestSpace>&)
    -> BoundaryIntegral<Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>;

  template <class NestedDerived, class FES>
  class BoundaryIntegral<ShapeFunctionBase<NestedDerived, FES, TestSpace>> final
    : public GaussianQuadrature<ShapeFunctionBase<NestedDerived, FES, TestSpace>>
  {
    public:
      using Integrand = ShapeFunctionBase<NestedDerived, FES, TestSpace>;
      using Parent = GaussianQuadrature<Integrand>;

      template <class LHSDerived, class RHSDerived>
      constexpr
      BoundaryIntegral(const FunctionBase<LHSDerived>& lhs, const ShapeFunctionBase<RHSDerived, FES, TestSpace>& rhs)
        : BoundaryIntegral(Dot(lhs, rhs))
      {}

      constexpr
      BoundaryIntegral(const Integrand& integrand)
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

    private:
      std::unique_ptr<Integrand> m_integrand;
  };

  template <class NestedDerived, class FES>
  BoundaryIntegral(const ShapeFunctionBase<NestedDerived, FES, TestSpace>&)
    -> BoundaryIntegral<ShapeFunctionBase<NestedDerived, FES, TestSpace>>;

  template <class LHSDerived, class RHSDerived, class FES>
  BoundaryIntegral(const FunctionBase<LHSDerived>&, const ShapeFunctionBase<RHSDerived, FES, TestSpace>&)
    -> BoundaryIntegral<ShapeFunctionBase<Dot<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, FES, TestSpace>>, FES, TestSpace>>;

  // template <class FES>
  // class BoundaryIntegral<Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>>
  //   : public Integral<Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>>
  // {
  //   public:
  //     using Parent =
  //       Integral<Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>>;
  //     using Integrand =
  //       Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>;
  //     using Parent::Parent;
  //     Integrator::Region getRegion() const override { return Integrator::Region::Boundary; }
  //     BoundaryIntegral* copy() const noexcept override { return new BoundaryIntegral(*this); }
  // };
  // template <class FES>
  // BoundaryIntegral(
  //     const Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>&,
  //     const ShapeFunction<FES, TestSpace>&)
  //   -> BoundaryIntegral<Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>>;
  // template <class FES>
  // BoundaryIntegral(
  //     const Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>&)
  //   -> BoundaryIntegral<Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>>;

  // template <class FES>
  // class BoundaryIntegral<Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>>
  //   : public Integral<Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>>
  // {
  //   public:
  //     using Parent =
  //       Integral<Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>>;
  //     using Integrand =
  //       Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>;
  //     using Parent::Parent;
  //     Integrator::Region getRegion() const override { return Integrator::Region::Boundary; }
  //     BoundaryIntegral* copy() const noexcept override { return new BoundaryIntegral(*this); }
  // };
  // template <class FES>
  // BoundaryIntegral(
  //     const Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>&,
  //     const Grad<ShapeFunction<FES, TestSpace>>&)
  //   -> BoundaryIntegral<Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>>;
  // template <class FES>
  // BoundaryIntegral(
  //     const Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>&)
  //   -> BoundaryIntegral<Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>>;

  // template <class FES>
  // class BoundaryIntegral<Dot<FunctionBase, ShapeFunction<FES, TestSpace>>>
  //   : public Integral<Dot<FunctionBase, ShapeFunction<FES, TestSpace>>>
  // {
  //   public:
  //     using Parent    = Integral<Dot<FunctionBase, ShapeFunction<FES, TestSpace>>>;
  //     using Integrand  = Dot<FunctionBase, ShapeFunction<FES, TestSpace>>;
  //     using Parent::Parent;
  //     Integrator::Region getRegion() const override { return Integrator::Region::Boundary; }
  //     BoundaryIntegral* copy() const noexcept override { return new BoundaryIntegral(*this); }
  // };
  // template <class FES>
  // BoundaryIntegral(const FunctionBase&, const ShapeFunction<FES, TestSpace>&)
  //   -> BoundaryIntegral<Dot<FunctionBase, ShapeFunction<FES, TestSpace>>>;
  // template <class FES>
  // BoundaryIntegral(const Dot<FunctionBase, ShapeFunction<FES, TestSpace>>&)
  //   -> BoundaryIntegral<Dot<FunctionBase, ShapeFunction<FES, TestSpace>>>;
}

#endif

