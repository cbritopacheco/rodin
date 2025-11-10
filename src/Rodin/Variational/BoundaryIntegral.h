/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file BoundaryIntegral.h
 * @brief Boundary integral classes for variational formulations.
 *
 * This file defines the BoundaryIntegral classes which represent integrals
 * over the domain boundary in finite element formulations. Boundary integrals
 * are essential for imposing natural (Neumann) boundary conditions and
 * boundary terms in weak formulations.
 *
 * ## Mathematical Foundation
 * Boundary integrals compute:
 * @f[
 *   \int_{\partial\Omega} f(x) \, ds
 * @f]
 * where @f$ \partial\Omega @f$ is the domain boundary and @f$ ds @f$ is the
 * surface measure.
 *
 * ## Applications
 * - **Neumann boundary conditions**: @f$ \int_{\Gamma_N} g \cdot v \, ds @f$
 * - **Robin boundary conditions**: @f$ \int_{\Gamma_R} (\alpha u - \beta) v \, ds @f$
 * - **Weak boundary terms**: Integration by parts produces boundary integrals
 * - **Flux calculations**: Computing normal fluxes across boundaries
 *
 * ## Usage Example
 * ```cpp
 * // Neumann BC: ∫_Γ g·v ds
 * auto neumann = BoundaryIntegral(g, v).on(2);  // On boundary attribute 2
 * 
 * // Robin BC: ∫_Γ (αu - β)v ds
 * auto robin = BoundaryIntegral(alpha * u - beta, v).on(3);
 * ```
 */
#ifndef RODIN_VARIATIONAL_BOUNDARYINTEGRAL_H
#define RODIN_VARIATIONAL_BOUNDARYINTEGRAL_H

#include <utility>

#include "ForwardDecls.h"
#include "Rodin/Geometry/Region.h"

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

      /**
       * @brief Constructs boundary integral from trial and test operators.
       * @param[in] lhs Trial operator @f$ A(u) @f$
       * @param[in] rhs Test operator @f$ B(v) @f$
       *
       * Represents the boundary integral:
       * @f[
       *   \int_{\partial\Omega} A(u) : B(v) \, d\sigma
       * @f]
       */
      BoundaryIntegral(const LHSType& lhs, const RHSType& rhs)
        : BoundaryIntegral(Dot(lhs, rhs))
      {}

      /**
       * @brief Constructs boundary integral from dot product.
       * @param[in] prod Dot product of trial and test operators
       */
      BoundaryIntegral(const IntegrandType& prod)
        : Parent(prod)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Boundary integral to copy
       */
      BoundaryIntegral(const BoundaryIntegral& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Boundary integral to move from
       */
      BoundaryIntegral(BoundaryIntegral&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Gets the integration region.
       * @return Region::Boundary indicating integration over boundary
       *
       * Boundary integrals are computed over boundary faces
       * @f$ F \in \mathcal{B}_h @f$ of the mesh.
       */
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Boundary;
      }

      /**
       * @brief Creates a polymorphic copy of this boundary integral.
       * @return Pointer to a new copy
       */
      BoundaryIntegral* copy() const noexcept override
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

      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Boundary;
      }

      BoundaryIntegral* copy() const noexcept override
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

