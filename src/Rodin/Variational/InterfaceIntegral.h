/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file InterfaceIntegral.h
 * @brief Interface integral classes for interior element faces.
 *
 * This file defines the InterfaceIntegral classes which represent integrals
 * over interior mesh faces (interfaces between elements) in discontinuous
 * Galerkin formulations. Interface integrals specifically exclude boundary
 * faces, focusing on element-to-element coupling.
 *
 * ## Mathematical Foundation
 * Interface integrals compute:
 * @f[
 *   \int_{\Gamma_h} f(x) \, ds
 * @f]
 * where @f$ \Gamma_h @f$ is the set of interior faces (excluding the domain
 * boundary @f$ \partial\Omega @f$).
 *
 * ## Discontinuous Galerkin Context
 * Interface integrals are fundamental to DG methods, enabling:
 * - **Element coupling**: Connect discontinuous solutions across elements
 * - **Flux consistency**: Ensure conservative numerical fluxes
 * - **Stability**: Add necessary stabilization terms
 *
 * ## Interior Penalty Method
 * A typical IP-DG formulation includes:
 * @f[
 *   a(u,v) = \int_\Omega \nabla u \cdot \nabla v \, dx
 *          - \int_{\Gamma_h} \{\!\{\nabla u\}\!\} \cdot [\![v]\!] \, ds
 *          - \int_{\Gamma_h} \{\!\{\nabla v\}\!\} \cdot [\![u]\!] \, ds
 *          + \int_{\Gamma_h} \frac{\sigma}{h} [\![u]\!] \cdot [\![v]\!] \, ds
 * @f]
 *
 * ## Usage Example
 * ```cpp
 * // Interior penalty terms (only on interior faces)
 * auto penalty = InterfaceIntegral(sigma / h * Jump(u), Jump(v));
 * auto consistency = InterfaceIntegral(Average(Grad(u)), Jump(v));
 * auto symmetry = InterfaceIntegral(Jump(u), Average(Grad(v)));
 * ```
 *
 * @see FaceIntegral, BoundaryIntegral, Jump, Average
 */
#ifndef RODIN_VARIATIONAL_INTERFACEINTEGRAL_H
#define RODIN_VARIATIONAL_INTERFACEINTEGRAL_H

#include <utility>

#include "Rodin/Geometry/Region.h"

#include "ForwardDecls.h"

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

      /**
       * @brief Constructs interface integral from trial and test operators.
       * @param lhs Trial operator (left-hand side)
       * @param rhs Test operator (right-hand side)
       * 
       * Creates the interior face integral @f$ \int_{\Gamma_h} A(u) : B(v) d\sigma @f$
       * over the set of interior faces @f$ \Gamma_h @f$ (excluding boundary).
       */
      InterfaceIntegral(const LHSType& lhs, const RHSType& rhs)
        : InterfaceIntegral(Dot(lhs, rhs))
      {}

      /**
       * @brief Constructs interface integral from dot product.
       * @param prod Dot product integrand
       */
      InterfaceIntegral(const IntegrandType& prod)
        : Parent(prod)
      {}

      /**
       * @brief Copy constructor.
       * @param other InterfaceIntegral to copy
       */
      InterfaceIntegral(const InterfaceIntegral& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other InterfaceIntegral to move
       */
      InterfaceIntegral(InterfaceIntegral&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Returns the integration region.
       * @return Region::Interface indicating integration over interior faces only
       * 
       * Interior faces are those shared between two elements, excluding
       * faces on the domain boundary @f$ \partial\Omega @f$.
       */
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Interface;
      }

      /**
       * @brief Creates a polymorphic copy.
       * @return Pointer to new InterfaceIntegral instance
       */
      InterfaceIntegral* copy() const noexcept override
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

      /**
       * @brief Constructs interface integral from function and test operator.
       * @tparam LHSDerived Type of left-hand side function
       * @tparam RHSDerived Type of right-hand side shape function
       * @param lhs Function operator (left-hand side)
       * @param rhs Test operator (right-hand side)
       * 
       * Creates the interior face integral @f$ \int_{\Gamma_h} f(x) \cdot A(v) d\sigma @f$.
       */
      template <class LHSDerived, class RHSDerived>
      constexpr
      InterfaceIntegral(
          const FunctionBase<LHSDerived>& lhs,
          const ShapeFunctionBase<RHSDerived, FES, TestSpace>& rhs)
        : InterfaceIntegral(Dot(lhs, rhs))
      {}

      /**
       * @brief Constructs interface integral from integrand.
       * @param integrand Test operator to integrate
       */
      constexpr
      InterfaceIntegral(const IntegrandType& integrand)
        : Parent(integrand)
      {}

      /**
       * @brief Copy constructor.
       * @param other InterfaceIntegral to copy
       */
      constexpr
      InterfaceIntegral(const InterfaceIntegral& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other InterfaceIntegral to move
       */
      constexpr
      InterfaceIntegral(InterfaceIntegral&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Returns the integration region.
       * @return Region::Interface indicating integration over interior faces only
       */
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Interface;
      }

      InterfaceIntegral* copy() const noexcept override
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

