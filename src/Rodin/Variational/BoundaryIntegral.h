/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_BOUNDARYINTEGRAL_H
#define RODIN_VARIATIONAL_BOUNDARYINTEGRAL_H

/**
 * @file
 * @brief Boundary integrals for variational formulations.
 *
 * This file defines boundary integrals over @f$ \partial\Omega @f$, the boundary
 * of a domain. Boundary integrals are fundamental for incorporating boundary
 * conditions and computing boundary fluxes in finite element methods.
 *
 * # Mathematical Foundation
 *
 * A boundary integral evaluates an expression over the domain boundary:
 * @f[
 *   \int_{\partial\Omega} f(x) \, ds
 * @f]
 * where @f$ ds @f$ denotes the surface measure on @f$ \partial\Omega @f$.
 *
 * For bilinear and linear forms:
 * - **Bilinear**: @f$ b(u,v) = \int_{\partial\Omega} a(u,v) \, ds @f$
 * - **Linear**: @f$ l(v) = \int_{\partial\Omega} f \cdot v \, ds @f$
 *
 * ## Common Applications
 * - **Neumann BC**: @f$ -\nabla u \cdot \mathbf{n} = g @f$ on @f$ \Gamma_N @f$
 * - **Robin BC**: @f$ \alpha u + \beta \nabla u \cdot \mathbf{n} = g @f$ on @f$ \Gamma_R @f$
 * - **Flux computation**: @f$ \int_{\partial\Omega} \mathbf{u} \cdot \mathbf{n} \, ds @f$
 * - **Boundary energy**: @f$ \int_{\partial\Omega} |u|^2 \, ds @f$
 *
 * # Usage Examples
 *
 * ## Neumann Boundary Condition
 * @code{.cpp}
 * using namespace Rodin;
 * using namespace Rodin::Variational;
 *
 * P1 Vh(mesh);
 * TrialFunction u(Vh);
 * TestFunction  v(Vh);
 *
 * // Neumann BC: -grad u · n = g on boundary
 * auto g = ScalarFunction([](auto& x) {
 *   return std::sin(x(0));
 * });
 *
 * Problem problem(u, v);
 * problem = Integral(Grad(u), Grad(v))
 *         - BoundaryIntegral(g * v);  // Neumann term
 *
 * problem.solve(solver);
 * @endcode
 *
 * ## Robin Boundary Condition
 * @code{.cpp}
 * // Robin BC: alpha*u + beta*grad(u)·n = g
 * double alpha = 1.0, beta = 0.1;
 *
 * auto g = ScalarFunction([](auto& x) { return 1.0; });
 *
 * Problem problem(u, v);
 * problem = Integral(Grad(u), Grad(v))
 *         + BoundaryIntegral(alpha * u * v)  // Robin term (u part)
 *         - BoundaryIntegral(g * v);          // Robin RHS
 *
 * problem.solve(solver);
 * @endcode
 *
 * ## Selective Boundary Regions (By Attribute)
 * @code{.cpp}
 * // Apply Neumann BC only on specific boundary regions
 * std::set<int> neumannAttrs = {1, 2};  // Boundary markers
 *
 * Problem problem(u, v);
 * problem = Integral(Grad(u), Grad(v))
 *         - BoundaryIntegral(g * v).over(neumannAttrs);
 *
 * problem.solve(solver);
 * @endcode
 *
 * ## Computing Boundary Flux
 * @code{.cpp}
 * // Compute flux through boundary
 * P1 Vh(mesh, 2);  // Vector field
 * GridFunction<P1> u(Vh);
 *
 * auto n = BoundaryNormal();
 *
 * // Total flux = int(u·n) ds
 * double total_flux = BoundaryIntegral(Dot(u, n)).compute();
 *
 * std::cout << "Total outward flux = " << total_flux << std::endl;
 * @endcode
 *
 * ## Surface Energy Integral
 * @code{.cpp}
 * // Compute L2 norm on boundary
 * P1 Vh(mesh);
 * GridFunction<P1> u(Vh);
 *
 * // ||u||_{L^2(\partial\Omega)}^2 = int_{\partial\Omega} |u|^2 ds
 * double boundary_norm_sq = BoundaryIntegral(u * u).compute();
 * double boundary_norm = std::sqrt(boundary_norm_sq);
 * @endcode
 *
 * ## Mixed Boundary Conditions
 * @code{.cpp}
 * // Dirichlet on some boundaries, Neumann on others
 * std::set<int> dirichletAttrs = {1, 2};
 * std::set<int> neumannAttrs = {3, 4};
 *
 * Problem problem(u, v);
 * problem = Integral(Grad(u), Grad(v))
 *         - BoundaryIntegral(g * v).over(neumannAttrs)
 *         + DirichletBC(u, 0.0).on(dirichletAttrs);
 *
 * problem.solve(solver);
 * @endcode
 *
 * ## Convective Boundary Condition (Heat Transfer)
 * @code{.cpp}
 * // Heat transfer with convection: k*grad(u)·n + h*(u - u_amb) = 0
 * double h = 10.0;       // Convection coefficient
 * double u_amb = 20.0;   // Ambient temperature
 *
 * Problem problem(u, v);
 * problem = Integral(Grad(u), Grad(v))
 *         + BoundaryIntegral(h * u * v)           // Convection term
 *         - BoundaryIntegral(h * u_amb * v);      // Ambient contribution
 *
 * problem.solve(solver);
 * @endcode
 *
 * ## Variational Derivative (Shape Optimization)
 * @code{.cpp}
 * // Compute boundary sensitivity
 * auto n = BoundaryNormal();
 *
 * // Shape derivative involves boundary integral
 * double sensitivity = BoundaryIntegral(
 *   Dot(Grad(u), Grad(u)) * Dot(perturbation, n)
 * ).compute();
 * @endcode
 *
 * @see Integral, FaceIntegral, InterfaceIntegral
 * @see BoundaryNormal, DirichletBC, Problem
 */

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

      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Boundary;
      }

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

