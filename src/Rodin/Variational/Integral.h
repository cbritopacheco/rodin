/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Integral.h
 * @brief Domain integral classes for variational formulations.
 *
 * This file defines the Integral classes which represent domain integrals
 * in finite element formulations. These integrals form the basis for defining
 * bilinear and linear forms in the weak formulation.
 *
 * ## Mathematical Foundation
 * Domain integrals compute:
 * @f[
 *   \int_\Omega f(x) \, dx
 * @f]
 * where @f$ \Omega @f$ is the computational domain and @f$ f @f$ is the integrand.
 *
 * ## Bilinear Form Integrals
 * For bilinear forms @f$ a(u,v) @f$:
 * @f[
 *   a(u,v) = \int_\Omega A(u) : B(v) \, dx
 * @f]
 * where @f$ A @f$ and @f$ B @f$ are operators on trial and test functions.
 *
 * ## Linear Form Integrals
 * For linear forms @f$ l(v) @f$:
 * @f[
 *   l(v) = \int_\Omega f \cdot v \, dx
 * @f]
 *
 * ## Usage Examples
 * ```cpp
 * // Stiffness matrix: ∫ ∇u·∇v dx
 * auto stiffness = Integral(Grad(u), Grad(v));
 * 
 * // Mass matrix: ∫ u·v dx
 * auto mass = Integral(u, v);
 * 
 * // Load vector: ∫ f·v dx
 * auto load = Integral(f, v);
 * ```
 */
#ifndef RODIN_VARIATIONAL_INTEGRAL_H
#define RODIN_VARIATIONAL_INTEGRAL_H

#include <cassert>
#include <utility>

#include "Rodin/Geometry/Region.h"

#include "Function.h"
#include "ForwardDecls.h"
#include "ShapeFunction.h"

namespace Rodin::Variational
{

  /**
   * @defgroup IntegralSpecializations Integral Template Specializations
   * @brief Template specializations of the Integral class.
   *
   * @see Integral
   */

  /**
   * @ingroup IntegralSpecializations
   * @brief Integration of the dot product of a trial and test operators.
   *
   * Given two operators defined over trial and test spaces @f$ U_h
   * @f$ and @f$ V_h @f$,
   * @f[
   *   A : U_h \rightarrow \mathbb{R}^{p \times q}, \quad B : V_h \rightarrow \mathbb{R}^{p \times q},
   * @f]
   * this class represents the integral of their dot product:
   * @f[
   *   \int_{\mathcal{T}_h} A(u) : B(v) \ dx \ .
   * @f]
   */
  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  class Integral<Dot<
          ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
          ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>> final
    : public QuadratureRule<Dot<
          ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
          ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>
  {
    public:
      /// Type of the left operand of the dot product
      using LHSType = ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>;

      /// Type of the right operand of the dot product
      using RHSType = ShapeFunctionBase<RHSDerived, TestFES, TestSpace>;

      /// Type of the integrand
      using IntegrandType = Dot<LHSType, RHSType>;

      /// Parent class
      using Parent = QuadratureRule<IntegrandType>;

      /**
       * @brief Integral of the dot product of trial and test operators
       *
       * Constructs an object representing the following integral:
       * @f[
       *   \int_\Omega A(u) : B(v) \ dx
       * @f]
       *
       * @param[in] lhs Trial operator @f$ A(u) @f$
       * @param[in] rhs Test operator @f$ B(v) @f$
       */
      Integral(const LHSType& lhs, const RHSType& rhs)
        : Integral(Dot(lhs, rhs))
      {}

      /**
       * @brief Integral of the dot product of trial and test operators
       *
       * Constructs the object representing the following integral:
       * @f[
       *   \int_\Omega A(u) : B(v) \ dx
       * @f]
       *
       * @param[in] prod Dot product instance
       */
      Integral(const IntegrandType& prod)
        : Parent(prod)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Integral to copy
       */
      Integral(const Integral& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Integral to move from
       */
      Integral(Integral&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Gets the integration region.
       * @return Region::Cells indicating integration over mesh elements
       *
       * Domain integrals are computed over all cells @f$ K \in \mathcal{T}_h @f$
       * in the mesh.
       */
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Cells;
      }

      /**
       * @brief Creates a polymorphic copy of this integral.
       * @return Pointer to a new copy
       */
      Integral* copy() const noexcept override
      {
        return new Integral(*this);
      }
  };

  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  Integral(const Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>&)
    -> Integral<Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>;

  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  Integral(const ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>&, const ShapeFunctionBase<RHSDerived, TestFES, TestSpace>&)
    -> Integral<Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>;

  /**
   * @ingroup IntegralSpecializations
   * @brief Integration of a test operator.
   *
   * Given an operator defined over a test space @f$ V_h @f$
   * @f[
   *   A : V_h \rightarrow \mathbb{R},
   * @f]
   * this class will represent its integral
   * @f[
   *   \int_{\mathcal{T}_h} A(v) \ dx \ .
   * @f]
   */
  template <class NestedDerived, class FES>
  class Integral<ShapeFunctionBase<NestedDerived, FES, TestSpace>> final
    : public QuadratureRule<ShapeFunctionBase<NestedDerived, FES, TestSpace>>
  {
    public:
      using IntegrandType = ShapeFunctionBase<NestedDerived, FES, TestSpace>;

      using Parent = QuadratureRule<IntegrandType>;

      template <class LHSDerived, class RHSDerived>
      Integral(const FunctionBase<LHSDerived>& lhs, const ShapeFunctionBase<RHSDerived, FES, TestSpace>& rhs)
        : Integral(Dot(lhs, rhs))
      {}

      Integral(const IntegrandType& integrand)
        : Parent(integrand)
      {}

      Integral(const Integral& other)
        : Parent(other)
      {}

      Integral(Integral&& other)
        : Parent(std::move(other))
      {}

      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Cells;
      }

      Integral* copy() const noexcept override
      {
        return new Integral(*this);
      }
  };

  template <class NestedDerived, class FES>
  Integral(const ShapeFunctionBase<NestedDerived, FES, TestSpace>&)
    -> Integral<ShapeFunctionBase<NestedDerived, FES, TestSpace>>;

  template <class LHSDerived, class RHSDerived, class FES>
  Integral(const FunctionBase<LHSDerived>&, const ShapeFunctionBase<RHSDerived, FES, TestSpace>&)
    -> Integral<ShapeFunctionBase<Dot<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, FES, TestSpace>>>>;

  /**
   * @ingroup IntegralSpecializations
   * @brief Integration of a GridFunction object.
   */
  template <class FES, class Data>
  class Integral<GridFunction<FES, Data>> final
    : public QuadratureRule<GridFunction<FES, Data>>
  {
    public:
      /// Type of integrand
      using IntegrandType = GridFunction<FES, Data>;

      /// Parent class
      using Parent = QuadratureRule<IntegrandType>;

      /**
       * @brief Constructs the integral object
       */
      Integral(const IntegrandType& u)
        : Parent(u)
      {
        assert(u.getFiniteElementSpace().getVectorDimension() == 1);
      }

      Integral(const Integral& other)
        : Parent(other)
      {}

      Integral(Integral&& other)
        : Parent(std::move(other))
      {}

      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Cells;
      }

      Integral* copy() const noexcept override
      {
        return new Integral(*this);
      }
  };

  template <class FES, class Data>
  Integral(const GridFunction<FES, Data>&) -> Integral<GridFunction<FES, Data>>;
}

#endif
