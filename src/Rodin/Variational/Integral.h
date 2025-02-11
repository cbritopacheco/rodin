/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_INTEGRAL_H
#define RODIN_VARIATIONAL_INTEGRAL_H

#include <cassert>
#include <set>
#include <utility>

#include "Rodin/FormLanguage/Base.h"

#include "Dot.h"
#include "Function.h"
#include "LinearForm.h"
#include "ForwardDecls.h"
#include "GridFunction.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "FiniteElement.h"
#include "MatrixFunction.h"
#include "QuadratureRule.h"
#include "LinearFormIntegrator.h"
#include "BilinearFormIntegrator.h"

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

      Integral(const Integral& other)
        : Parent(other)
      {}

      Integral(Integral&& other)
        : Parent(std::move(other))
      {}

      Integrator::Region getRegion() const override
      {
        return Integrator::Region::Cells;
      }

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

      Integrator::Region getRegion() const override
      {
        return Integrator::Region::Cells;
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
  template <class FES>
  class Integral<GridFunction<FES>> final
    : public QuadratureRule<GridFunction<FES>>
  {
    public:
      /// Type of integrand
      using IntegrandType = GridFunction<FES>;

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

      Integrator::Region getRegion() const override
      {
        return Integrator::Region::Cells;
      }

      Integral* copy() const noexcept override
      {
        return new Integral(*this);
      }
  };

  template <class FES>
  Integral(const GridFunction<FES>&) -> Integral<GridFunction<FES>>;
}

#endif
