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
   *   \int_{\mathcal{T}_h} A(u) : B(v) \ dx
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
      using LHS = ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>;

      /// Type of the right operand of the dot product
      using RHS = ShapeFunctionBase<RHSDerived, TestFES, TestSpace>;

      /// Type of the integrand
      using Integrand = Dot<LHS, RHS>;

      /// Parent class
      using Parent = QuadratureRule<Integrand>;

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
      Integral(const LHS& lhs, const RHS& rhs)
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
      Integral(const Integrand& prod)
        : Parent(prod)
      {}

      Integral(const Integral& other)
        : Parent(other)
      {}

      Integral(Integral&& other)
        : Parent(std::move(other))
      {}

      inline Integrator::Region getRegion() const override
      {
        return Integrator::Region::Domain;
      }

      inline Integral* copy() const noexcept override
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
   *   \int_\Omega A(v) \ dx \ .
   * @f]
   */
  template <class NestedDerived, class FES>
  class Integral<ShapeFunctionBase<NestedDerived, FES, TestSpace>> final
    : public QuadratureRule<ShapeFunctionBase<NestedDerived, FES, TestSpace>>
  {
    public:
      using Integrand = ShapeFunctionBase<NestedDerived, FES, TestSpace>;
      using Parent = QuadratureRule<Integrand>;

      template <class LHSDerived, class RHSDerived>
      Integral(const FunctionBase<LHSDerived>& lhs, const ShapeFunctionBase<RHSDerived, FES, TestSpace>& rhs)
        : Integral(Dot(lhs, rhs))
      {}

      Integral(const Integrand& integrand)
        : Parent(integrand)
      {}

      Integral(const Integral& other)
        : Parent(other)
      {}

      Integral(Integral&& other)
        : Parent(std::move(other))
      {}

      inline Integrator::Region getRegion() const override
      {
        return Integrator::Region::Domain;
      }

      inline Integral* copy() const noexcept override
      {
        return new Integral(*this);
      }
  };

  template <class NestedDerived, class FES>
  Integral(const ShapeFunctionBase<NestedDerived, FES, TestSpace>&)
    -> Integral<ShapeFunctionBase<NestedDerived, FES, TestSpace>>;

  template <class LHSDerived, class RHSDerived, class FES>
  Integral(const FunctionBase<LHSDerived>&, const ShapeFunctionBase<RHSDerived, FES, TestSpace>&)
    -> Integral<ShapeFunctionBase<Dot<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, FES, TestSpace>>, FES, TestSpace>>;

  /**
   * @ingroup IntegralSpecializations
   * @brief Integration of a GridFunction object.
   */
  template <class FES>
  class Integral<GridFunction<FES>> final : public FormLanguage::Base
  {
    public:
      /// Type of integrand
      using Integrand = GridFunction<FES>;

      /// Parent class
      using Parent = FormLanguage::Base;

      /**
       * @brief Constructs the integral object
       */
      Integral(const Integrand& u)
        : m_u(u),
          m_v(u.getFiniteElementSpace()),
          m_lf(m_v),
          m_assembled(false)
      {
        assert(u.getFiniteElementSpace().getVectorDimension() == 1);
        m_lf = Variational::Integral(m_v); // Prefix with namespace so CTAD kicks in.
      }

      Integral(const Integral& other)
        : Parent(other),
          m_u(other.m_u),
          m_v(other.m_u.get().getFiniteElementSpace()),
          m_lf(m_v),
          m_assembled(false)
      {}

      Integral(Integral&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_v(std::move(other.m_v)),
          m_lf(std::move(other.m_lf)),
          m_assembled(std::move(other.m_assembled))
      {}

      /**
       * @brief Integrates the expression and returns the value.
       *
       * Compute the value of the integral, caches it and returns it.
       *
       * @returns Value of integral
       */
      inline
      Scalar compute()
      {
        if (!m_assembled)
          m_lf.assemble();
        m_assembled = true;
        return m_value.emplace(m_lf(m_u));
      }

      /**
       * @brief Returns the value of the integral, computing it if necessary.
       *
       * If compute() has been called before, returns the value of the cached
       * value. Otherwise, it will call compute() and return the newly computed
       * value.
       *
       * @returns Value of integral
       */
      inline
      operator Scalar()
      {
        if (!m_value.has_value())
          return compute();
        else
          return m_value.value();
      }

      inline
      const std::optional<Scalar>& getValue() const
      {
        return m_value;
      }

      inline Integral* copy() const noexcept override
      {
        return new Integral(*this);
      }

    private:
      std::reference_wrapper<const GridFunction<FES>>   m_u;
      TestFunction<FES>                                 m_v;

      LinearForm<FES, Context::Serial, Math::Vector>    m_lf;
      bool m_assembled;

      std::optional<Scalar> m_value;
  };

  template <class FES>
  Integral(const GridFunction<FES>&) -> Integral<GridFunction<FES>>;

}

#endif
