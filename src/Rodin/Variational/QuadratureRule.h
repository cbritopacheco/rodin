/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file QuadratureRule.h
 * @brief Quadrature rule classes for numerical integration.
 *
 * This file defines the QuadratureRule classes which implement numerical
 * quadrature (integration) rules for evaluating integrals over mesh elements.
 * Quadrature rules are fundamental to finite element assembly.
 *
 * ## Mathematical Foundation
 * A quadrature rule approximates an integral as a weighted sum:
 * @f[
 *   \int_K f(x) \, dx \approx \sum_{i=1}^{n_q} w_i f(x_i)
 * @f]
 * where:
 * - @f$ x_i @f$ are quadrature points
 * - @f$ w_i @f$ are quadrature weights
 * - @f$ n_q @f$ is the number of quadrature points
 *
 * ## Accuracy
 * A quadrature rule is exact for polynomials up to degree @f$ p @f$ if:
 * @f[
 *   \int_K q(x) \, dx = \sum_{i=1}^{n_q} w_i q(x_i) \quad \forall q \in \mathbb{P}_p
 * @f]
 *
 * ## Common Quadrature Rules
 * - **Gauss-Legendre**: Optimal for polynomial integrands
 * - **Gauss-Lobatto**: Includes element vertices
 * - **Vertex quadrature**: Simple but low accuracy
 * - **Midpoint rule**: Single point at element center
 *
 * ## Applications in FEM
 * - Evaluating stiffness matrix entries: @f$ A_{ij} = \int_K \nabla \phi_i \cdot \nabla \phi_j \, dx @f$
 * - Evaluating load vector entries: @f$ b_i = \int_K f \phi_i \, dx @f$
 * - Computing functionals and error norms
 *
 * ## Usage Example
 * ```cpp
 * // Quadrature rule automatically selected based on element type
 * auto qr = QuadratureRule(integrand);
 * qr.setPolytope(element);
 * Real integral_value = qr.compute();
 * ```
 */
#ifndef RODIN_VARIATIONAL_QUADRATURERULE_H
#define RODIN_VARIATIONAL_QUADRATURERULE_H

#include "ForwardDecls.h"

#include "Rodin/QF/GenericPolytopeQuadrature.h"

#include "ShapeFunction.h"
#include "LinearFormIntegrator.h"
#include "BilinearFormIntegrator.h"

namespace Rodin::Variational
{
  /**
   * @defgroup QuadratureRuleSpecializations QuadratureRule Template Specializations
   * @brief Template specializations of the QuadratureRule class.
   *
   * @see QuadratureRule
   * @see RodinQuadrature
   */

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Quadrature rule for integrating functions on mesh polytopes.
   *
   * QuadratureRule manages the numerical integration of functions over mesh
   * elements using appropriate quadrature formulas for each polytope type.
   */
  template <class FunctionDerived>
  class QuadratureRule<FunctionBase<FunctionDerived>> final
    : public FormLanguage::Base
  {
    public:
      using IntegrandType = FunctionBase<FunctionDerived>;

      using IntegrandRangeType = typename FormLanguage::Traits<IntegrandType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<IntegrandRangeType>::ScalarType;

      using Parent = FormLanguage::Base;

      QuadratureRule(const IntegrandType& f)
        : m_integrand(f.copy()),
          m_polytope(nullptr)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_polytope(nullptr),
          m_integrand(other.m_integrand->copy()),
          m_qfgg(other.m_qfgg),
          m_qf(other.m_qf)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_integrand(std::move(other.m_integrand)),
          m_qfgg(std::move(other.m_qfgg)),
          m_qf(std::move(other.m_qf))
      {}

      const Geometry::Polytope& getPolytope() const
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope)
      {
        m_polytope = &polytope;
        if (!m_qf)
        {
          m_qfgg.emplace(polytope.getGeometry());
          m_qf = m_qfgg.value();
        }
        const auto& qf = m_qf.value().get();
        m_ps.clear();
        m_ps.reserve(qf.getSize());
        for (size_t i = 0; i < qf.getSize(); i++)
          m_ps.emplace_back(polytope, qf.getPoint(i));
        return *this;
      }

      ScalarType compute()
      {
        auto& res = m_value.emplace(0);
        const auto& qf = getQuadratureFormula();
        const auto& f = getIntegrand();
        assert(m_ps.size() == qf.getSize());
        for (size_t i = 0; i < m_ps.size(); i++)
          res += qf.getWeight(i) * m_ps[i].getDistortion() * f(m_ps[i]);
        return res;
      }

      const Optional<ScalarType>& getValue() const
      {
        return m_value;
      }

      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const QF::QuadratureFormulaBase& getQuadratureFormula() const
      {
        assert(m_qf);
        return m_qf.value().get();
      }

      QuadratureRule& setQuadratureFormula(const QF::QuadratureFormulaBase& qf)
      {
        m_qf = qf;
        return *this;
      }

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const Geometry::Polytope* m_polytope;
      Optional<const QF::GenericPolytopeQuadrature> m_qfgg;
      Optional<std::reference_wrapper<const QF::QuadratureFormulaBase>> m_qf;
      Optional<ScalarType> m_value;

      std::vector<Geometry::Point> m_ps;
  };

  /**
   * @ingroup IntegralSpecializations
   * @brief Integration of a GridFunction object.
   */
  template <class FES, class Data>
  class QuadratureRule<GridFunction<FES, Data>> : public Integrator
  {
    public:
      using FESType = FES;

      using IntegrandType = GridFunction<FESType, Data>;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using Parent = Integrator;

      /**
       * @brief Constructs the integral object from the given integrand.
       */
      QuadratureRule(const IntegrandType& u)
        : m_u(u),
          m_v(u.getFiniteElementSpace()),
          m_lf(m_v)
      {
        assert(u.getFiniteElementSpace().getVectorDimension() == 1);
      }

      /**
       * @brief Copy constructor.
       */
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_u(other.m_u),
          m_v(other.m_u.get().getFiniteElementSpace()),
          m_lf(m_v)
      {}

      /**
       * @brief Move constructor.
       */
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_v(std::move(other.m_v)),
          m_lf(std::move(other.m_lf))
      {}

      /**
       * @brief Integrates the expression and returns the value.
       *
       * Compute the value of the integral, caches it and returns it.
       *
       * @returns Value of integral
       */
      ScalarType compute()
      {
        switch (getRegion())
        {
          case Geometry::Region::Cells:
          {
            auto lfi = Variational::Integral(m_v);
            if (m_attrs.size() > 0)
              lfi.over(m_attrs);
            m_lf = lfi;
            m_lf.assemble();
            return m_value.emplace(m_lf(m_u.get()));
          }
          case Geometry::Region::Boundary:
          {
            auto lfi = Variational::BoundaryIntegral(m_v);
            if (m_attrs.size() > 0)
              lfi.over(m_attrs);
            m_lf = lfi;
            m_lf.assemble();
            return m_value.emplace(m_lf(m_u.get()));
          }
          case Geometry::Region::Faces:
          {
            auto lfi = Variational::FaceIntegral(m_v);
            if (m_attrs.size() > 0)
              lfi.over(m_attrs);
            m_lf = lfi;
            m_lf.assemble();
            return m_value.emplace(m_lf(m_u.get()));
          }
          case Geometry::Region::Interface:
          {
            auto lfi = Variational::InterfaceIntegral(m_v);
            if (m_attrs.size() > 0)
              lfi.over(m_attrs);
            m_lf = lfi;
            m_lf.assemble();
            return m_value.emplace(m_lf(m_u.get()));
          }
        }
        assert(false);
        return Math::nan<ScalarType>();
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
      operator ScalarType()
      {
        if (!m_value.has_value())
          return compute();
        else
          return m_value.value();
      }

      QuadratureRule& over(Geometry::Attribute attr)
      {
        return over(FlatSet<Geometry::Attribute>{attr});
      }

      QuadratureRule& over(const FlatSet<Geometry::Attribute>& attrs)
      {
        m_attrs = attrs;
        return *this;
      }

      const Optional<ScalarType>& getValue() const
      {
        return m_value;
      }

      Type getType() const override
      {
        return Integrator::Type::Linear;
      }

      virtual Geometry::Region getRegion() const = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::reference_wrapper<const GridFunction<FES, Data>>   m_u;
      TestFunction<FES>                                 m_v;

      FlatSet<Geometry::Attribute>                  m_attrs;
      LinearForm<FES, Math::Vector<ScalarType>>     m_lf;

      Optional<ScalarType> m_value;
  };

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Approximation of the integral of the the dot product between a
   * trial shape function and a test shape function.
   *
   * Represents the quadrature rule approximation of an integral:
   * @f[
   *  \int_{\mathcal{R}_h} \mathrm{IntegrandType} \ dx \approx \sum_{i = 1}^{n}
   *  w_i \ \mathrm{IntegrandType} (x_i)
   * @f]
   * where @f$ \mathcal{R}_h @f$ is some region of the mesh @f$ \mathcal{T}_h
   * @f$, the quadrature point @f$ x_i @f$ has an associated weight @f$ w_i @f$
   * and @f$ \mathrm{IntegrandType}(x_i) @f$ is the value of the integrand at the
   * quadrature point.
   */
  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
      ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<
            ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
            ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>::ScalarType>
  {
    public:
      using LHSType = ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>;

      using RHSType = ShapeFunctionBase<RHSDerived, TestFES, TestSpace>;

      using IntegrandType = Dot<LHSType, RHSType>;

      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      QuadratureRule(const LHSType& lhs, const RHSType& rhs)
        : QuadratureRule(Dot(lhs, rhs))
      {}

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_polytope(nullptr),
          m_set(false)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::move(other.m_qf)),
          m_ps(std::move(other.m_ps)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::move(other.m_set)),
          m_order(std::move(other.m_order)),
          m_geometry(std::move(other.m_geometry))
      {}

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) override
      {
        m_polytope = &polytope;
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& integrand = *m_integrand;
        const auto& trial = integrand.getLHS();
        const auto& test = integrand.getRHS();
        const auto& trialfes = trial.getFiniteElementSpace();
        const auto& testfes = test.getFiniteElementSpace();
        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe = testfes.getFiniteElement(d, idx);
        const size_t order = std::max(trialfe.getOrder(), testfe.getOrder());
        const auto& geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;
        if (recompute)
        {
          m_set = true;
          m_order = order;
          m_geometry = geometry;
          m_qf.reset(new QF::GenericPolytopeQuadrature(order, geometry));
          m_ps.clear();
          m_ps.reserve(m_qf->getSize());
          for (size_t i = 0; i < m_qf->getSize(); i++)
            m_ps.emplace_back(polytope, m_qf->getPoint(i));
        }
        else
        {
          for (size_t i = 0; i < m_qf->getSize(); i++)
            m_ps[i].setPolytope(polytope);
        }
        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        ScalarType res = 0;
        auto& integrand = *m_integrand;
        for (size_t i = 0; i < m_ps.size(); i++)
        {
          integrand.setPoint(m_ps[i]);
          res += m_qf->getWeight(i) * m_ps[i].getDistortion() * integrand(tr, te);
        }
        return res;
      }

      virtual Geometry::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      std::unique_ptr<QF::QuadratureFormulaBase> m_qf;
      std::vector<Geometry::Point> m_ps;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;
  };

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Approximation of the integral of a test shape function.
   */
  template <class NestedDerived, class FES>
  class QuadratureRule<ShapeFunctionBase<NestedDerived, FES, TestSpace>>
    : public LinearFormIntegratorBase<
        typename FormLanguage::Traits<ShapeFunctionBase<NestedDerived, FES, TestSpace>>::ScalarType>
  {
    public:
      using FESType = FES;

      using IntegrandType = ShapeFunctionBase<NestedDerived, FESType, TestSpace>;

      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent = LinearFormIntegratorBase<ScalarType>;

      template <class LHSDerived, class RHSDerived>
      constexpr
      QuadratureRule(const FunctionBase<LHSDerived>& lhs, const ShapeFunctionBase<RHSDerived, FES, TestSpace>& rhs)
        : QuadratureRule(Dot(lhs, rhs))
      {}

      constexpr
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLeaf()),
          m_integrand(integrand.copy()),
          m_polytope(nullptr),
          m_set(false)
      {}

      constexpr
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_polytope(other.m_polytope),
          m_set(false)
      {}

      constexpr
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::move(other.m_qf)),
          m_ps(std::move(other.m_ps)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::move(other.m_set)),
          m_order(std::move(other.m_order)),
          m_geometry(std::move(other.m_geometry))
      {}

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& integrand = getIntegrand();
        const auto& fes = integrand.getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, idx);
        const size_t order = fe.getOrder();
        const auto& geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;
        if (recompute)
        {
          m_set = true;
          m_order = order;
          m_geometry = geometry;
          m_qf.reset(new QF::GenericPolytopeQuadrature(order, geometry));
          m_ps.clear();
          m_ps.reserve(m_qf->getSize());
          for (size_t i = 0; i < m_qf->getSize(); i++)
            m_ps.emplace_back(polytope, m_qf->getPoint(i));
        }
        else
        {
          for (size_t i = 0; i < m_qf->getSize(); i++)
            m_ps[i].setPolytope(polytope);
        }
        return *this;
      }

      ScalarType integrate(size_t local) final override
      {
        ScalarType res = 0;
        auto& integrand = *m_integrand;
        for (size_t i = 0; i < m_ps.size(); i++)
        {
          integrand.setPoint(m_ps[i]);
          res += m_qf->getWeight(i) * m_ps[i].getDistortion() * integrand.getBasis(local);
        }
        return res;
      }

      virtual Geometry::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      std::unique_ptr<QF::QuadratureFormulaBase> m_qf;
      std::vector<Geometry::Point> m_ps;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;
  };
}

#endif
