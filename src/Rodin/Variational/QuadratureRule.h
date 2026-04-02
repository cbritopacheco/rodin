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
 * This file defines specializations of @ref Rodin::Variational::QuadratureRule
 * for integrating functions and shape-function expressions over mesh
 * polytopes.
 *
 * ## Mathematical foundation
 *
 * A quadrature rule approximates an integral by a weighted sum:
 * @f[
 *   \int_K f(x)\,dx \approx \sum_{q=1}^{n_q} w_q\,f(x_q),
 * @f]
 * where:
 * - @f$ x_q @f$ are quadrature points,
 * - @f$ w_q @f$ are quadrature weights,
 * - @f$ n_q @f$ is the number of quadrature points.
 *
 * In Rodin, the reference quadrature formula lives in the @ref Rodin::QF
 * module, while the mapped geometric points live in
 * @ref Rodin::Geometry::PolytopeQuadrature. This separation allows:
 * - canonical reuse of reference quadrature formulas,
 * - mesh-owned caching of mapped geometric quadrature points,
 * - simpler and more efficient integrator implementations.
 *
 * ## Usage
 *
 * Typical usage is:
 * @code{.cpp}
 * auto qr = QuadratureRule(integrand);
 * qr.setPolytope(cell);
 * const auto value = qr.compute();
 * @endcode
 *
 * For local element integrators, @ref setPolytope binds the integrator to a
 * concrete polytope, selects an appropriate quadrature formula, and retrieves
 * the corresponding cached @ref Rodin::Geometry::PolytopeQuadrature from the
 * mesh.
 */
#ifndef RODIN_VARIATIONAL_QUADRATURERULE_H
#define RODIN_VARIATIONAL_QUADRATURERULE_H

#include "ForwardDecls.h"

#include "Rodin/Geometry/PolytopeQuadrature.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"

#include "IntegrationPoint.h"
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
   * @brief Quadrature rule for integrating general functions on mesh polytopes.
   *
   * This specialization evaluates scalar-valued function integrals on a single
   * polytope using a quadrature formula and the corresponding mapped geometric
   * quadrature points cached on the mesh.
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

      /**
       * @brief Constructs the quadrature rule for the given integrand.
       * @param[in] f Integrand to be evaluated
       */
      QuadratureRule(const IntegrandType& f)
        : m_integrand(f.copy()),
          m_polytope(nullptr),
          m_qf(nullptr),
          m_quadrature(nullptr)
      {}

      /**
       * @brief Copy constructor.
       *
       * The cached value and polytope binding are not copied. The quadrature
       * formula pointer is copied because it refers to a canonical formula.
       */
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_polytope(nullptr),
          m_qf(other.m_qf),
          m_quadrature(nullptr)
      {}

      /**
       * @brief Move constructor.
       */
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::exchange(other.m_quadrature, nullptr)),
          m_value(std::move(other.m_value))
      {}

      /**
       * @brief Gets the currently bound polytope.
       * @returns Bound polytope
       */
      const Geometry::Polytope& getPolytope() const
      {
        assert(m_polytope);
        return *m_polytope;
      }

      /**
       * @brief Binds the quadrature rule to a polytope.
       * @param[in] polytope Polytope on which the integral is to be evaluated
       * @returns Reference to this object
       *
       * If no explicit quadrature formula was set previously, a default generic
       * polytope quadrature of order 1 is selected for the polytope geometry.
       * The corresponding mapped quadrature is then retrieved from the mesh.
       */
      QuadratureRule& setPolytope(const Geometry::Polytope& polytope)
      {
        m_polytope = &polytope;
        m_value.reset();

        if (!m_qf)
          m_qf = &QF::PolytopeQuadratureFormula::get(1, polytope.getGeometry());

        assert(m_qf);
        m_quadrature = &polytope.getQuadrature(*m_qf);
        return *this;
      }

      /**
       * @brief Computes the integral value on the currently bound polytope.
       * @returns Approximated integral value
       */
      ScalarType compute()
      {
        auto& res = m_value.emplace(0);
        const auto& qf = getQuadratureFormula();
        const auto& f = getIntegrand();

        assert(m_quadrature);
        const auto& q = *m_quadrature;

        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          res += qf.getWeight(qp) * p.getDistortion() * f(p);
        }

        return res;
      }

      /**
       * @brief Gets the cached value, if available.
       * @returns Cached integral value
       */
      const Optional<ScalarType>& getValue() const
      {
        return m_value;
      }

      /**
       * @brief Gets the integrand.
       * @returns Integrand reference
       */
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      /**
       * @brief Gets the quadrature formula currently used.
       * @returns Quadrature formula reference
       */
      const QF::QuadratureFormulaBase& getQuadratureFormula() const
      {
        assert(m_qf);
        return *m_qf;
      }

    private:
      std::unique_ptr<IntegrandType> m_integrand;                  ///< Integrand
      const Geometry::Polytope* m_polytope;                        ///< Bound polytope
      const QF::QuadratureFormulaBase* m_qf;                       ///< Reference quadrature formula
      const Geometry::PolytopeQuadrature* m_quadrature;            ///< Mapped geometric quadrature
      Optional<ScalarType> m_value;                                ///< Cached value
  };

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Integration of a scalar grid function over a mesh region.
   *
   * This specialization delegates the computation to linear forms over the
   * requested region. It is not a per-polytope quadrature rule and therefore
   * does not use @ref Geometry::PolytopeQuadrature directly.
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
       * @brief Constructs the integral object from the given grid function.
       * @param[in] u Grid function to integrate
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
       * @brief Computes the integral value over the selected region.
       * @returns Integral value
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
       * @brief Returns the integral value, computing it if necessary.
       * @returns Integral value
       */
      operator ScalarType()
      {
        if (!m_value.has_value())
          return compute();
        else
          return m_value.value();
      }

      /**
       * @brief Restricts integration to one attribute.
       * @param[in] attr Attribute
       * @returns Reference to this object
       */
      QuadratureRule& over(Geometry::Attribute attr)
      {
        return over(FlatSet<Geometry::Attribute>{ attr });
      }

      /**
       * @brief Restricts integration to a set of attributes.
       * @param[in] attrs Attributes
       * @returns Reference to this object
       */
      QuadratureRule& over(const FlatSet<Geometry::Attribute>& attrs)
      {
        m_attrs = attrs;
        return *this;
      }

      /**
       * @brief Gets the cached value, if available.
       * @returns Cached value
       */
      const Optional<ScalarType>& getValue() const
      {
        return m_value;
      }

      /**
       * @brief Gets the integrator type.
       * @returns Linear integrator type
       */
      Type getType() const override
      {
        return Integrator::Type::Linear;
      }

      /**
       * @brief Gets the region over which the integral is defined.
       * @returns Integration region
       */
      virtual Geometry::Region getRegion() const = 0;

      /**
       * @brief Polymorphic copy.
       * @returns Heap-allocated copy
       */
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::reference_wrapper<const GridFunction<FES, Data>> m_u;   ///< Integrated grid function
      TestFunction<FES> m_v;                                       ///< Auxiliary test function
      FlatSet<Geometry::Attribute> m_attrs;                        ///< Region attribute filter
      LinearForm<FES, Data> m_lf;                                  ///< Auxiliary linear form
      Optional<ScalarType> m_value;                                ///< Cached value
  };

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Local bilinear-form quadrature for the dot product of a trial and a
   * test shape function.
   *
   * This specialization assembles the local matrix associated with:
   * @f[
   *   \int_K \phi^{\mathrm{tr}} \cdot \phi^{\mathrm{te}} \, dx
   * @f]
   * or a more general integrand of the same expression-template type.
   *
   * The quadrature formula is chosen from the integrand order if available,
   * otherwise from the finite element orders.
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

      /**
       * @brief Constructs from explicit left and right shape functions.
       * @param[in] lhs Trial-side shape function
       * @param[in] rhs Test-side shape function
       */
      QuadratureRule(const LHSType& lhs, const RHSType& rhs)
        : QuadratureRule(Dot(lhs, rhs))
      {}

      /**
       * @brief Constructs from the full integrand expression.
       * @param[in] integrand Expression to integrate
       */
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Copy constructor.
       *
       * The bound polytope and mapped quadrature are not copied.
       */
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(other.m_qf),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Move constructor.
       */
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::exchange(other.m_quadrature, nullptr)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_mat(std::move(other.m_mat))
      {}

      /**
       * @brief Gets the integrand.
       * @returns Integrand reference
       */
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      /**
       * @brief Gets the currently bound polytope.
       * @returns Bound polytope
       */
      const Geometry::Polytope& getPolytope() const override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      /**
       * @brief Binds the integrator to a concrete polytope and assembles the
       * local matrix.
       * @param[in] polytope Target polytope
       * @returns Reference to this object
       */
      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) override
      {
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        assert(m_integrand);
        auto& integrand = *m_integrand;

        const auto& trial = integrand.getLHS();
        const auto& test  = integrand.getRHS();

        const auto& trialfes = trial.getFiniteElementSpace();
        const auto& testfes  = test.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes.getFiniteElement(d, idx);

        const auto geometry = polytope.getGeometry();

        const size_t order =
          integrand.getOrder(polytope).value_or(trialfe.getOrder() + testfe.getOrder());

        const bool recompute =
          !m_set || (m_order != order) || (m_geometry != geometry);

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;
          m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);
        }

        assert(m_qf);
        m_quadrature = &polytope.getQuadrature(*m_qf);

        const size_t ntr = trial.getDOFs(*m_polytope);
        const size_t nte = test.getDOFs(*m_polytope);

        m_mat.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_mat.setZero();

        // Eigen is assumed ColMajor. Columns are filled contiguously.
        ScalarType* __restrict M = m_mat.data();
        const Eigen::Index ld = m_mat.outerStride();

        assert(m_quadrature);
        const auto& q = *m_quadrature;

        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);

          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp))
            * static_cast<ScalarType>(p.getDistortion());

          const IntegrationPoint ip(p, *m_qf, qp);
          integrand.setIntegrationPoint(ip);

          for (size_t tr = 0; tr < ntr; ++tr)
          {
            ScalarType* __restrict col = M + static_cast<Eigen::Index>(tr) * ld;
            const auto& phi_tr = trial.getBasis(tr);

            for (size_t te = 0; te < nte; ++te)
            {
              const auto& phi_te = test.getBasis(te);
              col[static_cast<Eigen::Index>(te)] += wdet * Math::dot(phi_tr, phi_te);
            }
          }
        }

        return *this;
      }

      /**
       * @brief Gets one entry of the assembled local matrix.
       * @param[in] tr Trial local index
       * @param[in] te Test local index
       * @returns Local matrix coefficient
       */
      inline ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_mat(te, tr); // rows=test, cols=trial
      }

      /**
       * @brief Gets the geometric region supported by this integrator.
       * @returns Region type
       */
      virtual Geometry::Region getRegion() const override = 0;

      /**
       * @brief Polymorphic copy.
       * @returns Heap-allocated copy
       */
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;               ///< Integrand expression
      const QF::QuadratureFormulaBase* m_qf;                    ///< Reference quadrature formula
      const Geometry::PolytopeQuadrature* m_quadrature;         ///< Mapped geometric quadrature
      const Geometry::Polytope* m_polytope;                     ///< Bound polytope
      bool m_set;                                               ///< Whether formula selection data are initialized
      size_t m_order;                                           ///< Cached quadrature order
      Geometry::Polytope::Type m_geometry;                      ///< Cached geometry type
      Math::Matrix<ScalarType> m_mat;                           ///< Local matrix, rows=test, cols=trial
  };

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Local linear-form quadrature for a test shape function expression.
   *
   * This specialization assembles the local vector associated with:
   * @f[
   *   \int_K \phi^{\mathrm{te}} \, dx
   * @f]
   * or a more general integrand of the same expression-template type.
   */
  template <class NestedDerived, class FES>
  class QuadratureRule<ShapeFunctionBase<NestedDerived, FES, TestSpace>>
    : public LinearFormIntegratorBase<
        typename FormLanguage::Traits<
          ShapeFunctionBase<NestedDerived, FES, TestSpace>>::ScalarType>
  {
    public:
      using FESType = FES;
      using IntegrandType = ShapeFunctionBase<NestedDerived, FESType, TestSpace>;
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;
      using Parent = LinearFormIntegratorBase<ScalarType>;

      /**
       * @brief Constructs from a function-times-test-function expression.
       * @param[in] lhs Coefficient function
       * @param[in] rhs Test shape function
       */
      template <class LHSDerived, class RHSDerived>
      constexpr QuadratureRule(
          const FunctionBase<LHSDerived>& lhs,
          const ShapeFunctionBase<RHSDerived, FES, TestSpace>& rhs)
        : QuadratureRule(Dot(lhs, rhs))
      {}

      /**
       * @brief Constructs from the full integrand expression.
       * @param[in] integrand Expression to integrate
       */
      constexpr QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Copy constructor.
       *
       * The bound polytope and mapped quadrature are not copied.
       */
      constexpr QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(other.m_qf),
          m_quadrature(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Move constructor.
       */
      constexpr QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::exchange(other.m_quadrature, nullptr)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_vec(std::move(other.m_vec))
      {}

      /**
       * @brief Gets the integrand.
       * @returns Integrand reference
       */
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      /**
       * @brief Gets the currently bound polytope.
       * @returns Bound polytope
       */
      const Geometry::Polytope& getPolytope() const override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      /**
       * @brief Binds the integrator to a concrete polytope and assembles the
       * local vector.
       * @param[in] polytope Target polytope
       * @returns Reference to this object
       */
      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        assert(m_integrand);
        auto& integrand = *m_integrand;

        const auto& fes = integrand.getFiniteElementSpace();
        const auto& fe  = fes.getFiniteElement(d, idx);

        const auto geometry = polytope.getGeometry();

        const size_t order =
          integrand.getOrder(polytope).value_or(fe.getOrder());

        const bool recompute =
          !m_set || (m_order != order) || (m_geometry != geometry);

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;
          m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);
        }

        assert(m_qf);
        m_quadrature = &polytope.getQuadrature(*m_qf);

        const size_t nte = integrand.getDOFs(polytope);

        m_vec.resize(static_cast<Eigen::Index>(nte));
        m_vec.setZero();

        ScalarType* __restrict v = m_vec.data();

        assert(m_quadrature);
        const auto& q = *m_quadrature;

        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);

          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp))
            * static_cast<ScalarType>(p.getDistortion());

          const IntegrationPoint ip(p, *m_qf, qp);
          integrand.setIntegrationPoint(ip);

          for (size_t te = 0; te < nte; ++te)
            v[static_cast<Eigen::Index>(te)] += wdet * integrand.getBasis(te);
        }

        return *this;
      }

      /**
       * @brief Gets one entry of the assembled local vector.
       * @param[in] local Local test index
       * @returns Local vector coefficient
       */
      inline ScalarType integrate(size_t local) final override
      {
        return m_vec(local);
      }

      /**
       * @brief Gets the geometric region supported by this integrator.
       * @returns Region type
       */
      virtual Geometry::Region getRegion() const override = 0;

      /**
       * @brief Polymorphic copy.
       * @returns Heap-allocated copy
       */
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;               ///< Integrand expression
      const QF::QuadratureFormulaBase* m_qf;                    ///< Reference quadrature formula
      const Geometry::PolytopeQuadrature* m_quadrature;         ///< Mapped geometric quadrature
      const Geometry::Polytope* m_polytope;                     ///< Bound polytope
      bool m_set;                                               ///< Whether formula selection data are initialized
      size_t m_order;                                           ///< Cached quadrature order
      Geometry::Polytope::Type m_geometry;                      ///< Cached geometry type
      Math::Vector<ScalarType> m_vec;                           ///< Local vector
  };
}

#endif
