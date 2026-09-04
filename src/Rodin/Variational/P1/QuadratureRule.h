/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file QuadratureRule.h
 * @brief Optimized quadrature rules for P1 finite element spaces.
 *
 * This file provides specialized quadrature rules for P1 spaces that exploit
 * the constant gradient property of P1 basis functions to achieve exact
 * integration with minimal quadrature points.
 *
 * ## Quadrature Strategies
 * - **Centroid quadrature**: Uses single point at element barycenter
 * - **Exact for P1 bilinear forms**: @f$ \int \nabla \phi_i \cdot \nabla \phi_j @f$
 * - **Reduced integration**: Enables efficient assembly
 *
 * ## Supported Integrands
 * - Scalar P1 shape functions: @f$ \int v \, dx @f$
 * - Dot products: @f$ \int f \cdot v \, dx @f$
 * - Mass forms: @f$ \int (Au) \cdot v \, dx @f$
 * - Stiffness forms: @f$ \int \nabla u \cdot \nabla v \, dx @f$
 * - Anisotropic stiffness: @f$ \int (A\nabla u) \cdot \nabla v \, dx @f$
 * - Jacobian forms: @f$ \int \mathbf{J}u : \mathbf{J}v \, dx @f$
 * - Potential operators for boundary elements
 *
 * These specializations mirror the generic quadrature rules in
 * `Rodin/Variational/QuadratureRule.h` but collapse to a single centroid
 * evaluation when possible. Mixed trial/test P1 spaces are supported for the
 * mass, stiffness, and Jacobian forms.
 *
 * @see <a href="class_rodin_1_1_variational_1_1_p1.html">P1</a>
 * @see <a href="_quadrature_formula_8h.html">QuadratureFormula</a>
 * @see <a href="_integral_8h.html">Integral</a>
 */
#ifndef RODIN_VARIATIONAL_P1_QUADRATURERULE_H
#define RODIN_VARIATIONAL_P1_QUADRATURERULE_H

#include "Rodin/FormLanguage/Traits.h"
#include "Rodin/Geometry/Region.h"
#include "Rodin/Math/Common.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Variational/IntegrationPoint.h"
#include "Rodin/Variational/ShapeFunction.h"
#include "Rodin/QF/Centroid.h"
#include "Rodin/QF/GaussLegendre.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"

#include "P1.h"
#include "Rodin/Math/Traits.h"


namespace Rodin::Variational
{
  /// @cond RODIN_DOXYGEN_INTERNAL
  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Integration of a P1 ShapeFunction.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int v \ dx \: ,
   * @f]
   * where @f$ v \in \mathbb{P}_1 @f$.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int v \ dx :
   * \texttt{QuadratureRule}}
   * {\vdash v : \mathbb{P}_1}
   * @f]
   */
  template <class NestedDerived, class Range, class Mesh>
  class QuadratureRule<
    ShapeFunctionBase<
      ShapeFunction<NestedDerived, P1<Range, Mesh>, TestSpace>, P1<Range, Mesh>, TestSpace>>
    : public LinearFormIntegratorBase<
        typename FormLanguage::Traits<
          ShapeFunctionBase<
            ShapeFunction<NestedDerived, P1<Range, Mesh>, TestSpace>,
            P1<Range, Mesh>,
            TestSpace>>
        ::ScalarType>
  {
    public:
      /// @brief Reports this handler as an optimized specialization.
      static constexpr bool Specialized = true;

      /// @brief Finite element space type.
      using FESType = P1<Range, Mesh>;

      /// @brief Integrand expression type.
      using IntegrandType =
        ShapeFunctionBase<ShapeFunction<NestedDerived, FESType, TestSpace>>;

      using IntegrandRangeType = typename FormLanguage::Traits<IntegrandType>::RangeType;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      /// @brief Parent class type.
      using Parent = LinearFormIntegratorBase<ScalarType>;

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_quadrature(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_quadrature(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(QuadratureRule&& other)
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

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        static_assert(std::is_same_v<IntegrandRangeType, ScalarType>);

        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& fes = integrand.getFiniteElementSpace();
        const auto& fe  = fes.getFiniteElement(d, idx);

        const size_t order = this->getOrder(polytope).value_or(
          integrand.getOrder(polytope).value_or(fe.getOrder()));

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

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
        assert(nte == fe.getCount());

        m_vec.resize(static_cast<Eigen::Index>(nte));
        m_vec.setZero();

        assert(m_quadrature);
        const auto& q = *m_quadrature;
        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());
          const auto& rc = m_qf->getPoint(qp);
          for (size_t local = 0; local < fe.getCount(); ++local)
            m_vec(local) += wdet * fe.getBasis(local)(rc);
        }

        return *this;
      }

      ScalarType integrate(size_t local) final override
      {
        return m_vec(local);
      }

      virtual Geometry::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      const Geometry::PolytopeQuadrature* m_quadrature;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Math::Vector<ScalarType> m_vec;
  };

  /**
   * @ingroup RodinCTAD
   */
  template <class NestedDerived, class Range, class Mesh>
  QuadratureRule(const ShapeFunctionBase<ShapeFunction<NestedDerived, P1<Range, Mesh>, TestSpace>>&)
    -> QuadratureRule<ShapeFunctionBase<ShapeFunction<NestedDerived, P1<Range, Mesh>, TestSpace>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Integration of the Dot product of some coefficient function and a
   * P1 ShapeFunction.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int f \cdot v \ dx \: ,
   * @f]
   * where @f$ v \in \mathbb{P}_1 @f$.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int f \cdot v \ dx :
   * \texttt{QuadratureRule}}
   * {\vdash v : \mathbb{P}_1}
   * @f]
   */
  template <class LHSDerived, class RHSDerived, class Range, class Mesh>
  class QuadratureRule<
    ShapeFunctionBase<
      Dot<
        FunctionBase<LHSDerived>,
        ShapeFunctionBase<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>, P1<Range, Mesh>, TestSpace>>,
          P1<Range, Mesh>, TestSpace>>
    : public LinearFormIntegratorBase<
        typename FormLanguage::Traits<
          ShapeFunctionBase<
            Dot<
              FunctionBase<LHSDerived>,
              ShapeFunctionBase<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>>
        ::ScalarType>
  {
    public:
      /// @brief Reports this handler as an optimized specialization.
      static constexpr bool Specialized = true;

      /// @brief Finite element space type.
      using FESType =
        P1<Range, Mesh>;

      /// @brief Left-hand side operand type.
      using LHSType =
        FunctionBase<LHSDerived>;

      /// @brief Right-hand side operand type.
      using RHSType =
        ShapeFunctionBase<ShapeFunction<RHSDerived, FESType, TestSpace>, FESType, TestSpace>;

      /// @brief Range type of the left-hand side operand.
      using LHSRangeType =
        typename FormLanguage::Traits<LHSType>::RangeType;

      /// @brief Range type of the right-hand side operand.
      using RHSRangeType =
        typename FormLanguage::Traits<RHSType>::RangeType;

      /// @brief Integrand expression type.
      using IntegrandType =
        ShapeFunctionBase<Dot<LHSType, RHSType>>;

      using IntegrandRangeType =
        typename FormLanguage::Traits<IntegrandType>::RangeType;

      /// @brief Scalar value type.
      using ScalarType =
        typename FormLanguage::Traits<IntegrandType>::ScalarType;

      /// @brief Parent class type.
      using Parent =
        LinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_quadrature(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_quadrature(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(QuadratureRule&& other)
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

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        static thread_local LHSRangeType s_v;
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& f   = integrand.getDerived().getLHS();
        const auto& fes = integrand.getFiniteElementSpace();
        const auto& fe  = fes.getFiniteElement(d, idx);

        const size_t order = this->getOrder(polytope).value_or(
          integrand.getOrder(polytope).value_or(fe.getOrder()));

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

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
        assert(nte == fe.getCount());

        m_vec.resize(static_cast<Eigen::Index>(nte));
        m_vec.setZero();

        if constexpr (std::is_same_v<RHSRangeType, ScalarType>)
        {
          assert(m_quadrature);
          const auto& q = *m_quadrature;
          for (size_t qp = 0; qp < q.getSize(); ++qp)
          {
            const auto& p = q.getPoint(qp);
            const IntegrationPoint ip(p, m_qf, qp);
            const ScalarType wdet =
              static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

            const auto& rc = m_qf->getPoint(qp);
            const ScalarType fval = f(ip);

            for (size_t local = 0; local < nte; ++local)
              m_vec(local) += wdet * fval * fe.getBasis(local)(rc);
          }
        }
        else if constexpr (FormLanguage::IsVectorRange<RHSRangeType>::Value)
        {
          assert(m_quadrature);
          const auto& q = *m_quadrature;
          for (size_t qp = 0; qp < q.getSize(); ++qp)
          {
            const auto& p = q.getPoint(qp);
            const IntegrationPoint ip(p, m_qf, qp);
            const ScalarType wdet =
              static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

            const auto& rc = m_qf->getPoint(qp);
            s_v = f(ip);

            for (size_t local = 0; local < nte; ++local)
              m_vec(local) += wdet * Math::dot(s_v, fe.getBasis(local)(rc));
          }
        }
        else
        {
          static_assert(
            std::is_same_v<RHSRangeType, ScalarType>
            || FormLanguage::IsVectorRange<RHSRangeType>::Value,
            "Unsupported P1 Integral(f.v) range type.");
        }

        return *this;
      }

      ScalarType integrate(size_t local) final override
      {
        return m_vec(local);
      }

      virtual Geometry::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      const Geometry::PolytopeQuadrature* m_quadrature;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Math::Vector<ScalarType> m_vec;
  };

  /**
   * @ingroup RodinCTAD
   */
  template <class LHSDerived, class RHSDerived, class Range, class Mesh>
  QuadratureRule(
      const ShapeFunctionBase<
        Dot<FunctionBase<LHSDerived>,
        ShapeFunctionBase<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>&)
    -> QuadratureRule<ShapeFunctionBase<
        Dot<
          FunctionBase<LHSDerived>,
          ShapeFunctionBase<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Integration of the isotropic Dot product of two instances of the
   * P1 ShapeFunction.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int u \cdot v \ dx \: ,
   * @f]
   * where @f$ u \in \mathbb{P}_1 @f$ and @f$ v \in \mathbb{P}_1 @f$.
   *
   * The specialization evaluates the integrand at the element centroid and
   * supports distinct P1 trial and test spaces.
   */
  template <
    class LHSDerived, class RHSDerived,
    class LHSRange, class RHSRange, class LHSMesh, class RHSMesh>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>,
        P1<LHSRange, LHSMesh>, TrialSpace>,
      ShapeFunctionBase<
        ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>,
        P1<RHSRange, RHSMesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<
            ShapeFunctionBase<
              ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>>,
            ShapeFunctionBase<
              ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>>>>::ScalarType>
  {
    public:
      /// @brief Reports this handler as an optimized specialization.
      static constexpr bool Specialized = true;

      using LHSFESType = P1<LHSRange, LHSMesh>;

      using RHSFESType = P1<RHSRange, RHSMesh>;

      /// @brief Left-hand side operand type.
      using LHSType =
        ShapeFunctionBase<ShapeFunction<LHSDerived, LHSFESType, TrialSpace>>;

      /// @brief Right-hand side operand type.
      using RHSType =
        ShapeFunctionBase<ShapeFunction<RHSDerived, RHSFESType, TestSpace>>;

      /// @brief Integrand expression type.
      using IntegrandType = Dot<LHSType, RHSType>;

      /// @brief Range type of the left-hand side operand.
      using LHSRangeType =
        typename FormLanguage::Traits<LHSType>::RangeType;

      /// @brief Range type of the right-hand side operand.
      using RHSRangeType =
        typename FormLanguage::Traits<RHSType>::RangeType;

      /// @brief Scalar value type.
      using ScalarType =
        typename FormLanguage::Traits<IntegrandType>::ScalarType;

      /// @brief Parent class type.
      using Parent =
        LocalBilinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);

      constexpr
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy())
      {}

      constexpr QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_basis(other.m_basis)
      {}

      constexpr QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_polytope(std::move(other.m_polytope)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::exchange(other.m_quadrature, nullptr)),
          m_matrix(std::move(other.m_matrix)),
          m_basis(std::move(other.m_basis)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::move(other.m_geometry))
      {}

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        return m_polytope.value().get();
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;
        const auto& geometry = polytope.getGeometry();

        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(
          polytope.getDimension(), polytope.getIndex());
        const auto& testfe = testfes.getFiniteElement(
          polytope.getDimension(), polytope.getIndex());

        // The integrand knows its own degree, coefficient included; the sum of
        // the two shape function orders is only the degree of the product of
        // the bases, and using it alone under-integrates any form carrying a
        // non-constant coefficient. The generic handler asks the integrand
        // first and this must agree with it, an optimized path being obliged
        // to compute the same thing faster, not something else.
        const size_t order =
          this->getOrder(polytope).value_or(getIntegrand().getOrder(polytope).value_or(
            trialfe.getOrder() + testfe.getOrder()));
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);

          // Reference P1 basis values depend on the geometry and quadrature
          // formula, but not on the physical mesh cell.
          const P1Element<ScalarType> scalarfe(geometry);
          const size_t nv = scalarfe.getCount();
          m_basis.resize(
            static_cast<Eigen::Index>(m_qf->getSize()), static_cast<Eigen::Index>(nv));
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
          {
            const auto& rc = m_qf->getPoint(qp);
            for (size_t a = 0; a < nv; ++a)
              m_basis(qp, a) = scalarfe.getBasis(a)(rc);
          }
        }

        assert(m_qf);
        m_quadrature = &polytope.getQuadrature(*m_qf);

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);

        assert(ntr == trialfe.getCount());
        assert(nte == testfe.getCount());

        const bool symmetric =
          (&trialfes.getMesh() == &testfes.getMesh()) && (ntr == nte);
        const size_t vdim = [&]() {
          if constexpr (FormLanguage::IsVectorRange<LHSRangeType>::Value)
          {
            const size_t trialVdim = trialfes.getVectorDimension();
            assert(trialVdim == testfes.getVectorDimension());
            return trialVdim;
          }
          else
          {
            return size_t(1);
          }
        }();
        assert(vdim > 0);
        assert(ntr % vdim == 0);
        assert(nte % vdim == 0);
        assert(static_cast<size_t>(m_basis.cols()) == ntr / vdim);
        assert(static_cast<size_t>(m_basis.cols()) == nte / vdim);

        m_matrix.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_matrix.setZero();

        assert(m_quadrature);
        const auto& q = *m_quadrature;
        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          if constexpr (std::is_same_v<LHSRangeType, ScalarType>)
          {
            if (symmetric)
            {
              for (size_t ib = 0; ib < nte; ++ib)
              {
                const ScalarType phi_te = m_basis(qp, ib);

                {
                  const ScalarType kii = wdet * phi_te * m_basis(qp, ib);
                  m_matrix(ib, ib) += kii;
                }

                for (size_t ia = 0; ia < ib; ++ia)
                {
                  const ScalarType kij = wdet * phi_te * m_basis(qp, ia);
                  m_matrix(ib, ia) += kij;
                }
              }
            }
            else
            {
              for (size_t ib = 0; ib < nte; ++ib)
              {
                const ScalarType phi_te = m_basis(qp, ib);
                for (size_t ia = 0; ia < ntr; ++ia)
                {
                  const ScalarType kij = wdet * phi_te * m_basis(qp, ia);
                  m_matrix(ib, ia) += kij;
                }
              }
            }
          }
          else if constexpr (FormLanguage::IsVectorRange<LHSRangeType>::Value)
          {
            if (symmetric)
            {
              for (size_t ib = 0; ib < nte; ++ib)
              {
                const size_t vb = ib / vdim;
                const size_t cb = ib % vdim;
                const ScalarType phi_te = m_basis(qp, vb);

                {
                  const ScalarType kii = wdet * phi_te * phi_te;
                  m_matrix(ib, ib) += kii;
                }

                for (size_t ia = 0; ia < ib; ++ia)
                {
                  if (ia % vdim == cb)
                  {
                    const ScalarType phi_tr = m_basis(qp, ia / vdim);
                    m_matrix(ib, ia) += wdet * phi_te * phi_tr;
                  }
                }
              }
            }
            else
            {
              for (size_t ib = 0; ib < nte; ++ib)
              {
                const size_t vb = ib / vdim;
                const size_t cb = ib % vdim;
                const ScalarType phi_te = m_basis(qp, vb);
                for (size_t ia = 0; ia < ntr; ++ia)
                {
                  if (ia % vdim == cb)
                  {
                    const ScalarType phi_tr = m_basis(qp, ia / vdim);
                    m_matrix(ib, ia) += wdet * phi_te * phi_tr;
                  }
                }
              }
            }
          }
        }

        if (symmetric)
        {
          m_matrix.template triangularView<Eigen::Upper>() =
            m_matrix.transpose().template triangularView<Eigen::Upper>();
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_matrix(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      const QF::QuadratureFormulaBase* m_qf = nullptr;
      const Geometry::PolytopeQuadrature* m_quadrature;

      Math::Matrix<ScalarType> m_matrix;
      Math::Matrix<ScalarType> m_basis;

      bool m_set = false;
      size_t m_order = 0;
      Optional<Geometry::Polytope::Type> m_geometry;
  };

  template <class LHSDerived, class RHSDerived, class Range, class Mesh>
  QuadratureRule(
    const Dot<
      ShapeFunctionBase<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>,
      ShapeFunctionBase<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>&)
    -> QuadratureRule<
        Dot<
          ShapeFunctionBase<
            ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>,
          ShapeFunctionBase<
            ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Integration of the anisotropic Dot product of two instances of the
   * P1 ShapeFunction.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int (A u) \cdot v \ dx \: ,
   * @f]
   * where @f$ u \in \mathbb{P}_1 @f$ and @f$ v \in \mathbb{P}_1 @f$, and @f$
   * A @f$ is a coefficient function.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int (A u) \cdot v \ dx :
   * \texttt{QuadratureRule}}
   * {\vdash u, v : \mathbb{P}_1}
   * @f]
   *
   * The specialization evaluates the integrand at the element centroid and
   * supports distinct P1 trial and test spaces while keeping the coefficient
   * evaluation to a single point. This specialization only applies to P1
   * finite element spaces by construction.
   */
  template <
    class CoefficientDerived, class LHSDerived, class RHSDerived,
    class LHSRange, class RHSRange, class LHSMesh, class RHSMesh>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Mult<
          FunctionBase<CoefficientDerived>,
          ShapeFunctionBase<
            ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>,
            P1<LHSRange, LHSMesh>, TrialSpace>>,
        P1<LHSRange, LHSMesh>, TrialSpace>,
      ShapeFunctionBase<
        ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>,
        P1<RHSRange, RHSMesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<
            ShapeFunctionBase<
              Mult<
                FunctionBase<CoefficientDerived>,
                ShapeFunctionBase<
                  ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>>>>,
            ShapeFunctionBase<
              ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>>>>
        ::ScalarType>
  {
    public:
      /// @brief Reports this handler as an optimized specialization.
      static constexpr bool Specialized = true;

      using LHSFESType = P1<LHSRange, LHSMesh>;

      using RHSFESType = P1<RHSRange, RHSMesh>;

      using CoefficientType = FunctionBase<CoefficientDerived>;

      using MultiplicandType =
        ShapeFunctionBase<ShapeFunction<LHSDerived, LHSFESType, TrialSpace>>;

      /// @brief Left-hand side operand type.
      using LHSType =
        ShapeFunctionBase<Mult<CoefficientType, MultiplicandType>>;

      /// @brief Right-hand side operand type.
      using RHSType =
        ShapeFunctionBase<ShapeFunction<RHSDerived, RHSFESType, TestSpace>>;

      /// @brief Integrand expression type.
      using IntegrandType = Dot<LHSType, RHSType>;

      using CoefficientRangeType =
        typename FormLanguage::Traits<CoefficientType>::RangeType;

      using MultiplicandRangeType =
        typename FormLanguage::Traits<MultiplicandType>::RangeType;

      /// @brief Range type of the left-hand side operand.
      using LHSRangeType =
        typename FormLanguage::Traits<LHSType>::RangeType;

      /// @brief Range type of the right-hand side operand.
      using RHSRangeType =
        typename FormLanguage::Traits<RHSType>::RangeType;

      /// @brief Scalar value type.
      using ScalarType =
        typename FormLanguage::Traits<IntegrandType>::ScalarType;

      /// @brief Parent class type.
      using Parent =
        LocalBilinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_quadrature(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_quadrature(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::exchange(other.m_quadrature, nullptr)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_matrix(std::move(other.m_matrix)),
          m_basis(std::move(other.m_basis))
      {}

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& coeff = lhs.getDerived().getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes.getFiniteElement(d, idx);

        const size_t order =
          this->getOrder(polytope).value_or(getIntegrand().getOrder(polytope).value_or(
            trialfe.getOrder() + testfe.getOrder()));

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);

          // Reference P1 basis values are shared by every physical cell with
          // this geometry and quadrature formula.
          const P1Element<ScalarType> scalarfe(geometry);
          const size_t nv = scalarfe.getCount();
          m_basis.resize(
            static_cast<Eigen::Index>(m_qf->getSize()), static_cast<Eigen::Index>(nv));
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
          {
            const auto& rc = m_qf->getPoint(qp);
            for (size_t a = 0; a < nv; ++a)
              m_basis(qp, a) = scalarfe.getBasis(a)(rc);
          }
        }

        assert(m_qf);
        m_quadrature = &polytope.getQuadrature(*m_qf);

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);

        assert(ntr == trialfe.getCount());
        assert(nte == testfe.getCount());

        const size_t vdim = [&]() {
          if constexpr (FormLanguage::IsVectorRange<MultiplicandRangeType>::Value)
          {
            const size_t trialVdim = trialfes.getVectorDimension();
            assert(trialVdim == testfes.getVectorDimension());
            return trialVdim;
          }
          else
          {
            return size_t(1);
          }
        }();
        assert(vdim > 0);
        assert(ntr % vdim == 0);
        assert(nte % vdim == 0);
        assert(static_cast<size_t>(m_basis.cols()) == ntr / vdim);
        assert(static_cast<size_t>(m_basis.cols()) == nte / vdim);

        m_matrix.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_matrix.setZero();

        assert(m_quadrature);
        const auto& q = *m_quadrature;
        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const IntegrationPoint ip(p, m_qf, qp);
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          if constexpr (std::is_same_v<CoefficientRangeType, ScalarType>)
          {
            const ScalarType csv = coeff.getValue(ip);

            if constexpr (std::is_same_v<MultiplicandRangeType, ScalarType>)
            {
              for (size_t ib = 0; ib < nte; ++ib)
              {
                const ScalarType phi_te = m_basis(qp, ib);
                for (size_t ia = 0; ia < ntr; ++ia)
                {
                  const ScalarType kij = wdet * csv * phi_te * m_basis(qp, ia);
                  m_matrix(ib, ia) += kij;
                }
              }
            }
            else
            {
              for (size_t ib = 0; ib < nte; ++ib)
              {
                const size_t vb = ib / vdim;
                const size_t cb = ib % vdim;
                const ScalarType phi_te = m_basis(qp, vb);
                for (size_t ia = 0; ia < ntr; ++ia)
                {
                  if (ia % vdim == cb)
                  {
                    const ScalarType phi_tr = m_basis(qp, ia / vdim);
                    m_matrix(ib, ia) += wdet * csv * phi_te * phi_tr;
                  }
                }
              }
            }
          }
          else if constexpr (FormLanguage::IsMatrixRange<CoefficientRangeType>::Value)
          {
            static_assert(FormLanguage::IsVectorRange<MultiplicandRangeType>::Value);
            static_assert(FormLanguage::IsVectorRange<RHSRangeType>::Value);

            Math::SpatialMatrix<ScalarType> cmv;
            cmv = coeff.getValue(ip);

            for (size_t ib = 0; ib < nte; ++ib)
            {
              const size_t vb = ib / vdim;
              const size_t cb = ib % vdim;
              const ScalarType phi_te = m_basis(qp, vb);
              for (size_t ia = 0; ia < ntr; ++ia)
              {
                const size_t ca = ia % vdim;
                const ScalarType phi_tr = m_basis(qp, ia / vdim);
                m_matrix(ib, ia) += wdet * phi_te * cmv(cb, ca) * phi_tr;
              }
            }
          }
          else
          {
            assert(false);
          }
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_matrix(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      const Geometry::PolytopeQuadrature* m_quadrature;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Math::Matrix<ScalarType> m_matrix;
      Math::Matrix<ScalarType> m_basis;
  };

  template <class CoefficientDerived, class LHSDerived, class RHSDerived, class Number, class Mesh>
  QuadratureRule(const
    Dot<
      ShapeFunctionBase<
        Mult<
          FunctionBase<CoefficientDerived>,
          ShapeFunctionBase<ShapeFunction<LHSDerived, P1<Number, Mesh>, TrialSpace>>>>,
      ShapeFunctionBase<
        ShapeFunction<RHSDerived, P1<Number, Mesh>, TestSpace>>>&)
  ->
  QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Mult<
          FunctionBase<CoefficientDerived>,
          ShapeFunctionBase<ShapeFunction<LHSDerived, P1<Number, Mesh>, TrialSpace>>>>,
      ShapeFunctionBase<

        ShapeFunction<RHSDerived, P1<Number, Mesh>, TestSpace>>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Integration of the isotropic Dot product of two instances of the P1
   * Grad of ShapeFunction.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int \nabla u \cdot \nabla v \ dx \: ,
   * @f]
   * where @f$ u \in \mathbb{P}_1 @f$ and @f$ v \in \mathbb{P}_1 @f$, and @f$ A
   * @f$ is a coefficient function.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int \nabla u \cdot \nabla v \ dx :
   * \texttt{QuadratureRule}}
   * {\vdash u, v : \mathbb{P}_1}
   * @f]
   *
   * A single centroid evaluation is used and mixed P1 trial/test spaces are
   * accommodated.
   */
  template <class LHSDerived, class RHSDerived, class LHSRange, class RHSRange, class LHSMesh, class RHSMesh>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Grad<ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>>, P1<LHSRange, LHSMesh>, TrialSpace>,
      ShapeFunctionBase<
        Grad<ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>>, P1<RHSRange, RHSMesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<
            ShapeFunctionBase<
              Grad<ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>>, P1<LHSRange, LHSMesh>, TrialSpace>,
            ShapeFunctionBase<
              Grad<ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>>, P1<RHSRange, RHSMesh>, TestSpace>>>
        ::ScalarType>
  {
    public:
      /// @brief Reports this handler as an optimized specialization.
      static constexpr bool Specialized = true;

      using LHSFESType = P1<LHSRange, LHSMesh>;

      using RHSFESType = P1<RHSRange, RHSMesh>;

      /// @brief Left-hand side operand type.
      using LHSType =
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, LHSFESType, TrialSpace>>>;

      using LHSOperandType =
        ShapeFunction<LHSDerived, LHSFESType, TrialSpace>;

      using LHSOperandRangeType =
        typename FormLanguage::Traits<LHSOperandType>::RangeType;

      /// @brief Right-hand side operand type.
      using RHSType =
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, RHSFESType, TestSpace>>>;

      using RHSOperandType =
        ShapeFunction<RHSDerived, RHSFESType, TestSpace>;

      using RHSOperandRangeType =
        typename FormLanguage::Traits<RHSOperandType>::RangeType;

      /// @brief Integrand expression type.
      using IntegrandType = Dot<LHSType, RHSType>;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      /// @brief Parent class type.
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSOperandRangeType, ScalarType>);
      static_assert(std::is_same_v<RHSOperandRangeType, ScalarType>);

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_quadrature(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_quadrature(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::exchange(other.m_quadrature, nullptr)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_refGrad(std::move(other.m_refGrad)),
          m_matrix(std::move(other.m_matrix))
      {}

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes  = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes.getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe.getOrder();
        const size_t order =
          this->getOrder(polytope).value_or(getIntegrand().getOrder(polytope).value_or(
            ((k_tr == 0 || k_te == 0) ? 0 : (k_tr + k_te - 2))));

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);

          const size_t n = trialfe.getCount();

          // Cached per quadrature point, not at the first one: reference
          // gradients are constant only on a simplex, so freezing them at one
          // point is exact only there, or when the rule has a single point.
          const size_t nqp = m_qf->getSize();
          m_refGrad.assign(nqp, {});
          for (size_t qp = 0; qp < nqp; ++qp)
          {
            const auto& rc = m_qf->getPoint(qp);
            m_refGrad[qp].resize(n);
            for (size_t local = 0; local < n; ++local)
            {
              auto& g = m_refGrad[qp][local];
              g.resize(d);
              const auto& basis = trialfe.getBasis(local);
              for (size_t j = 0; j < d; ++j)
                g(j) = basis.template getDerivative<1>(j)(rc);
            }
          }

          m_matrix.resize(n, n);
        }

        assert(m_qf);
        m_quadrature = &polytope.getQuadrature(*m_qf);

        const size_t n = m_refGrad.empty() ? 0 : m_refGrad.front().size();

        m_matrix.setZero();

        assert(m_quadrature);
        const auto& q = *m_quadrature;
        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto& Jinv = p.getJacobianInverse();
          const auto G = Jinv * Jinv.transpose();

          const auto& refGrad = m_refGrad[qp];
          for (size_t i = 0; i < n; ++i)
          {
            const auto Ggi = G * refGrad[i];
            m_matrix(i, i) += wdet * Math::dot(refGrad[i], Ggi);
            for (size_t j = 0; j < i; ++j)
              m_matrix(i, j) += wdet * Math::dot(refGrad[j], Ggi);
          }
        }

        m_matrix.template triangularView<Eigen::Upper>() =
          m_matrix.template triangularView<Eigen::Lower>().transpose();

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_matrix(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      const Geometry::PolytopeQuadrature* m_quadrature;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      /// @brief Reference gradients per quadrature point and dof.
      std::vector<std::vector<Math::SpatialVector<ScalarType>>> m_refGrad;
      Math::Matrix<ScalarType> m_matrix;
  };

  template <class LHSDerived, class RHSDerived, class Range, class Mesh>
  QuadratureRule(
      const Dot<
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>>,
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>&)
    -> QuadratureRule<Dot<
          ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>>,
          ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Integration of the anisotropic Dot product of two instances of the
   * P1 Grad of ShapeFunction.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int (A \nabla u) \cdot \nabla v \ dx \: ,
   * @f]
   * where @f$ u \in \mathbb{P}_1 @f$ and @f$ v \in \mathbb{P}_1 @f$, and @f$ A
   * @f$ is a coefficient function.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int (A \nabla u) \cdot \nabla v \ dx :
   * \texttt{QuadratureRule}}
   * {\vdash u, v : \mathbb{P}_1}
   * @f]
   *
   * The rule evaluates the coefficient and gradients at the centroid and
   * permits distinct P1 trial and test spaces.
   */
  template <
    class CoefficientDerived, class LHSDerived, class RHSDerived,
    class LHSRange, class RHSRange, class LHSMesh, class RHSMesh>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Mult<
          FunctionBase<CoefficientDerived>,
          ShapeFunctionBase<
            Grad<ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>>,
            P1<LHSRange, LHSMesh>, TrialSpace>>,
        P1<LHSRange, LHSMesh>, TrialSpace>,
      ShapeFunctionBase<
        Grad<ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>>,
        P1<RHSRange, RHSMesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<
            ShapeFunctionBase<
              Mult<
                FunctionBase<CoefficientDerived>,
                ShapeFunctionBase<
                  Grad<ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>>,
                  P1<LHSRange, LHSMesh>, TrialSpace>>,
              P1<LHSRange, LHSMesh>, TrialSpace>,
            ShapeFunctionBase<
              Grad<ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>>,
              P1<RHSRange, RHSMesh>, TestSpace>>>
        ::ScalarType>
  {
    public:
      /// @brief Reports this handler as an optimized specialization.
      static constexpr bool Specialized = true;

      using LHSFESType = P1<LHSRange, LHSMesh>;

      using RHSFESType = P1<RHSRange, RHSMesh>;

      using CoefficientType = FunctionBase<CoefficientDerived>;

      using MultiplicandType =
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, LHSFESType, TrialSpace>>>;

      using MultiplicandOperandType =
        ShapeFunction<LHSDerived, LHSFESType, TrialSpace>;

      using CoefficientRangeType =
        typename FormLanguage::Traits<CoefficientType>::RangeType;

      using MultiplicandRangeType =
        typename FormLanguage::Traits<MultiplicandType>::RangeType;

      /// @brief Left-hand side operand type.
      using LHSType =
        ShapeFunctionBase<Mult<CoefficientType, MultiplicandType>>;

      using MultiplicandOperandRangeType =
        typename FormLanguage::Traits<MultiplicandOperandType>::RangeType;

      /// @brief Right-hand side operand type.
      using RHSType =
        ShapeFunctionBase<
          Grad<ShapeFunction<RHSDerived, RHSFESType, TestSpace>>>;

      using RHSOperandType =
        ShapeFunction<RHSDerived, RHSFESType, TestSpace>;

      using RHSOperandRangeType =
        typename FormLanguage::Traits<RHSOperandType>::RangeType;

      /// @brief Integrand expression type.
      using IntegrandType = Dot<LHSType, RHSType>;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      /// @brief Parent class type.
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<MultiplicandOperandRangeType, ScalarType>);
      static_assert(std::is_same_v<RHSOperandRangeType, ScalarType>);

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_quadrature(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_quadrature(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::exchange(other.m_quadrature, nullptr)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_trialRefGrad(std::move(other.m_trialRefGrad)),
          m_testRefGrad(std::move(other.m_testRefGrad)),
          m_matrix(std::move(other.m_matrix))
      {}

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& coeff = lhs.getDerived().getLHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes.getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe.getOrder();
        const size_t order =
          this->getOrder(polytope).value_or(getIntegrand().getOrder(polytope).value_or(
            ((k_tr == 0 || k_te == 0) ? 0 : (k_tr + k_te - 2))));

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);

          // Cached per quadrature point, not at the first one: reference
          // gradients are constant only on a simplex, so freezing them at one
          // point is exact only there, or when the rule has a single point.
          const size_t nqp = m_qf->getSize();
          m_trialRefGrad.assign(nqp, {});
          m_testRefGrad.assign(nqp, {});
          for (size_t qp = 0; qp < nqp; ++qp)
          {
            const auto& rc = m_qf->getPoint(qp);

            m_trialRefGrad[qp].resize(trialfe.getCount());
            for (size_t local = 0; local < trialfe.getCount(); ++local)
            {
              auto& g = m_trialRefGrad[qp][local];
              g.resize(d);
              const auto& basis = trialfe.getBasis(local);
              for (size_t j = 0; j < d; ++j)
                g(j) = basis.template getDerivative<1>(j)(rc);
            }

            if (trialfes == testfes)
            {
              m_testRefGrad[qp] = m_trialRefGrad[qp];
            }
            else
            {
              m_testRefGrad[qp].resize(testfe.getCount());
              for (size_t local = 0; local < testfe.getCount(); ++local)
              {
                auto& g = m_testRefGrad[qp][local];
                g.resize(d);
                const auto& basis = testfe.getBasis(local);
                for (size_t j = 0; j < d; ++j)
                  g(j) = basis.template getDerivative<1>(j)(rc);
              }
            }
          }

          m_matrix.resize(testfe.getCount(), trialfe.getCount());
        }

        assert(m_qf);
        m_quadrature = &polytope.getQuadrature(*m_qf);

        m_matrix.setZero();

        assert(m_quadrature);
        const auto& q = *m_quadrature;
        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const IntegrationPoint ip(p, m_qf, qp);
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto& Jinv = p.getJacobianInverse();
          const auto G = Jinv * Jinv.transpose();
          const auto& trialRefGrad = m_trialRefGrad[qp];
          const auto& testRefGrad = m_testRefGrad[qp];

          if constexpr (std::is_same_v<CoefficientRangeType, ScalarType>)
          {
            const ScalarType csv = coeff.getValue(ip);

            if (trialfes == testfes)
            {
              const size_t n = trialRefGrad.size();
              for (size_t i = 0; i < n; ++i)
              {
                const auto Ggi = G * trialRefGrad[i];
                m_matrix(i, i) += wdet * csv * Math::dot(trialRefGrad[i], Ggi);
                for (size_t j = 0; j < i; ++j)
                  m_matrix(i, j) += wdet * csv * Math::dot(trialRefGrad[j], Ggi);
              }
            }
            else
            {
              const size_t ntr = trialRefGrad.size();
              const size_t nte = testRefGrad.size();

              for (size_t te = 0; te < nte; ++te)
              {
                const auto Ggte = G * testRefGrad[te];
                for (size_t tr = 0; tr < ntr; ++tr)
                  m_matrix(te, tr) += wdet * csv * Math::dot(trialRefGrad[tr], Ggte);
              }
            }
          }
          else if constexpr (FormLanguage::IsMatrixRange<CoefficientRangeType>::Value)
          {
            Math::SpatialMatrix<ScalarType> cmv;
            cmv = coeff.getValue(ip);

            if (trialfes == testfes)
            {
              const size_t n = trialRefGrad.size();
              for (size_t i = 0; i < n; ++i)
              {
                const auto AGgi = cmv * (G * trialRefGrad[i]);
                m_matrix(i, i) += wdet * Math::dot(AGgi, trialRefGrad[i]);
                for (size_t j = 0; j < i; ++j)
                  m_matrix(i, j) += wdet * Math::dot(AGgi, trialRefGrad[j]);
              }

              for (size_t i = 0; i < n; ++i)
                for (size_t j = i + 1; j < n; ++j)
                  m_matrix(i, j) +=
                    wdet * Math::dot(cmv * (G * trialRefGrad[j]), trialRefGrad[i]);
            }
            else
            {
              const size_t ntr = trialRefGrad.size();
              const size_t nte = testRefGrad.size();

              for (size_t te = 0; te < nte; ++te)
                for (size_t tr = 0; tr < ntr; ++tr)
                  m_matrix(te, tr) +=
                    wdet * Math::dot(cmv * (G * trialRefGrad[tr]), testRefGrad[te]);
            }
          }
          else
          {
            assert(false);
          }
        }

        if constexpr (std::is_same_v<CoefficientRangeType, ScalarType>)
        {
          if (trialfes == testfes)
          {
            m_matrix.template triangularView<Eigen::Upper>() =
              m_matrix.template triangularView<Eigen::Lower>().transpose();
          }
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_matrix(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      const Geometry::PolytopeQuadrature* m_quadrature;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      /// @brief Trial reference gradients per quadrature point and dof.
      std::vector<std::vector<Math::SpatialVector<ScalarType>>> m_trialRefGrad;
      /// @brief Test reference gradients per quadrature point and dof.
      std::vector<std::vector<Math::SpatialVector<ScalarType>>> m_testRefGrad;

      Math::Matrix<ScalarType> m_matrix;
  };

  /**
   * @ingroup RodinCTAD
   */
  template <class LHSFunctionDerived, class LHSDerived, class RHSDerived, class Range, class Mesh>
  QuadratureRule(
      const Dot<
        ShapeFunctionBase<
          Mult<
            FunctionBase<LHSFunctionDerived>,
            ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>>>>,
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>&)
    -> QuadratureRule<
          Dot<ShapeFunctionBase<
           Mult<
             FunctionBase<LHSFunctionDerived>,
             ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>>>>,
          ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Specialization for @f$\int c \, (u \cdot v) \ dx@f$ in the case of P1 shape functions.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int c \, (u \cdot v) \ dx \: ,
   * @f]
   * where @f$ u, v \in \mathbb{P}_1 @f$ and @f$ c @f$ is a coefficient function.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int c \, (u \cdot v) \ dx : \texttt{QuadratureRule}}
   * {\vdash u, v : \mathbb{P}_1}
   * @f]
   */
  template <
    class CoefficientDerived, class LHSDerived, class RHSDerived,
    class LHSRange, class RHSRange, class LHSMesh, class RHSMesh>
  class QuadratureRule<
    Mult<
      FunctionBase<CoefficientDerived>,
      Dot<
        ShapeFunctionBase<
          ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>,
          P1<LHSRange, LHSMesh>, TrialSpace>,
        ShapeFunctionBase<
          ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>,
          P1<RHSRange, RHSMesh>, TestSpace>>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<
            ShapeFunctionBase<
              ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>>,
            ShapeFunctionBase<
              ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>>>>::ScalarType>
  {
    public:
      /// @brief Reports this handler as an optimized specialization.
      static constexpr bool Specialized = true;

      using LHSFESType = P1<LHSRange, LHSMesh>;
      using RHSFESType = P1<RHSRange, RHSMesh>;

      using CoefficientType = FunctionBase<CoefficientDerived>;

      /// @brief Left-hand side operand type.
      using LHSType =
        ShapeFunctionBase<ShapeFunction<LHSDerived, LHSFESType, TrialSpace>>;

      /// @brief Right-hand side operand type.
      using RHSType =
        ShapeFunctionBase<ShapeFunction<RHSDerived, RHSFESType, TestSpace>>;

      using InnerIntegrandType = Dot<LHSType, RHSType>;
      /// @brief Integrand expression type.
      using IntegrandType = Mult<CoefficientType, InnerIntegrandType>;

      /// @brief Range type of the left-hand side operand.
      using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;
      /// @brief Range type of the right-hand side operand.
      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;
      using ScalarType   = typename FormLanguage::Traits<InnerIntegrandType>::ScalarType;

      /// @brief Parent class type.
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getRHS().getLHS().getLeaf(), integrand.getRHS().getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_quadrature(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_quadrature(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::exchange(other.m_quadrature, nullptr)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_matrix(std::move(other.m_matrix))
      {}

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        const auto& integrand = getIntegrand();
        const auto& coeff = integrand.getLHS();
        const auto& inner = integrand.getRHS();
        const auto& lhs = inner.getLHS();
        const auto& rhs = inner.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes.getFiniteElement(d, idx);

        const size_t order =
          this->getOrder(polytope).value_or(getIntegrand().getOrder(polytope).value_or(
            trialfe.getOrder() + testfe.getOrder()));

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);
        }

        assert(m_qf);
        m_quadrature = &polytope.getQuadrature(*m_qf);

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);

        assert(ntr == trialfe.getCount());
        assert(nte == testfe.getCount());

        m_matrix.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_matrix.setZero();

        assert(m_quadrature);
        const auto& q = *m_quadrature;
        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const IntegrationPoint ip(p, m_qf, qp);
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto& rc = m_qf->getPoint(qp);

          const ScalarType csv = coeff(ip);

          if constexpr (std::is_same_v<LHSRangeType, ScalarType>)
          {
            for (size_t ib = 0; ib < nte; ++ib)
            {
              const ScalarType phi_te = testfe.getBasis(ib)(rc);
              for (size_t ia = 0; ia < ntr; ++ia)
              {
                const ScalarType kij =
                  wdet * csv * phi_te * trialfe.getBasis(ia)(rc);
                m_matrix(ib, ia) += kij;
              }
            }
          }
          else if constexpr (FormLanguage::IsVectorRange<LHSRangeType>::Value)
          {
            for (size_t ib = 0; ib < nte; ++ib)
            {
              const auto phi_te = testfe.getBasis(ib)(rc);
              for (size_t ia = 0; ia < ntr; ++ia)
              {
                const auto phi_tr = trialfe.getBasis(ia)(rc);
                const ScalarType kij = wdet * csv * Math::dot(phi_te, phi_tr);
                m_matrix(ib, ia) += kij;
              }
            }
          }
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_matrix(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      const Geometry::PolytopeQuadrature* m_quadrature;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Math::Matrix<ScalarType> m_matrix;
  };

  template <class CoefficientDerived, class LHSDerived, class RHSDerived, class Range, class Mesh>
  QuadratureRule(
    const Mult<
      FunctionBase<CoefficientDerived>,
      Dot<
        ShapeFunctionBase<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>,
        ShapeFunctionBase<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>&)
    -> QuadratureRule<
        Mult<
          FunctionBase<CoefficientDerived>,
          Dot<
            ShapeFunctionBase<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>,
            ShapeFunctionBase<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Specialization for @f$\int (\nabla \cdot u)\, q \ dx@f$ in the case of P1 shape functions.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int (\nabla \cdot u)\, q \ dx \: ,
   * @f]
   * where @f$ u \in \mathbb{P}_1 @f$ (vector-valued) and @f$ q \in \mathbb{P}_1 @f$ (scalar-valued).
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int (\nabla \cdot u)\, q \ dx : \texttt{QuadratureRule}}
   * {\vdash u, q : \mathbb{P}_1}
   * @f]
   */
  template <
    class LHSDerived, class RHSDerived,
    class LHSMesh, class RHSMesh>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Div<ShapeFunction<LHSDerived, P1<Math::SpatialVector<Real>, LHSMesh>, TrialSpace>>,
        P1<Math::SpatialVector<Real>, LHSMesh>, TrialSpace>,
      ShapeFunctionBase<
        ShapeFunction<RHSDerived, P1<Real, RHSMesh>, TestSpace>,
        P1<Real, RHSMesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<typename FormLanguage::Traits<P1<Real, LHSMesh>>::ScalarType>
  {
    public:
      /// @brief Reports this handler as an optimized specialization.
      static constexpr bool Specialized = true;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<P1<Real, LHSMesh>>::ScalarType;
      using TrialFESType = P1<Math::SpatialVector<Real>, LHSMesh>;
      using TestFESType  = P1<Real, RHSMesh>;

      /// @brief Left-hand side operand type.
      using LHSType = ShapeFunctionBase<
        Div<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>>,
        TrialFESType, TrialSpace>;

      /// @brief Right-hand side operand type.
      using RHSType = ShapeFunctionBase<
        ShapeFunction<RHSDerived, TestFESType, TestSpace>,
        TestFESType, TestSpace>;

      /// @brief Integrand expression type.
      using IntegrandType = Dot<LHSType, RHSType>;
      /// @brief Parent class type.
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_quadrature(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_quadrature(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::exchange(other.m_quadrature, nullptr)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_refGrad(std::move(other.m_refGrad)),
          m_testBasis(std::move(other.m_testBasis)),
          m_matrix(std::move(other.m_matrix))
      {}

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes  = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes.getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe.getOrder();
        const size_t order =
          this->getOrder(polytope).value_or(getIntegrand().getOrder(polytope).value_or(
            ((k_tr == 0 || k_te == 0) ? 0 : (k_tr + k_te - 2))));

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);

          // Cached per quadrature point, not at the first one: the basis
          // values vary within every element, and the reference gradients
          // vary within every element that is not a simplex. Freezing either
          // at one point is exact only when the rule has a single point.
          const size_t nqp = m_qf->getSize();

          m_testBasis.assign(nqp, {});
          m_refGrad.assign(nqp, {});
          for (size_t qp = 0; qp < nqp; ++qp)
          {
            const auto& rc = m_qf->getPoint(qp);

            m_testBasis[qp].resize(testfe.getCount());
            for (size_t i = 0; i < testfe.getCount(); ++i)
              m_testBasis[qp][i] = testfe.getBasis(i)(rc);

            m_refGrad[qp].resize(trialfe.getCount());
            for (size_t i = 0; i < trialfe.getCount(); ++i)
            {
              auto& refg = m_refGrad[qp][i];
              refg.resize(trialfes.getVectorDimension());
              for (size_t comp = 0; comp < trialfes.getVectorDimension(); ++comp)
              {
                refg[comp].resize(d);
                const auto& basis = trialfe.getBasis(i);
                for (size_t j = 0; j < d; ++j)
                  refg[comp](j) = basis.template getDerivative<1>(comp, j)(rc);
              }
            }
          }

          m_matrix.resize(testfe.getCount(), trialfe.getCount());
        }

        assert(m_qf);
        m_quadrature = &polytope.getQuadrature(*m_qf);

        m_matrix.setZero();

        assert(m_quadrature);
        const auto& q = *m_quadrature;
        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto& Jinv = p.getJacobianInverse();
          const auto& refGrad = m_refGrad[qp];
          const auto& testBasis = m_testBasis[qp];
          const size_t vdim = refGrad.empty() ? 0 : refGrad.front().size();

          for (size_t i = 0; i < refGrad.size(); ++i)
          {
            ScalarType div = 0;
            for (size_t comp = 0; comp < std::min(vdim, d); ++comp)
            {
              ScalarType physComp = 0;
              for (size_t j = 0; j < d; ++j)
                physComp += Jinv(comp, j) * refGrad[i][comp](j);
              div += physComp;
            }

            const size_t nte = testBasis.size();
            for (size_t te = 0; te < nte; ++te)
              m_matrix(te, i) += wdet * testBasis[te] * div;
          }
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_matrix(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      const Geometry::PolytopeQuadrature* m_quadrature;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      /// @brief Reference gradients per quadrature point, dof and component.
      std::vector<std::vector<std::vector<Math::SpatialVector<ScalarType>>>> m_refGrad;
      /// @brief Test basis values per quadrature point and dof.
      std::vector<std::vector<ScalarType>> m_testBasis;
      Math::Matrix<ScalarType> m_matrix;
  };

  template <class LHSDerived, class RHSDerived, class LHSMesh, class RHSMesh>
  using P1DivTrialIntegrand =
    Dot<
      ShapeFunctionBase<
        Div<ShapeFunction<LHSDerived, P1<Math::SpatialVector<Real>, LHSMesh>, TrialSpace>>,
        P1<Math::SpatialVector<Real>, LHSMesh>, TrialSpace>,
      ShapeFunctionBase<
        ShapeFunction<RHSDerived, P1<Real, RHSMesh>, TestSpace>,
        P1<Real, RHSMesh>, TestSpace>>;

  template <class LHSDerived, class RHSDerived, class LHSMesh, class RHSMesh>
  QuadratureRule(const P1DivTrialIntegrand<LHSDerived, RHSDerived, LHSMesh, RHSMesh>&)
    -> QuadratureRule<P1DivTrialIntegrand<LHSDerived, RHSDerived, LHSMesh, RHSMesh>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Specialization for @f$\int p \, (\nabla \cdot v)\, dx@f$ in the case of P1 shape functions.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int p \, (\nabla \cdot v)\, dx \: ,
   * @f]
   * where @f$ p \in \mathbb{P}_1 @f$ (scalar-valued) and @f$ v \in \mathbb{P}_1 @f$ (vector-valued).
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int p \, (\nabla \cdot v)\, dx : \texttt{QuadratureRule}}
   * {\vdash p, v : \mathbb{P}_1}
   * @f]
   */
  template <
    class LHSDerived, class RHSDerived,
    class LHSMesh, class RHSMesh>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        ShapeFunction<LHSDerived, P1<Real, LHSMesh>, TrialSpace>,
        P1<Real, LHSMesh>, TrialSpace>,
      ShapeFunctionBase<
        Div<ShapeFunction<RHSDerived, P1<Math::SpatialVector<Real>, RHSMesh>, TestSpace>>,
        P1<Math::SpatialVector<Real>, RHSMesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<typename FormLanguage::Traits<P1<Real, LHSMesh>>::ScalarType>
  {
    public:
      /// @brief Reports this handler as an optimized specialization.
      static constexpr bool Specialized = true;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<P1<Real, LHSMesh>>::ScalarType;
      using TrialFESType = P1<Real, LHSMesh>;
      using TestFESType  = P1<Math::SpatialVector<Real>, RHSMesh>;

      /// @brief Left-hand side operand type.
      using LHSType = ShapeFunctionBase<
        ShapeFunction<LHSDerived, TrialFESType, TrialSpace>,
        TrialFESType, TrialSpace>;

      /// @brief Right-hand side operand type.
      using RHSType = ShapeFunctionBase<
        Div<ShapeFunction<RHSDerived, TestFESType, TestSpace>>,
        TestFESType, TestSpace>;

      /// @brief Integrand expression type.
      using IntegrandType = Dot<LHSType, RHSType>;
      /// @brief Parent class type.
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_quadrature(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_quadrature(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::exchange(other.m_quadrature, nullptr)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_refGrad(std::move(other.m_refGrad)),
          m_trialBasis(std::move(other.m_trialBasis)),
          m_matrix(std::move(other.m_matrix))
      {}

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes  = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes.getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe.getOrder();
        const size_t order =
          this->getOrder(polytope).value_or(getIntegrand().getOrder(polytope).value_or(
            ((k_tr == 0 || k_te == 0) ? 0 : (k_tr + k_te - 2))));

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);

          // Cached per quadrature point, not at the first one: the basis
          // values vary within every element, and the reference gradients
          // vary within every element that is not a simplex. Freezing either
          // at one point is exact only when the rule has a single point.
          const size_t nqp = m_qf->getSize();

          m_trialBasis.assign(nqp, {});
          m_refGrad.assign(nqp, {});
          for (size_t qp = 0; qp < nqp; ++qp)
          {
            const auto& rc = m_qf->getPoint(qp);

            m_trialBasis[qp].resize(trialfe.getCount());
            for (size_t i = 0; i < trialfe.getCount(); ++i)
              m_trialBasis[qp][i] = trialfe.getBasis(i)(rc);

            m_refGrad[qp].resize(testfe.getCount());
            for (size_t i = 0; i < testfe.getCount(); ++i)
            {
              auto& refg = m_refGrad[qp][i];
              refg.resize(testfes.getVectorDimension());
              for (size_t comp = 0; comp < testfes.getVectorDimension(); ++comp)
              {
                refg[comp].resize(d);
                const auto& basis = testfe.getBasis(i);
                for (size_t j = 0; j < d; ++j)
                  refg[comp](j) = basis.template getDerivative<1>(comp, j)(rc);
              }
            }
          }

          m_matrix.resize(testfe.getCount(), trialfe.getCount());
        }

        assert(m_qf);
        m_quadrature = &polytope.getQuadrature(*m_qf);

        m_matrix.setZero();

        assert(m_quadrature);
        const auto& q = *m_quadrature;
        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto& Jinv = p.getJacobianInverse();
          const auto& refGrad = m_refGrad[qp];
          const auto& trialBasis = m_trialBasis[qp];
          const size_t vdim = refGrad.empty() ? 0 : refGrad.front().size();

          for (size_t i = 0; i < refGrad.size(); ++i)
          {
            ScalarType div = 0;
            for (size_t comp = 0; comp < std::min(vdim, d); ++comp)
            {
              ScalarType physComp = 0;
              for (size_t j = 0; j < d; ++j)
                physComp += Jinv(comp, j) * refGrad[i][comp](j);
              div += physComp;
            }

            const size_t ntr = trialBasis.size();
            for (size_t tr = 0; tr < ntr; ++tr)
              m_matrix(i, tr) += wdet * trialBasis[tr] * div;
          }
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_matrix(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      const Geometry::PolytopeQuadrature* m_quadrature;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      /// @brief Reference gradients per quadrature point, dof and component.
      std::vector<std::vector<std::vector<Math::SpatialVector<ScalarType>>>> m_refGrad;
      /// @brief Trial basis values per quadrature point and dof.
      std::vector<std::vector<ScalarType>> m_trialBasis;
      Math::Matrix<ScalarType> m_matrix;
  };

  template <class LHSDerived, class RHSDerived, class LHSMesh, class RHSMesh>
  using P1DivTestIntegrand =
    Dot<
      ShapeFunctionBase<
        ShapeFunction<LHSDerived, P1<Real, LHSMesh>, TrialSpace>,
        P1<Real, LHSMesh>, TrialSpace>,
      ShapeFunctionBase<
        Div<ShapeFunction<RHSDerived, P1<Math::SpatialVector<Real>, RHSMesh>, TestSpace>>,
        P1<Math::SpatialVector<Real>, RHSMesh>, TestSpace>>;

  template <class LHSDerived, class RHSDerived, class LHSMesh, class RHSMesh>
  QuadratureRule(const P1DivTestIntegrand<LHSDerived, RHSDerived, LHSMesh, RHSMesh>&)
    -> QuadratureRule<P1DivTestIntegrand<LHSDerived, RHSDerived, LHSMesh, RHSMesh>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Integration of the isotropic Frobenius inner product two instances
   * of the P1 Jacobian of ShapeFunction.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int \mathbf{J} \: u : \mathbf{J} \: v \ dx \: ,
   * @f]
   * where @f$ u \in \mathbb{P}_1 @f$ and @f$ v \in
   * \mathbb{P}_1 @f$, and @f$ A @f$ is coefficient function.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int \mathbf{J} \: u : \mathbf{J} \: v \ dx :
   * \texttt{QuadratureRule}}
   * {\vdash u, v : \mathbb{P}_1}
   * @f]
   */
  template <class LHSDerived, class RHSDerived, class LHSRange, class RHSRange, class LHSMesh, class RHSMesh>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Jacobian<ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>>,
          P1<LHSRange, LHSMesh>, TrialSpace>,
      ShapeFunctionBase<
        Jacobian<ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>>,
          P1<RHSRange, RHSMesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<
            ShapeFunctionBase<
              Jacobian<ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>>,
                P1<LHSRange, LHSMesh>, TrialSpace>,
            ShapeFunctionBase<
              Jacobian<ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>>,
                P1<RHSRange, RHSMesh>, TestSpace>>>
        ::ScalarType>
  {
    public:
      /// @brief Reports this handler as an optimized specialization.
      static constexpr bool Specialized = true;

      using LHSFESType = P1<LHSRange, LHSMesh>;

      using RHSFESType = P1<RHSRange, RHSMesh>;

      /// @brief Left-hand side operand type.
      using LHSType =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<LHSDerived, LHSFESType, TrialSpace>>>;

      using LHSOperandType =
        ShapeFunction<LHSDerived, LHSFESType, TrialSpace>;

      using LHSOperandRangeType =
        typename FormLanguage::Traits<LHSOperandType>::RangeType;

      /// @brief Right-hand side operand type.
      using RHSType =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<RHSDerived, RHSFESType, TestSpace>>>;

      using RHSOperandType =
        ShapeFunction<RHSDerived, RHSFESType, TestSpace>;

      using RHSOperandRangeType =
        typename FormLanguage::Traits<RHSOperandType>::RangeType;

      /// @brief Integrand expression type.
      using IntegrandType = Dot<LHSType, RHSType>;

      using IntegrandRangeType =
        typename FormLanguage::Traits<IntegrandType>::RangeType;

      /// @brief Scalar value type.
      using ScalarType =
        typename FormLanguage::Traits<IntegrandType>::ScalarType;

      /// @brief Parent class type.
      using Parent =
        LocalBilinearFormIntegratorBase<ScalarType>;

      static_assert(FormLanguage::IsVectorRange<LHSOperandRangeType>::Value);
      static_assert(FormLanguage::IsVectorRange<RHSOperandRangeType>::Value);

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_quadrature(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_quadrature(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::exchange(other.m_quadrature, nullptr)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_trialRefJac(std::move(other.m_trialRefJac)),
          m_testRefJac(std::move(other.m_testRefJac)),
          m_matrix(std::move(other.m_matrix))
      {}

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes.getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe.getOrder();
        const size_t order =
          this->getOrder(polytope).value_or(getIntegrand().getOrder(polytope).value_or(
            ((k_tr == 0 || k_te == 0) ? 0 : (k_tr + k_te - 2))));

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);

          // Cached per quadrature point, not at the first one: reference
          // derivatives are constant only on a simplex, so freezing them at
          // one point is exact only there, or for a single-point rule.
          const size_t nqp = m_qf->getSize();
          m_trialRefJac.assign(nqp, {});
          m_testRefJac.assign(nqp, {});
          for (size_t qp = 0; qp < nqp; ++qp)
          {
            const auto& rc = m_qf->getPoint(qp);

            m_trialRefJac[qp].resize(trialfe.getCount());
            for (size_t local = 0; local < trialfe.getCount(); ++local)
            {
              auto& J = m_trialRefJac[qp][local];
              J.resize(trialfes.getVectorDimension(), d);
              const auto& basis = trialfe.getBasis(local);
              for (size_t i = 0; i < trialfes.getVectorDimension(); ++i)
                for (size_t j = 0; j < d; ++j)
                  J(i, j) = basis.template getDerivative<1>(i, j)(rc);
            }

            if (trialfes == testfes)
            {
              m_testRefJac[qp] = m_trialRefJac[qp];
            }
            else
            {
              m_testRefJac[qp].resize(testfe.getCount());
              for (size_t local = 0; local < testfe.getCount(); ++local)
              {
                auto& J = m_testRefJac[qp][local];
                J.resize(testfes.getVectorDimension(), d);
                const auto& basis = testfe.getBasis(local);
                for (size_t i = 0; i < testfes.getVectorDimension(); ++i)
                  for (size_t j = 0; j < d; ++j)
                    J(i, j) = basis.template getDerivative<1>(i, j)(rc);
              }
            }
          }

          m_matrix.resize(testfe.getCount(), trialfe.getCount());
        }

        assert(m_qf);
        m_quadrature = &polytope.getQuadrature(*m_qf);

        m_matrix.setZero();

        assert(m_quadrature);
        const auto& q = *m_quadrature;
        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto& Jinv = p.getJacobianInverse();
          const auto& trialRefJac = m_trialRefJac[qp];
          const auto& testRefJac = m_testRefJac[qp];

          if (trialfes == testfes)
          {
            const size_t n = trialRefJac.size();
            for (size_t i = 0; i < n; ++i)
            {
              const auto Ji = trialRefJac[i] * Jinv;
              m_matrix(i, i) += wdet * Ji.squaredNorm();

              for (size_t j = 0; j < i; ++j)
                m_matrix(i, j) += wdet * Math::dot(trialRefJac[j] * Jinv, Ji);
            }
          }
          else
          {
            const size_t ntr = trialRefJac.size();
            const size_t nte = testRefJac.size();

            for (size_t te = 0; te < nte; ++te)
            {
              const auto Jte = testRefJac[te] * Jinv;
              for (size_t tr = 0; tr < ntr; ++tr)
                m_matrix(te, tr) += wdet * Math::dot(trialRefJac[tr] * Jinv, Jte);
            }
          }
        }

        if (trialfes == testfes)
        {
          m_matrix.template triangularView<Eigen::Upper>() =
            m_matrix.template triangularView<Eigen::Lower>().transpose();
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_matrix(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      const Geometry::PolytopeQuadrature* m_quadrature;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      /// @brief Reference Jacobians per quadrature point and dof.
      std::vector<std::vector<Math::SpatialMatrix<ScalarType>>> m_trialRefJac,
        m_testRefJac;

      Math::Matrix<ScalarType> m_matrix;
  };

  /**
   * @ingroup RodinCTAD
   */
  template <class LHSDerived, class RHSDerived, class Range, class Mesh>
  QuadratureRule(const
    Dot<
      ShapeFunctionBase<
        Jacobian<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>>,
      ShapeFunctionBase<
        Jacobian<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>&)
    ->
      QuadratureRule<
        Dot<
          ShapeFunctionBase<
            Jacobian<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>>,
          ShapeFunctionBase<
            Jacobian<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Integration of the anisotropic Frobenius inner product two
   * instances of the P1 Jacobian of ShapeFunction.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int (A \: \mathbf{J} \: u) : \mathbf{J} \: v \ dx \: ,
   * @f]
   * where @f$ u \in \mathbb{P}_1 @f$ and @f$ v \in
   * \mathbb{P}_1 @f$, and @f$ A @f$ is a coefficient function.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int (A \: \mathbf{J} \: u) : \mathbf{J} \: v \ dx :
   * \texttt{QuadratureRule}}
   * {\vdash u, v : \mathbb{P}_1}
   * @f]
   *
   * Uses a single centroid quadrature point and supports different P1 trial
   * and test spaces.
   */
  template <
    class CoefficientDerived, class LHSDerived, class RHSDerived,
    class LHSRange, class RHSRange, class LHSMesh, class RHSMesh>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Mult<
          FunctionBase<CoefficientDerived>,
          ShapeFunctionBase<
            Jacobian<ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>>,
            P1<LHSRange, LHSMesh>, TrialSpace>>,
        P1<LHSRange, LHSMesh>, TrialSpace>,
      ShapeFunctionBase<
        Jacobian<ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>>,
        P1<RHSRange, RHSMesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<
            ShapeFunctionBase<
              Mult<
                FunctionBase<CoefficientDerived>,
                ShapeFunctionBase<
                  Jacobian<ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>>,
                  P1<LHSRange, LHSMesh>, TrialSpace>>,
              P1<LHSRange, LHSMesh>, TrialSpace>,
            ShapeFunctionBase<
              Jacobian<ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>>,
              P1<RHSRange, RHSMesh>, TestSpace>>>
        ::ScalarType>
  {
    public:
      /// @brief Reports this handler as an optimized specialization.
      static constexpr bool Specialized = true;

      using LHSFESType = P1<LHSRange, LHSMesh>;

      using RHSFESType = P1<RHSRange, RHSMesh>;

      using CoefficientType = FunctionBase<CoefficientDerived>;

      using CoefficientRangeType =
        typename FormLanguage::Traits<CoefficientType>::RangeType;

      using MultiplicandType =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<LHSDerived, LHSFESType, TrialSpace>>>;

      /// @brief Left-hand side operand type.
      using LHSType =
        ShapeFunctionBase<Mult<CoefficientType, MultiplicandType>>;

      using LHSOperandType =
        ShapeFunction<LHSDerived, LHSFESType, TrialSpace>;

      using LHSOperandRangeType =
        typename FormLanguage::Traits<LHSOperandType>::RangeType;

      /// @brief Right-hand side operand type.
      using RHSType =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>>>;

      using RHSOperandType =
        ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>;

      using RHSOperandRangeType =
        typename FormLanguage::Traits<RHSOperandType>::RangeType;

      /// @brief Integrand expression type.
      using IntegrandType = Dot<LHSType, RHSType>;

      using IntegrandRangeType =
        typename FormLanguage::Traits<IntegrandType>::RangeType;

      /// @brief Scalar value type.
      using ScalarType =
        typename FormLanguage::Traits<IntegrandType>::ScalarType;

      /// @brief Parent class type.
      using Parent =
        LocalBilinearFormIntegratorBase<ScalarType>;

      static_assert(FormLanguage::IsVectorRange<LHSOperandRangeType>::Value);
      static_assert(FormLanguage::IsVectorRange<RHSOperandRangeType>::Value);

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_quadrature(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_quadrature(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::exchange(other.m_quadrature, nullptr)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_trialRefJac(std::move(other.m_trialRefJac)),
          m_testRefJac(std::move(other.m_testRefJac)),
          m_matrix(std::move(other.m_matrix))
      {}

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& coeff = lhs.getDerived().getLHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes.getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe.getOrder();
        const size_t order =
          this->getOrder(polytope).value_or(getIntegrand().getOrder(polytope).value_or(
            ((k_tr == 0 || k_te == 0) ? 0 : (k_tr + k_te - 2))));

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);

          // Cached per quadrature point, not at the first one: reference
          // derivatives are constant only on a simplex, so freezing them at
          // one point is exact only there, or for a single-point rule.
          const size_t nqp = m_qf->getSize();
          m_trialRefJac.assign(nqp, {});
          m_testRefJac.assign(nqp, {});
          for (size_t qp = 0; qp < nqp; ++qp)
          {
            const auto& rc = m_qf->getPoint(qp);

            m_trialRefJac[qp].resize(trialfe.getCount());
            for (size_t local = 0; local < trialfe.getCount(); ++local)
            {
              auto& J = m_trialRefJac[qp][local];
              J.resize(trialfes.getVectorDimension(), d);
              const auto& basis = trialfe.getBasis(local);
              for (size_t i = 0; i < trialfes.getVectorDimension(); ++i)
                for (size_t j = 0; j < d; ++j)
                  J(i, j) = basis.template getDerivative<1>(i, j)(rc);
            }

            if (trialfes == testfes)
            {
              m_testRefJac[qp] = m_trialRefJac[qp];
            }
            else
            {
              m_testRefJac[qp].resize(testfe.getCount());
              for (size_t local = 0; local < testfe.getCount(); ++local)
              {
                auto& J = m_testRefJac[qp][local];
                J.resize(testfes.getVectorDimension(), d);
                const auto& basis = testfe.getBasis(local);
                for (size_t i = 0; i < testfes.getVectorDimension(); ++i)
                  for (size_t j = 0; j < d; ++j)
                    J(i, j) = basis.template getDerivative<1>(i, j)(rc);
              }
            }
          }

          m_matrix.resize(testfe.getCount(), trialfe.getCount());
        }

        assert(m_qf);
        m_quadrature = &polytope.getQuadrature(*m_qf);

        m_matrix.setZero();

        assert(m_quadrature);
        const auto& q = *m_quadrature;
        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const IntegrationPoint ip(p, m_qf, qp);
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto& Jinv = p.getJacobianInverse();
          const auto& trialRefJac = m_trialRefJac[qp];
          const auto& testRefJac = m_testRefJac[qp];

          if constexpr (std::is_same_v<CoefficientRangeType, ScalarType>)
          {
            const ScalarType csv = coeff.getValue(ip);

            if (trialfes == testfes)
            {
              const size_t n = trialRefJac.size();
              for (size_t i = 0; i < n; ++i)
              {
                const auto Ji = trialRefJac[i] * Jinv;
                m_matrix(i, i) += wdet * csv * Ji.squaredNorm();

                for (size_t j = 0; j < i; ++j)
                  m_matrix(i, j) += wdet * csv * Math::dot(trialRefJac[j] * Jinv, Ji);
              }
            }
            else
            {
              const size_t ntr = trialRefJac.size();
              const size_t nte = testRefJac.size();

              for (size_t te = 0; te < nte; ++te)
              {
                const auto Jte = testRefJac[te] * Jinv;
                for (size_t tr = 0; tr < ntr; ++tr)
                  m_matrix(te, tr) += wdet * csv * Math::dot(trialRefJac[tr] * Jinv, Jte);
              }
            }
          }
          else if constexpr (FormLanguage::IsMatrixRange<CoefficientRangeType>::Value)
          {
            Math::SpatialMatrix<ScalarType> cmv;
            cmv = coeff.getValue(ip);

            if (trialfes == testfes)
            {
              const size_t n = trialRefJac.size();
              for (size_t i = 0; i < n; ++i)
              {
                const auto Ji = trialRefJac[i] * Jinv;
                m_matrix(i, i) += wdet * Math::dot(cmv * Ji, Ji);

                for (size_t j = 0; j < i; ++j)
                  m_matrix(i, j) += wdet * Math::dot(cmv * (trialRefJac[j] * Jinv), Ji);
              }

              for (size_t i = 0; i < n; ++i)
                for (size_t j = i + 1; j < n; ++j)
                  m_matrix(i, j) += wdet *
                    Math::dot(cmv * (trialRefJac[j] * Jinv), trialRefJac[i] * Jinv);
            }
            else
            {
              const size_t ntr = trialRefJac.size();
              const size_t nte = testRefJac.size();

              for (size_t te = 0; te < nte; ++te)
              {
                const auto Jte = testRefJac[te] * Jinv;
                for (size_t tr = 0; tr < ntr; ++tr)
                  m_matrix(te, tr) +=
                    wdet * Math::dot(cmv * (trialRefJac[tr] * Jinv), Jte);
              }
            }
          }
          else
          {
            assert(false);
          }
        }

        if constexpr (std::is_same_v<CoefficientRangeType, ScalarType>)
        {
          if (trialfes == testfes)
          {
            m_matrix.template triangularView<Eigen::Upper>() =
              m_matrix.template triangularView<Eigen::Lower>().transpose();
          }
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_matrix(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      const Geometry::PolytopeQuadrature* m_quadrature;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      /// @brief Reference Jacobians per quadrature point and dof.
      std::vector<std::vector<Math::SpatialMatrix<ScalarType>>> m_trialRefJac,
        m_testRefJac;

      Math::Matrix<ScalarType> m_matrix;
  };

  /**
   * @ingroup RodinCTAD
   */
  template <class LHSFunctionDerived, class LHSDerived, class RHSDerived, class Mesh>
  QuadratureRule(const
    Dot<
      ShapeFunctionBase<
        Mult<
          FunctionBase<LHSFunctionDerived>,
          ShapeFunctionBase<Jacobian<ShapeFunction<LHSDerived, P1<Math::SpatialVector<Real>, Mesh>, TrialSpace>>>>>,
      ShapeFunctionBase<
        Jacobian<ShapeFunction<RHSDerived, P1<Math::SpatialVector<Real>, Mesh>, TestSpace>>>>&)
  ->
  QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Mult<
          FunctionBase<LHSFunctionDerived>,
          ShapeFunctionBase<
            Jacobian<ShapeFunction<LHSDerived, P1<Math::SpatialVector<Real>, Mesh>, TrialSpace>>>>>,
      ShapeFunctionBase<
        Jacobian<ShapeFunction<RHSDerived, P1<Math::SpatialVector<Real>, Mesh>, TestSpace>>>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Specialization for @f$ \int (\nabla u \cdot f) \cdot v @f$ in the case of P1 shape functions
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int (\mathbf{J} \: u \cdot f) \cdot v \ dx \: ,
   * @f]
   * where @f$ u \in \mathbb{P}_1 @f$ and @f$ v \in \mathbb{P}_1 @f$, and
   * @f$ f @f$ is a coefficient function.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int (\mathbf{J} \: u \cdot f) \cdot v \ dx :
   * \texttt{QuadratureRule}}
   * {\vdash u, v : \mathbb{P}_1}
   * @f]
   *
   * Evaluates at the element centroid. Supports distinct P1 trial and test
   * spaces. Implements the linearized convection term for Navier-Stokes.
   */
  template <
    class CoefficientDerived, class LHSDerived, class RHSDerived,
    class LHSRange, class RHSRange, class LHSMesh, class RHSMesh>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Mult<
          ShapeFunctionBase<
            Jacobian<ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>>,
            P1<LHSRange, LHSMesh>, TrialSpace>,
          FunctionBase<CoefficientDerived>>,
        P1<LHSRange, LHSMesh>, TrialSpace>,
      ShapeFunctionBase<
        ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>,
        P1<RHSRange, RHSMesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<
            ShapeFunctionBase<
              Mult<
                ShapeFunctionBase<
                  Jacobian<ShapeFunction<LHSDerived, P1<LHSRange, LHSMesh>, TrialSpace>>,
                  P1<LHSRange, LHSMesh>, TrialSpace>,
                FunctionBase<CoefficientDerived>>,
              P1<LHSRange, LHSMesh>, TrialSpace>,
            ShapeFunctionBase<
              ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>,
              P1<RHSRange, RHSMesh>, TestSpace>>>::ScalarType>
  {
    public:
      /// @brief Reports this handler as an optimized specialization.
      static constexpr bool Specialized = true;

      using TrialFESType = P1<LHSRange, LHSMesh>;
      using TestFESType  = P1<RHSRange, RHSMesh>;

      using TrialSFType =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>>,
          TrialFESType, TrialSpace>;

      using CoefficientType = FunctionBase<CoefficientDerived>;

      /// @brief Left-hand side operand type.
      using LHSType =
        ShapeFunctionBase<
          Mult<TrialSFType, CoefficientType>,
          TrialFESType, TrialSpace>;

      /// @brief Right-hand side operand type.
      using RHSType =
        ShapeFunctionBase<
          ShapeFunction<RHSDerived, TestFESType, TestSpace>,
          TestFESType, TestSpace>;

      /// @brief Integrand expression type.
      using IntegrandType = Dot<LHSType, RHSType>;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;
      /// @brief Parent class type.
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_quadrature(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_quadrature(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::exchange(other.m_quadrature, nullptr)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_refGrad(std::move(other.m_refGrad)),
          m_basis(std::move(other.m_basis)),
          m_matrix(std::move(other.m_matrix))
      {}

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = rhs.getFiniteElementSpace().getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe.getOrder();
        const size_t order =
          this->getOrder(polytope).value_or(getIntegrand().getOrder(polytope).value_or(
            ((k_tr == 0 || k_te == 0) ? 0 : (k_tr + k_te - 2))));

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        // The coefficient is the RHS of the Mult node
        const auto& coeff = lhs.getDerived().getRHS();

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);

          // Cached per quadrature point, not at the first one: basis values
          // vary within every element and reference gradients vary within
          // every element that is not a simplex, so freezing either at one
          // point is exact only for a single-point rule.
          const size_t nqp = m_qf->getSize();
          m_refGrad.assign(nqp, {});
          m_basis.assign(nqp, {});
          for (size_t qp = 0; qp < nqp; ++qp)
          {
            const auto& rc = m_qf->getPoint(qp);

            if constexpr (std::is_same_v<LHSRange, ScalarType>)
            {
              const size_t n = trialfe.getCount();
              m_refGrad[qp].resize(n);
              for (size_t a = 0; a < n; ++a)
              {
                m_refGrad[qp][a].resize(d);
                const auto& basisFn = trialfe.getBasis(a);
                for (size_t j = 0; j < d; ++j)
                  m_refGrad[qp][a](j) = basisFn.template getDerivative<1>(j)(rc);
              }

              m_basis[qp].resize(n);
              for (size_t b = 0; b < n; ++b)
                m_basis[qp][b] = testfe.getBasis(b)(rc);
            }
            else
            {
              const size_t vdim = trialfes.getVectorDimension();
              const size_t nVertices = trialfe.getCount() / vdim;

              m_refGrad[qp].resize(nVertices);
              for (size_t v = 0; v < nVertices; ++v)
              {
                m_refGrad[qp][v].resize(d);
                const auto& basisFn = trialfe.getBasis(v * vdim);
                for (size_t j = 0; j < d; ++j)
                  m_refGrad[qp][v](j) = basisFn.template getDerivative<1>(0, j)(rc);
              }

              m_basis[qp].resize(nVertices);
              for (size_t v = 0; v < nVertices; ++v)
              {
                const auto bv = testfe.getBasis(v * vdim)(rc);
                m_basis[qp][v] = bv(0);
              }
            }
          }
        }

        assert(m_qf);
        m_quadrature = &polytope.getQuadrature(*m_qf);

        const size_t n = m_refGrad.empty() ? 0 : m_refGrad.front().size();
        const size_t vdim = trialfes.getVectorDimension();

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);
        m_matrix.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_matrix.setZero();

        assert(m_quadrature);
        const auto& q = *m_quadrature;
        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const IntegrationPoint ip(p, m_qf, qp);
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto& Jinv = p.getJacobianInverse();
          const auto fval = coeff.getValue(ip);
          const auto& refGrad = m_refGrad[qp];
          const auto& basis = m_basis[qp];

          for (size_t a = 0; a < n; ++a)
          {
            const auto physGrad = Jinv.transpose() * refGrad[a];
            const ScalarType gradDotF = Math::dot(physGrad, fval);

            for (size_t b = 0; b < n; ++b)
            {
              const ScalarType val = wdet * gradDotF * basis[b];

              for (size_t c = 0; c < vdim; ++c)
              {
                const size_t row = b * vdim + c;
                const size_t col = a * vdim + c;
                m_matrix(
                  static_cast<Eigen::Index>(row),
                  static_cast<Eigen::Index>(col)) += val;
              }
            }
          }
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_matrix(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      const Geometry::PolytopeQuadrature* m_quadrature;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      /// @brief Reference gradients per quadrature point and dof.
      std::vector<std::vector<Math::SpatialVector<ScalarType>>> m_refGrad;
      /// @brief Test basis values per quadrature point and dof.
      std::vector<std::vector<ScalarType>> m_basis;

      Math::Matrix<ScalarType> m_matrix;
  };

  /**
   * @ingroup RodinCTAD
   */
  template <
    class CoefficientDerived, class LHSDerived, class RHSDerived,
    class Range, class Mesh>
  QuadratureRule(
    const Dot<
      ShapeFunctionBase<
        Mult<
          ShapeFunctionBase<
            Jacobian<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>,
            P1<Range, Mesh>, TrialSpace>,
          FunctionBase<CoefficientDerived>>,
        P1<Range, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>,
        P1<Range, Mesh>, TestSpace>>&)
    -> QuadratureRule<
      Dot<
        ShapeFunctionBase<
          Mult<
            ShapeFunctionBase<
              Jacobian<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>,
              P1<Range, Mesh>, TrialSpace>,
            FunctionBase<CoefficientDerived>>,
          P1<Range, Mesh>, TrialSpace>,
        ShapeFunctionBase<
          ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>,
          P1<Range, Mesh>, TestSpace>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Specialization for @f$\int \mathcal{K}(u) \cdot v@f$ in the case of P1 trial/test shape functions.
   *
   * This class represents the CTAD for the expression:
   * @f[
   * \int \mathcal{K}(u) \cdot v \: ,
   * @f]
   * where @f$ u, v \in \mathbb{P}_1 @f$ and @f$\mathcal{K}@f$ is a potential kernel operator.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a well formed type
   * of QuadratureRule.
   * @f[
   * \dfrac
   * {\vdash \int \mathcal{K}(u) \cdot v : \texttt{QuadratureRule}}
   * {\vdash u, v : \mathbb{P}_1}
   * @f]
   *
   * When the trial and test polytopes coincide, triangles use the optimized
   * six-point construction; all other geometries fall back to a single centroid
   * pairing consistent with the generic QuadratureRule.
   */
  template <class Kernel, class Range, class Mesh, class LHSDerived, class RHSDerived>
  class QuadratureRule<
    Dot<
      Potential<Kernel, ShapeFunctionBase<LHSDerived, P1<Range, Mesh>, TrialSpace>>,
      ShapeFunctionBase<RHSDerived, P1<Range, Mesh>, TestSpace>>>
    : public GlobalBilinearFormIntegratorBase<typename FormLanguage::Traits<Range>::ScalarType>
  {
    public:
      /// @brief Reports this handler as an optimized specialization.
      static constexpr bool Specialized = true;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<Range>::ScalarType;

      using KernelType = Kernel;

      using TrialFESType = P1<Range, Mesh>;

      using TestFESType = P1<Range, Mesh>;

      /// @brief Left-hand side operand type.
      using LHSType =
        Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFESType, TrialSpace>>;

      /// @brief Right-hand side operand type.
      using RHSType =
        ShapeFunctionBase<RHSDerived, TestFESType, TestSpace>;

      /// @brief Range type of the left-hand side operand.
      using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;

      /// @brief Range type of the right-hand side operand.
      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      /// @brief Integrand expression type.
      using IntegrandType = Dot<LHSType, RHSType>;

      /// @brief Parent class type.
      using Parent = GlobalBilinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);

      constexpr
      QuadratureRule(const LHSType& lhs, const RHSType& rhs)
        : QuadratureRule(Dot(lhs, rhs))
      {}

      constexpr
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getOperand().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qfs(Geometry::Polytope::Type::Segment)
      {}

      constexpr
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qfs(other.m_qfs)
      {}

      constexpr
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_trp(std::move(other.m_trp)),
          m_tep(std::move(other.m_tep)),
          m_qfs(std::move(other.m_qfs)),
          m_qftr(std::move(other.m_qftr)),
          m_qfte(std::move(other.m_qfte)),
          m_weight(std::move(other.m_weight)),
          m_distortion(std::move(other.m_distortion)),
          m_sk(std::move(other.m_sk)),
          m_mk(std::move(other.m_mk)),
          m_trv(std::move(other.m_trv)),
          m_tev(std::move(other.m_tev)),
          m_k0(std::move(other.m_k0)),
          m_k1(std::move(other.m_k1)),
          m_k2(std::move(other.m_k2)),
          m_k3(std::move(other.m_k3)),
          m_k4(std::move(other.m_k4)),
          m_k5(std::move(other.m_k5)),
          m_matrix(std::move(other.m_matrix))
      {}

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      QuadratureRule& setPolytope(
        const Geometry::Polytope& trp, const Geometry::Polytope& tep) override
      {
        m_trp = trp;
        m_tep = tep;

        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& kernel = lhs.getKernel();
        const auto& trialfes = lhs.getOperand().getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        const auto& trg = trp.getGeometry();
        const auto& teg = tep.getGeometry();

        const auto& trialfe = trialfes.getFiniteElement(
          trp.getDimension(), trp.getIndex());
        const auto& testfe = testfes.getFiniteElement(
          tep.getDimension(), tep.getIndex());

        // The regularised integrand is smooth but not polynomial for the
        // singular kernels this exists for, so the order is an accuracy
        // choice rather than a degree that can be inferred. Seven is exact
        // for the Jacobian and the first-order bases over it, and is the
        // default when the caller states nothing.
        {
          const size_t order = this->getOrder(trp).value_or(7);
          const size_t n = std::max<size_t>(1, (order + 2) / 2);
          if (m_collapsed.size() != n)
          {
            std::vector<Real> nodes, weights;
            QF::GaussLegendre::gl1dUnit(n, nodes, weights);
            m_collapsed.clear();
            m_collapsed.reserve(n);
            for (size_t i = 0; i < n; ++i)
              m_collapsed.emplace_back(nodes[i], weights[i]);
          }
        }

        if constexpr (std::is_same_v<Range, ScalarType>)
        {

          if (trp == tep)
          {
            const auto& polytope = trp;
            switch (polytope.getGeometry())
            {
              case Geometry::Polytope::Type::Triangle:
              {
                Math::SpatialVector<ScalarType> rx0(2), rz0(2),
                                                rx1(2), rz1(2),
                                                rx2(2), rz2(2),
                                                rx3(2), rz3(2),
                                                rx4(2), rz4(2),
                                                rx5(2), rz5(2);

                // Sauter-Schwab for coincident panels: the four-dimensional
                // integral is regularised onto the unit cube in the collapsed
                // variables, and what remains is smooth and integrated by a
                // tensor Gauss rule --- the same rule in each direction, as
                // the construction intends.
                //
                // Evaluating it at one point instead, as this did, does not
                // integrate even the Jacobian: xi^3 eta1^2 eta2 has midpoint
                // value 1/64 against an exact 1/24 over the cube, so a
                // constant kernel came out at three eighths of its value on
                // every coincident pair. The transformation was never the
                // problem.
                m_matrix.setZero(testfe.getCount(), trialfe.getCount());
                for (size_t i3 = 0; i3 < m_collapsed.size(); ++i3)
                {
                  const ScalarType eta3 = m_collapsed[i3].first;
                  const Real w3 = m_collapsed[i3].second;
                  for (size_t i2 = 0; i2 < m_collapsed.size(); ++i2)
                  {
                    const ScalarType eta2 = m_collapsed[i2].first;
                    const Real w2 = m_collapsed[i2].second;
                    for (size_t i1 = 0; i1 < m_collapsed.size(); ++i1)
                    {
                      const ScalarType eta1 = m_collapsed[i1].first;
                      const Real w1 = m_collapsed[i1].second;
                      for (size_t i0 = 0; i0 < m_collapsed.size(); ++i0)
                      {
                        const ScalarType xi = m_collapsed[i0].first;
                        const Real w0 = m_collapsed[i0].second;
                        const Real jacobian = xi * xi * xi * eta1 * eta1 * eta2;
                        const Real factor = w0 * w1 * w2 * w3 * jacobian;

                        // The six sub-domains of the identical-panel case,
                        // in the variables of Sauter and Schwab. Written as
                        // they appear there: the first coordinate carries
                        // 1 - xi and the second carries xi, which a symmetric
                        // one-point rule cannot tell apart and a real rule
                        // can. Getting these the wrong way round leaves the
                        // Jacobian integrating correctly, so the entries sum
                        // to the right total while being distributed wrongly
                        // among the basis functions.
                        const ScalarType a = xi - xi * eta1;
                        const ScalarType b = xi - xi * eta1 * eta2 * eta3;
                        const ScalarType c = xi - xi * eta1 * eta2;
                        const ScalarType d = xi * eta1 * (1 - eta2);
                        const ScalarType e = xi * eta1 * (1 - eta2 * eta3);
                        const ScalarType f = xi * eta1 * (1 - eta2 + eta2 * eta3);
                        const ScalarType g = xi - xi * eta1 + xi * eta1 * eta2;

                        rx0[0] = 1 - xi;
                        rx0[1] = g;
                        rz0[0] = 1 - b;
                        rz0[1] = a;

                        rx1[0] = 1 - b;
                        rx1[1] = a;
                        rz1[0] = 1 - xi;
                        rz1[1] = g;

                        rx2[0] = 1 - xi;
                        rx2[1] = f;
                        rz2[0] = 1 - c;
                        rz2[1] = d;

                        rx3[0] = 1 - c;
                        rx3[1] = d;
                        rz3[0] = 1 - xi;
                        rz3[1] = f;

                        rx4[0] = 1 - b;
                        rx4[1] = e;
                        rz4[0] = 1 - xi;
                        rz4[1] = d;

                        rx5[0] = 1 - xi;
                        rx5[1] = d;
                        rz5[0] = 1 - b;
                        rz5[1] = e;

                        const Geometry::Point x0(polytope, rx0);
                        const Geometry::Point x1(polytope, rx1);
                        const Geometry::Point x2(polytope, rx2);
                        const Geometry::Point x3(polytope, rx3);
                        const Geometry::Point x4(polytope, rx4);
                        const Geometry::Point x5(polytope, rx5);

                        const Geometry::Point z0(polytope, rz0);
                        const Geometry::Point z1(polytope, rz1);
                        const Geometry::Point z2(polytope, rz2);
                        const Geometry::Point z3(polytope, rz3);
                        const Geometry::Point z4(polytope, rz4);
                        const Geometry::Point z5(polytope, rz5);

                        const Real s0 = x0.getDistortion() * z0.getDistortion();
                        const Real s1 = x1.getDistortion() * z1.getDistortion();
                        const Real s2 = x2.getDistortion() * z2.getDistortion();
                        const Real s3 = x3.getDistortion() * z3.getDistortion();
                        const Real s4 = x4.getDistortion() * z4.getDistortion();
                        const Real s5 = x5.getDistortion() * z5.getDistortion();

                        assert(std::isfinite(s0));
                        assert(std::isfinite(s1));
                        assert(std::isfinite(s2));
                        assert(std::isfinite(s3));
                        assert(std::isfinite(s4));
                        assert(std::isfinite(s5));

                        for (size_t l = 0; l < testfe.getCount(); ++l)
                        {
                          const auto& teb = testfe.getBasis(l);
                          for (size_t m = 0; m < trialfe.getCount(); ++m)
                          {
                            const auto& trb = trialfe.getBasis(m);
                            m_matrix(l, m) +=
                              factor * s0 * kernel(x0, z0) * trb(rx0) * teb(rz0);
                            m_matrix(l, m) +=
                              factor * s1 * kernel(x1, z1) * trb(rx1) * teb(rz1);
                            m_matrix(l, m) +=
                              factor * s2 * kernel(x2, z2) * trb(rx2) * teb(rz2);
                            m_matrix(l, m) +=
                              factor * s3 * kernel(x3, z3) * trb(rx3) * teb(rz3);
                            m_matrix(l, m) +=
                              factor * s4 * kernel(x4, z4) * trb(rx4) * teb(rz4);
                            m_matrix(l, m) +=
                              factor * s5 * kernel(x5, z5) * trb(rx5) * teb(rz5);
                          }
                        }
                      }
                    }
                  }
                }
                m_distortion = 1;
                m_weight = 1;
                break;
              }
              default:
              {
                m_qftr.emplace(polytope.getGeometry());
                assert(m_qftr->getSize() == 1);
                m_qfte.emplace(polytope.getGeometry());
                assert(m_qfte->getSize() == 1);

                const auto& rx = m_qftr->getPoint(0);
                const auto& ry = m_qfte->getPoint(0);

                const Geometry::Point x(polytope, std::cref(rx));
                const Geometry::Point y(polytope, std::cref(ry));

                m_distortion = x.getDistortion() * y.getDistortion();
                m_weight = m_qftr->getWeight(0) * m_qfte->getWeight(0);
                m_matrix.resize(testfe.getCount(), trialfe.getCount());

                m_sk = kernel(x, y);
                for (size_t l = 0; l < testfe.getCount(); ++l)
                {
                  const ScalarType teb = testfe.getBasis(l)(ry);
                  for (size_t m = 0; m < trialfe.getCount(); ++m)
                  {
                    const ScalarType trb = trialfe.getBasis(m)(rx);
                    m_matrix(l, m) = m_sk * trb * teb;
                  }
                }
                break;
              }
            }
          }
          else
          {
            m_qftr.emplace(trg);
            assert(m_qftr->getSize() == 1);
            m_qfte.emplace(teg);
            assert(m_qfte->getSize() == 1);

            const auto& rx = m_qftr->getPoint(0);
            const auto& ry = m_qfte->getPoint(0);

            const Geometry::Point x(trp, std::cref(rx));
            const Geometry::Point y(tep, std::cref(ry));

            m_distortion = x.getDistortion() * y.getDistortion();
            m_weight = m_qftr->getWeight(0) * m_qfte->getWeight(0);
            m_matrix.resize(testfe.getCount(), trialfe.getCount());

            m_sk = kernel(x, y);
            for (size_t l = 0; l < testfe.getCount(); ++l)
            {
              const ScalarType teb = testfe.getBasis(l)(ry);
              for (size_t m = 0; m < trialfe.getCount(); ++m)
              {
                const ScalarType trb = trialfe.getBasis(m)(rx);
                m_matrix(l, m) = m_sk * trb * teb;
              }
            }
          }
        }
        else if constexpr (FormLanguage::IsVectorRange<Range>::Value)
        {

          if (trp == tep)
          {
            const auto& polytope = trp;
            switch (polytope.getGeometry())
            {
              case Geometry::Polytope::Type::Triangle:
              {
                Math::SpatialVector<ScalarType> rx0(2), rz0(2),
                                                rx1(2), rz1(2),
                                                rx2(2), rz2(2),
                                                rx3(2), rz3(2),
                                                rx4(2), rz4(2),
                                                rx5(2), rz5(2);

                // The vector-valued counterpart of the scalar case above, and
                // the same two corrections: a tensor Gauss rule over the
                // collapsed variables rather than a single point, and the
                // sub-domain maps written as Sauter and Schwab give them.
                m_matrix.setZero(testfe.getCount(), trialfe.getCount());
                for (size_t i3 = 0; i3 < m_collapsed.size(); ++i3)
                {
                  const ScalarType eta3 = m_collapsed[i3].first;
                  const Real w3 = m_collapsed[i3].second;
                  for (size_t i2 = 0; i2 < m_collapsed.size(); ++i2)
                  {
                    const ScalarType eta2 = m_collapsed[i2].first;
                    const Real w2 = m_collapsed[i2].second;
                    for (size_t i1 = 0; i1 < m_collapsed.size(); ++i1)
                    {
                      const ScalarType eta1 = m_collapsed[i1].first;
                      const Real w1 = m_collapsed[i1].second;
                      for (size_t i0 = 0; i0 < m_collapsed.size(); ++i0)
                      {
                        const ScalarType xi = m_collapsed[i0].first;
                        const Real w0 = m_collapsed[i0].second;
                        const Real jacobian = xi * xi * xi * eta1 * eta1 * eta2;
                        const Real factor = w0 * w1 * w2 * w3 * jacobian;

                        const ScalarType a = xi - xi * eta1;
                        const ScalarType b = xi - xi * eta1 * eta2 * eta3;
                        const ScalarType c = xi - xi * eta1 * eta2;
                        const ScalarType d = xi * eta1 * (1 - eta2);
                        const ScalarType e = xi * eta1 * (1 - eta2 * eta3);
                        const ScalarType f = xi * eta1 * (1 - eta2 + eta2 * eta3);
                        const ScalarType g = xi - xi * eta1 + xi * eta1 * eta2;

                        rx0[0] = 1 - xi;
                        rx0[1] = g;
                        rz0[0] = 1 - b;
                        rz0[1] = a;

                        rx1[0] = 1 - b;
                        rx1[1] = a;
                        rz1[0] = 1 - xi;
                        rz1[1] = g;

                        rx2[0] = 1 - xi;
                        rx2[1] = f;
                        rz2[0] = 1 - c;
                        rz2[1] = d;

                        rx3[0] = 1 - c;
                        rx3[1] = d;
                        rz3[0] = 1 - xi;
                        rz3[1] = f;

                        rx4[0] = 1 - b;
                        rx4[1] = e;
                        rz4[0] = 1 - xi;
                        rz4[1] = d;

                        rx5[0] = 1 - xi;
                        rx5[1] = d;
                        rz5[0] = 1 - b;
                        rz5[1] = e;

                        const Geometry::Point x0(polytope, rx0);
                        const Geometry::Point x1(polytope, rx1);
                        const Geometry::Point x2(polytope, rx2);
                        const Geometry::Point x3(polytope, rx3);
                        const Geometry::Point x4(polytope, rx4);
                        const Geometry::Point x5(polytope, rx5);

                        const Geometry::Point z0(polytope, rz0);
                        const Geometry::Point z1(polytope, rz1);
                        const Geometry::Point z2(polytope, rz2);
                        const Geometry::Point z3(polytope, rz3);
                        const Geometry::Point z4(polytope, rz4);
                        const Geometry::Point z5(polytope, rz5);

                        const Real s0 = x0.getDistortion() * z0.getDistortion();
                        const Real s1 = x1.getDistortion() * z1.getDistortion();
                        const Real s2 = x2.getDistortion() * z2.getDistortion();
                        const Real s3 = x3.getDistortion() * z3.getDistortion();
                        const Real s4 = x4.getDistortion() * z4.getDistortion();
                        const Real s5 = x5.getDistortion() * z5.getDistortion();

                        assert(std::isfinite(s0));
                        assert(std::isfinite(s1));
                        assert(std::isfinite(s2));
                        assert(std::isfinite(s3));
                        assert(std::isfinite(s4));
                        assert(std::isfinite(s5));

                        kernel(m_k0, x0, z0);
                        kernel(m_k1, x1, z1);
                        kernel(m_k2, x2, z2);
                        kernel(m_k3, x3, z3);
                        kernel(m_k4, x4, z4);
                        kernel(m_k5, x5, z5);

                        for (size_t l = 0; l < testfe.getCount(); ++l)
                        {
                          const auto& teb = testfe.getBasis(l);
                          for (size_t m = 0; m < trialfe.getCount(); ++m)
                          {
                            const auto& trb = trialfe.getBasis(m);
                            const auto trv0 = trb(rx0);
                            const auto tev0 = teb(rz0);
                            const auto trv1 = trb(rx1);
                            const auto tev1 = teb(rz1);
                            const auto trv2 = trb(rx2);
                            const auto tev2 = teb(rz2);
                            const auto trv3 = trb(rx3);
                            const auto tev3 = teb(rz3);
                            const auto trv4 = trb(rx4);
                            const auto tev4 = teb(rz4);
                            const auto trv5 = trb(rx5);
                            const auto tev5 = teb(rz5);

                            m_matrix(l, m) += factor * s0 * (m_k0 * trv0).dot(tev0);
                            m_matrix(l, m) += factor * s1 * (m_k1 * trv1).dot(tev1);
                            m_matrix(l, m) += factor * s2 * (m_k2 * trv2).dot(tev2);
                            m_matrix(l, m) += factor * s3 * (m_k3 * trv3).dot(tev3);
                            m_matrix(l, m) += factor * s4 * (m_k4 * trv4).dot(tev4);
                            m_matrix(l, m) += factor * s5 * (m_k5 * trv5).dot(tev5);
                          }
                        }
                      }
                    }
                  }
                }
                m_distortion = 1;
                m_weight = 1;
                break;
              }
              default:
              {
                m_qftr.emplace(polytope.getGeometry());
                assert(m_qftr->getSize() == 1);
                m_qfte.emplace(polytope.getGeometry());
                assert(m_qfte->getSize() == 1);

                const auto& rx = m_qftr->getPoint(0);
                const auto& ry = m_qfte->getPoint(0);

                const Geometry::Point x(polytope, std::cref(rx));
                const Geometry::Point y(polytope, std::cref(ry));

                m_distortion = x.getDistortion() * y.getDistortion();
                m_weight = m_qftr->getWeight(0) * m_qfte->getWeight(0);
                m_matrix.resize(testfe.getCount(), trialfe.getCount());

                kernel(m_mk, x, y);
                for (size_t l = 0; l < testfe.getCount(); ++l)
                {
                  m_tev = testfe.getBasis(l)(ry);
                  for (size_t m = 0; m < trialfe.getCount(); ++m)
                  {
                    m_trv = trialfe.getBasis(m)(rx);
                    m_matrix(l, m) = (m_mk * m_trv).dot(m_tev);
                  }
                }
                break;
              }
            }
          }
          else
          {
            m_qftr.emplace(trg);
            assert(m_qftr->getSize() == 1);
            m_qfte.emplace(teg);
            assert(m_qfte->getSize() == 1);

            const auto& rx = m_qftr->getPoint(0);
            const auto& ry = m_qfte->getPoint(0);

            const Geometry::Point x(trp, std::cref(rx));
            const Geometry::Point y(tep, std::cref(ry));

            m_distortion = x.getDistortion() * y.getDistortion();
            m_weight = m_qftr->getWeight(0) * m_qfte->getWeight(0);
            m_matrix.resize(testfe.getCount(), trialfe.getCount());

            kernel(m_mk, x, y);
            for (size_t l = 0; l < testfe.getCount(); ++l)
            {
              m_tev = testfe.getBasis(l)(ry);
              for (size_t m = 0; m < trialfe.getCount(); ++m)
              {
                m_trv = trialfe.getBasis(m)(rx);
                m_matrix(l, m) = (m_mk * m_trv).dot(m_tev);
              }
            }
          }
        }
        else
        {
          assert(false);
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) override
      {
        return m_distortion * m_weight * m_matrix(te, tr);
      }

      Geometry::Region getTrialRegion() const override
      {
        return getIntegrand().getLHS().getRegion();
      }

      virtual Geometry::Region getTestRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_trp;
      Optional<std::reference_wrapper<const Geometry::Polytope>> m_tep;

      /**
       * @brief Points and weights on @f$ [0, 1] @f$ for the collapsed
       * variables of a Sauter--Schwab transformation.
       *
       * One rule, applied in each of the four directions as a tensor product.
       * Its size follows from the integrator's order, so setOrder() controls
       * how finely the regularised integrand is resolved --- which is the only
       * knob that matters once the transformation has removed the singularity.
       */
      std::vector<std::pair<Real, Real>> m_collapsed;

      const QF::Centroid m_qfs;
      Optional<QF::Centroid> m_qftr;
      Optional<QF::Centroid> m_qfte;

      Real m_weight;
      Real m_distortion;

      ScalarType m_sk;
      Math::SpatialMatrix<ScalarType> m_mk;

      Math::SpatialVector<ScalarType> m_trv, m_tev;
      Math::SpatialMatrix<ScalarType> m_k0, m_k1, m_k2, m_k3, m_k4, m_k5;

      Math::Matrix<ScalarType> m_matrix;
  };
  /// @endcond
}

#endif
