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
 * @see P1, QuadratureFormula, Integral
 */
#ifndef RODIN_VARIATIONAL_P1_QUADRATURERULE_H
#define RODIN_VARIATIONAL_P1_QUADRATURERULE_H

#include "Rodin/FormLanguage/Traits.h"
#include "Rodin/Geometry/Region.h"
#include "Rodin/Math/Common.h"
#include "Rodin/Variational/ShapeFunction.h"
#include "Rodin/QF/Centroid.h"
#include "Rodin/QF/GenericPolytopeQuadrature.h"

#include "P1.h"
#include "P1Element.h"

namespace Rodin::Variational
{
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
      using FESType = P1<Range, Mesh>;

      using IntegrandType =
        ShapeFunctionBase<ShapeFunction<NestedDerived, FESType, TestSpace>>;

      using IntegrandRangeType = typename FormLanguage::Traits<FESType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent = LinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<IntegrandRangeType, ScalarType>);

      constexpr
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLeaf()),
          m_integrand(integrand.copy()),
          m_distortion(Math::nan<Real>()),
          m_weight(Math::nan<Real>()),
          m_set(false)
      {}

      constexpr
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_distortion(Math::nan<Real>()),
          m_weight(Math::nan<Real>()),
          m_set(false)
      {}

      constexpr
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_polytope(std::move(other.m_polytope)),
          m_qf(std::move(other.m_qf)),
          m_p(std::move(other.m_p)),
          m_distortion(std::move(other.m_distortion)),
          m_weight(std::move(other.m_weight)),
          m_basis(std::move(other.m_basis)),
          m_set(std::move(other.m_set)),
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
        const bool recompute = !m_set || m_geometry != geometry;

        if (recompute)
        {
          m_set = true;
          m_geometry = geometry;

          const P1Element<ScalarType> fe(geometry);
          const size_t n = fe.getCount();

          m_qf.emplace(geometry);
          assert(m_qf->getSize() == 1);

          m_p.emplace(polytope, m_qf->getPoint(0));
          m_weight = m_qf->getWeight(0);

          const auto& rc = m_qf->getPoint(0);

          m_basis.resize(n);
          for (size_t i = 0; i < n; ++i)
            m_basis[i] = fe.getBasis(i)(rc);
        }

        assert(m_p);
        auto& p = *m_p;
        p.setPolytope(polytope);
        m_distortion = p.getDistortion();

        return *this;
      }

      ScalarType integrate(size_t local) final override
      {
        return m_weight * m_distortion * m_basis[local];
      }

      virtual Geometry::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Optional<QF::Centroid> m_qf;
      Optional<Geometry::Point> m_p;

      Real m_distortion;
      Real m_weight;
      std::vector<ScalarType> m_basis;

      bool m_set;
      Optional<Geometry::Polytope::Type> m_geometry;
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
      using FESType =
        P1<Range, Mesh>;

      using LHSType =
        FunctionBase<LHSDerived>;

      using RHSType =
        ShapeFunctionBase<ShapeFunction<RHSDerived, FESType, TestSpace>, FESType, TestSpace>;

      using LHSRangeType =
        typename FormLanguage::Traits<LHSType>::RangeType;

      using RHSRangeType =
        typename FormLanguage::Traits<RHSType>::RangeType;

      using IntegrandType =
        ShapeFunctionBase<Dot<LHSType, RHSType>>;

      using IntegrandRangeType =
        typename FormLanguage::Traits<IntegrandType>::RangeType;

      using ScalarType =
        typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent =
        LinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);

      constexpr
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLeaf()),
          m_integrand(integrand.copy()),
          m_set(false)
      {}

      constexpr
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_set(false)
      {}

      constexpr
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_polytope(std::move(other.m_polytope)),
          m_qf(std::move(other.m_qf)),
          m_p(std::move(other.m_p)),
          m_distortion(std::move(other.m_distortion)),
          m_weight(std::move(other.m_weight)),
          m_basis(std::move(other.m_basis)),
          m_set(std::move(other.m_set)),
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
        static thread_local LHSRangeType s_v;
        m_polytope = polytope;
        const auto& geometry = polytope.getGeometry();
        const auto& integrand = getIntegrand();
        const auto& f = integrand.getDerived().getLHS();
        const auto& fes = integrand.getFiniteElementSpace();
        const bool recompute = !m_set || m_geometry != geometry;
        P1Element<RHSRangeType> fe;
        if constexpr (std::is_same_v<RHSRangeType, ScalarType>)
          fe = P1Element<RHSRangeType>(geometry);
        else if constexpr (std::is_same_v<RHSRangeType, Math::Vector<ScalarType>>)
          fe = P1Element<RHSRangeType>(geometry, fes.getVectorDimension());
        else
          assert(false);
        if (recompute)
        {
          m_set = true;
          m_geometry = geometry;
          m_qf.emplace(geometry);
          assert(m_qf->getSize() == 1);
          m_p.emplace(polytope, m_qf->getPoint(0));
          m_weight = m_qf->getWeight(0);
          m_basis.resize(fe.getCount());
          for (size_t local = 0; local < fe.getCount(); local++)
            m_basis[local] = fe.getBasis(local)(m_qf->getPoint(0));
          m_dot.resize(fe.getCount());
        }
        assert(m_p);
        auto& p = *m_p;
        p.setPolytope(polytope);
        s_v = f(p);
        for (size_t local = 0; local < fe.getCount(); local++)
          m_dot[local] = Math::dot(s_v, m_basis[local]);
        m_distortion = p.getDistortion();
        return *this;
      }

      ScalarType integrate(size_t local) final override
      {
        return m_weight * m_distortion * m_dot[local];
      }

      virtual Geometry::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:

      std::unique_ptr<IntegrandType> m_integrand;

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Optional<QF::Centroid> m_qf;
      Optional<Geometry::Point> m_p;

      Real m_weight;
      Real m_distortion;

      std::vector<RHSRangeType> m_basis;
      std::vector<ScalarType> m_dot;

      bool m_set;
      Optional<Geometry::Polytope::Type> m_geometry;
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
      using LHSFESType = P1<LHSRange, LHSMesh>;

      using RHSFESType = P1<RHSRange, RHSMesh>;

      using LHSType =
        ShapeFunctionBase<ShapeFunction<LHSDerived, LHSFESType, TrialSpace>>;

      using RHSType =
        ShapeFunctionBase<ShapeFunction<RHSDerived, RHSFESType, TestSpace>>;

      using IntegrandType = Dot<LHSType, RHSType>;

      using LHSRangeType =
        typename FormLanguage::Traits<LHSType>::RangeType;

      using RHSRangeType =
        typename FormLanguage::Traits<RHSType>::RangeType;

      using ScalarType =
        typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent =
        LocalBilinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);

      constexpr
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy())
      {}

      constexpr
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy())
      {}

      constexpr
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_polytope(std::move(other.m_polytope)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_ps(std::move(other.m_ps)),
          m_matrix(std::move(other.m_matrix)),
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

        // P1 × P1 product is degree 2 — need at least order-2 quadrature.
        const size_t order = trialfe.getOrder() + testfe.getOrder();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::GenericPolytopeQuadrature::get(order, geometry);

          m_ps.clear();
          m_ps.reserve(m_qf->getSize());
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps.emplace_back(polytope, m_qf->getPoint(qp));
        }
        else
        {
          assert(m_qf);
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps[qp].setPolytope(polytope);
        }

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);

        const P1Element<ScalarType> trialScalarFE(geometry);
        const P1Element<ScalarType> testScalarFE(geometry);
        const size_t scalarCountTr = trialScalarFE.getCount();
        const size_t scalarCountTe = testScalarFE.getCount();

        assert(scalarCountTr > 0 && ntr % scalarCountTr == 0);
        assert(scalarCountTe > 0 && nte % scalarCountTe == 0);
        const size_t vdim = ntr / scalarCountTr;
        assert(vdim == nte / scalarCountTe);

        const bool symmetric =
          (&trialfes.getMesh() == &testfes.getMesh()) && (ntr == nte);

        m_matrix.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_matrix.setZero();

        for (size_t qp = 0; qp < m_ps.size(); ++qp)
        {
          const auto& p = m_ps[qp];
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto& rc = m_qf->getPoint(qp);

          if (symmetric)
          {
            for (size_t ib = 0; ib < scalarCountTe; ++ib)
            {
              const ScalarType phi_te = testScalarFE.getBasis(ib)(rc);

              {
                const ScalarType kii = wdet * phi_te * trialScalarFE.getBasis(ib)(rc);
                for (size_t c = 0; c < vdim; ++c)
                  m_matrix(ib * vdim + c, ib * vdim + c) += kii;
              }

              for (size_t ia = 0; ia < ib; ++ia)
              {
                const ScalarType kij = wdet * phi_te * trialScalarFE.getBasis(ia)(rc);
                for (size_t c = 0; c < vdim; ++c)
                  m_matrix(ib * vdim + c, ia * vdim + c) += kij;
              }
            }
          }
          else
          {
            for (size_t ib = 0; ib < scalarCountTe; ++ib)
            {
              const ScalarType phi_te = testScalarFE.getBasis(ib)(rc);
              for (size_t ia = 0; ia < scalarCountTr; ++ia)
              {
                const ScalarType kij = wdet * phi_te * trialScalarFE.getBasis(ia)(rc);
                for (size_t c = 0; c < vdim; ++c)
                  m_matrix(ib * vdim + c, ia * vdim + c) += kij;
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
      std::vector<Geometry::Point> m_ps;

      Math::Matrix<ScalarType> m_matrix;

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
      using LHSFESType = P1<LHSRange, LHSMesh>;

      using RHSFESType = P1<RHSRange, RHSMesh>;

      using CoefficientType = FunctionBase<CoefficientDerived>;

      using MultiplicandType =
        ShapeFunctionBase<ShapeFunction<LHSDerived, LHSFESType, TrialSpace>>;

      using LHSType =
        ShapeFunctionBase<Mult<CoefficientType, MultiplicandType>>;

      using RHSType =
        ShapeFunctionBase<ShapeFunction<RHSDerived, RHSFESType, TestSpace>>;

      using IntegrandType = Dot<LHSType, RHSType>;

      using CoefficientRangeType =
        typename FormLanguage::Traits<CoefficientType>::RangeType;

      using MultiplicandRangeType =
        typename FormLanguage::Traits<MultiplicandType>::RangeType;

      using LHSRangeType =
        typename FormLanguage::Traits<LHSType>::RangeType;

      using RHSRangeType =
        typename FormLanguage::Traits<RHSType>::RangeType;

      using ScalarType =
        typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent =
        LocalBilinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);

      constexpr
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_distortion(Math::nan<Real>()),
          m_weight(Math::nan<Real>()),
          m_set(false)
      {}

      constexpr
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_distortion(Math::nan<Real>()),
          m_weight(Math::nan<Real>()),
          m_set(false)
      {}

      constexpr
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_polytope(std::move(other.m_polytope)),
          m_qf(std::move(other.m_qf)),
          m_p(std::move(other.m_p)),
          m_distortion(std::move(other.m_distortion)),
          m_weight(std::move(other.m_weight)),
          m_cmv(std::move(other.m_cmv)),
          m_trialBasis(std::move(other.m_trialBasis)),
          m_testBasis(std::move(other.m_testBasis)),
          m_matrix(std::move(other.m_matrix)),
          m_set(std::move(other.m_set)),
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
        const bool recompute = !m_set || m_geometry != geometry;

        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& coeff = lhs.getDerived().getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        if (recompute)
        {
          m_set = true;
          m_geometry = geometry;

          m_qf.emplace(geometry);
          assert(m_qf->getSize() == 1);

          m_p.emplace(polytope, m_qf->getPoint(0));
          m_weight = m_qf->getWeight(0);

          const auto& rc = m_qf->getPoint(0);

          if constexpr (std::is_same_v<MultiplicandRangeType, ScalarType>)
          {
            const P1Element<MultiplicandRangeType> trialfe(geometry);
            m_trialBasis.resize(trialfe.getCount());
            for (size_t i = 0; i < trialfe.getCount(); ++i)
              m_trialBasis[i] = trialfe.getBasis(i)(rc);
          }
          else if constexpr (std::is_same_v<MultiplicandRangeType, Math::Vector<ScalarType>>)
          {
            const P1Element<MultiplicandRangeType> trialfe(geometry, trialfes.getVectorDimension());
            m_trialBasis.resize(trialfe.getCount());
            for (size_t i = 0; i < trialfe.getCount(); ++i)
              m_trialBasis[i] = trialfe.getBasis(i)(rc);
          }
          else
          {
            assert(false);
          }

          if (trialfes == testfes)
          {
            m_testBasis = m_trialBasis;
          }
          else
          {
            if constexpr (std::is_same_v<RHSRangeType, ScalarType>)
            {
              const P1Element<RHSRangeType> testfe(geometry);
              m_testBasis.resize(testfe.getCount());
              for (size_t i = 0; i < testfe.getCount(); ++i)
                m_testBasis[i] = testfe.getBasis(i)(rc);
            }
            else if constexpr (std::is_same_v<RHSRangeType, Math::Vector<ScalarType>>)
            {
              const P1Element<RHSRangeType> testfe(geometry, testfes.getVectorDimension());
              m_testBasis.resize(testfe.getCount());
              for (size_t i = 0; i < testfe.getCount(); ++i)
                m_testBasis[i] = testfe.getBasis(i)(rc);
            }
            else
            {
              assert(false);
            }
          }

          m_matrix.resize(m_testBasis.size(), m_trialBasis.size());
        }

        assert(m_p);
        auto& p = *m_p;
        p.setPolytope(polytope);
        m_distortion = p.getDistortion();

        if constexpr (std::is_same_v<CoefficientRangeType, ScalarType>)
        {
          const ScalarType csv = coeff.getValue(p);

          if (trialfes == testfes)
          {
            const size_t n = m_trialBasis.size();
            for (size_t i = 0; i < n; ++i)
            {
              m_matrix(i, i) = csv * Math::dot(m_trialBasis[i], m_trialBasis[i]);
              for (size_t j = 0; j < i; ++j)
                m_matrix(i, j) = csv * Math::dot(m_trialBasis[j], m_trialBasis[i]);
            }

            m_matrix.template triangularView<Eigen::Upper>() =
              m_matrix.template triangularView<Eigen::Lower>().transpose();
          }
          else
          {
            const size_t ntr = m_trialBasis.size();
            const size_t nte = m_testBasis.size();

            for (size_t te = 0; te < nte; ++te)
              for (size_t tr = 0; tr < ntr; ++tr)
                m_matrix(te, tr) = csv * Math::dot(m_testBasis[te], m_trialBasis[tr]);
          }
        }
        else if constexpr (std::is_same_v<CoefficientRangeType, Math::Matrix<ScalarType>>)
        {
          static_assert(std::is_same_v<MultiplicandRangeType, Math::Vector<ScalarType>>);
          static_assert(std::is_same_v<RHSRangeType, Math::Vector<ScalarType>>);

          coeff.getValue(m_cmv, p);

          const size_t ntr = m_trialBasis.size();
          const size_t nte = m_testBasis.size();

          for (size_t te = 0; te < nte; ++te)
            for (size_t tr = 0; tr < ntr; ++tr)
              m_matrix(te, tr) = Math::dot(m_cmv * m_trialBasis[tr], m_testBasis[te]);
        }
        else
        {
          assert(false);
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_weight * m_distortion * m_matrix(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Optional<QF::Centroid> m_qf;
      Optional<Geometry::Point> m_p;

      Real m_distortion;
      Real m_weight;

      Math::Matrix<ScalarType> m_cmv;

      std::vector<MultiplicandRangeType> m_trialBasis;
      std::vector<RHSRangeType> m_testBasis;

      Math::Matrix<ScalarType> m_matrix;

      bool m_set;
      Optional<Geometry::Polytope::Type> m_geometry;
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
      using LHSFESType = P1<LHSRange, LHSMesh>;

      using RHSFESType = P1<RHSRange, RHSMesh>;

      using LHSType =
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, LHSFESType, TrialSpace>>>;

      using LHSOperandType =
        ShapeFunction<LHSDerived, LHSFESType, TrialSpace>;

      using LHSOperandRangeType =
        typename FormLanguage::Traits<LHSOperandType>::RangeType;

      using RHSType =
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, RHSFESType, TestSpace>>>;

      using RHSOperandType =
        ShapeFunction<RHSDerived, RHSFESType, TestSpace>;

      using RHSOperandRangeType =
        typename FormLanguage::Traits<RHSOperandType>::RangeType;

      using IntegrandType = Dot<LHSType, RHSType>;

      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSOperandRangeType, ScalarType>);
      static_assert(std::is_same_v<RHSOperandRangeType, ScalarType>);

      constexpr
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_set(false)
      {}

      constexpr
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_set(false)
      {}

      constexpr
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_polytope(std::move(other.m_polytope)),
          m_qf(std::move(other.m_qf)),
          m_p(std::move(other.m_p)),
          m_weight(std::move(other.m_weight)),
          m_distortion(std::move(other.m_distortion)),
          m_grad(std::move(other.m_grad)),
          m_matrix(std::move(other.m_matrix)),
          m_set(std::move(other.m_set)),
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
        const bool recompute = !m_set || m_geometry != geometry;

        if (recompute)
        {
          m_set = true;
          m_geometry = geometry;

          const size_t d = polytope.getDimension();
          const P1Element<ScalarType> fe(geometry);
          const size_t n = fe.getCount();

          m_qf.emplace(geometry);
          assert(m_qf->getSize() == 1);

          m_p.emplace(polytope, m_qf->getPoint(0));
          m_weight = m_qf->getWeight(0);

          const auto& rc = m_qf->getPoint(0);

          m_grad.resize(n);
          m_matrix.resize(n, n);

          for (size_t local = 0; local < n; ++local)
          {
            auto& g = m_grad[local];
            g.resize(d);
            const auto& basis = fe.getBasis(local);
            for (size_t j = 0; j < d; ++j)
              g(j) = basis.template getDerivative<1>(j)(rc);
          }
        }

        assert(m_p);
        auto& p = *m_p;
        p.setPolytope(polytope);

        const auto& Jinv = p.getJacobianInverse();
        const auto G = Jinv * Jinv.transpose();

        const size_t n = m_grad.size();

        for (size_t i = 0; i < n; ++i)
        {
          const auto Ggi = G * m_grad[i];
          m_matrix(i, i) = Math::dot(m_grad[i], Ggi);
          for (size_t j = 0; j < i; ++j)
            m_matrix(i, j) = Math::dot(m_grad[j], Ggi);
        }

        m_matrix.template triangularView<Eigen::Upper>() =
          m_matrix.template triangularView<Eigen::Lower>().transpose();

        m_distortion = p.getDistortion();
        return *this;
      }


      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_weight * m_distortion * m_matrix(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Optional<QF::Centroid> m_qf;
      Optional<Geometry::Point> m_p;

      Real m_weight;
      Real m_distortion;

      std::vector<Math::SpatialVector<ScalarType>> m_grad;
      Math::Matrix<ScalarType> m_matrix;

      bool m_set;
      Optional<Geometry::Polytope::Type> m_geometry;
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

      using LHSType =
        ShapeFunctionBase<Mult<CoefficientType, MultiplicandType>>;

      using MultiplicandOperandRangeType =
        typename FormLanguage::Traits<MultiplicandOperandType>::RangeType;

      using RHSType =
        ShapeFunctionBase<
          Grad<ShapeFunction<RHSDerived, RHSFESType, TestSpace>>>;

      using RHSOperandType =
        ShapeFunction<RHSDerived, RHSFESType, TestSpace>;

      using RHSOperandRangeType =
        typename FormLanguage::Traits<RHSOperandType>::RangeType;

      using IntegrandType = Dot<LHSType, RHSType>;

      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<MultiplicandOperandRangeType, ScalarType>);
      static_assert(std::is_same_v<RHSOperandRangeType, ScalarType>);

      constexpr
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_weight(Math::nan<Real>()),
          m_distortion(Math::nan<Real>()),
          m_set(false)
      {}

      constexpr
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_weight(Math::nan<Real>()),
          m_distortion(Math::nan<Real>()),
          m_set(false)
      {}

      constexpr
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_polytope(std::move(other.m_polytope)),
          m_qf(std::move(other.m_qf)),
          m_p(std::move(other.m_p)),
          m_weight(std::move(other.m_weight)),
          m_distortion(std::move(other.m_distortion)),
          m_cmv(std::move(other.m_cmv)),
          m_trialGrad(std::move(other.m_trialGrad)),
          m_testGrad(std::move(other.m_testGrad)),
          m_matrix(std::move(other.m_matrix)),
          m_set(std::move(other.m_set)),
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
        const bool recompute = !m_set || m_geometry != geometry;

        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& coeff = lhs.getDerived().getLHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        if (recompute)
        {
          m_set = true;
          m_geometry = geometry;

          const size_t d = polytope.getDimension();

          const P1Element<ScalarType> trialfe(geometry);
          const P1Element<ScalarType> testfe(geometry);

          m_qf.emplace(geometry);
          assert(m_qf->getSize() == 1);

          m_p.emplace(polytope, m_qf->getPoint(0));
          m_weight = m_qf->getWeight(0);

          const auto& rc = m_qf->getPoint(0);

          m_trialGrad.resize(trialfe.getCount());
          for (size_t local = 0; local < trialfe.getCount(); ++local)
          {
            auto& g = m_trialGrad[local];
            g.resize(d);
            const auto& basis = trialfe.getBasis(local);
            for (size_t j = 0; j < d; ++j)
              g(j) = basis.template getDerivative<1>(j)(rc);
          }

          if (trialfes == testfes)
          {
            m_testGrad = m_trialGrad;
          }
          else
          {
            m_testGrad.resize(testfe.getCount());
            for (size_t local = 0; local < testfe.getCount(); ++local)
            {
              auto& g = m_testGrad[local];
              g.resize(d);
              const auto& basis = testfe.getBasis(local);
              for (size_t j = 0; j < d; ++j)
                g(j) = basis.template getDerivative<1>(j)(rc);
            }
          }

          m_matrix.resize(m_testGrad.size(), m_trialGrad.size());
        }

        assert(m_p);
        auto& p = *m_p;
        p.setPolytope(polytope);
        m_distortion = p.getDistortion();

        const auto& Jinv = p.getJacobianInverse();
        const auto G = Jinv * Jinv.transpose();

        if constexpr (std::is_same_v<CoefficientRangeType, ScalarType>)
        {
          const ScalarType csv = coeff.getValue(p);

          if (trialfes == testfes)
          {
            const size_t n = m_trialGrad.size();
            for (size_t i = 0; i < n; ++i)
            {
              const auto Ggi = G * m_trialGrad[i];
              m_matrix(i, i) = csv * Math::dot(m_trialGrad[i], Ggi);
              for (size_t j = 0; j < i; ++j)
                m_matrix(i, j) = csv * Math::dot(m_trialGrad[j], Ggi);
            }

            m_matrix.template triangularView<Eigen::Upper>() =
              m_matrix.template triangularView<Eigen::Lower>().transpose();
          }
          else
          {
            const size_t ntr = m_trialGrad.size();
            const size_t nte = m_testGrad.size();

            for (size_t te = 0; te < nte; ++te)
            {
              const auto Ggte = G * m_testGrad[te];
              for (size_t tr = 0; tr < ntr; ++tr)
                m_matrix(te, tr) = csv * Math::dot(m_trialGrad[tr], Ggte);
            }
          }
        }
        else if constexpr (std::is_same_v<CoefficientRangeType, Math::Matrix<ScalarType>>)
        {
          coeff.getValue(m_cmv, p);

          if (trialfes == testfes)
          {
            const size_t n = m_trialGrad.size();
            for (size_t i = 0; i < n; ++i)
            {
              const auto AGgi = m_cmv * (G * m_trialGrad[i]);
              m_matrix(i, i) = Math::dot(AGgi, m_trialGrad[i]);
              for (size_t j = 0; j < i; ++j)
                m_matrix(i, j) = Math::dot(AGgi, m_trialGrad[j]);
            }

            // Keep the safe generic behavior only if you know A is symmetric.
            // Here we fill the full matrix instead of mirroring.
            for (size_t i = 0; i < n; ++i)
              for (size_t j = i + 1; j < n; ++j)
                m_matrix(i, j) = Math::dot(m_cmv * (G * m_trialGrad[j]), m_trialGrad[i]);
          }
          else
          {
            const size_t ntr = m_trialGrad.size();
            const size_t nte = m_testGrad.size();

            for (size_t te = 0; te < nte; ++te)
            {
              const auto Ggte = G * m_testGrad[te];
              for (size_t tr = 0; tr < ntr; ++tr)
                m_matrix(te, tr) = Math::dot(m_cmv * (G * m_trialGrad[tr]), m_testGrad[te]);
            }
          }
        }
        else
        {
          assert(false);
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_weight * m_distortion * m_matrix(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Optional<QF::Centroid> m_qf;
      Optional<Geometry::Point> m_p;

      Real m_weight;
      Real m_distortion;

      Math::Matrix<ScalarType> m_cmv;

      std::vector<Math::SpatialVector<ScalarType>> m_trialGrad;
      std::vector<Math::SpatialVector<ScalarType>> m_testGrad;

      Math::Matrix<ScalarType> m_matrix;

      bool m_set;
      Optional<Geometry::Polytope::Type> m_geometry;
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
      using LHSFESType = P1<LHSRange, LHSMesh>;
      using RHSFESType = P1<RHSRange, RHSMesh>;

      using CoefficientType = FunctionBase<CoefficientDerived>;

      using LHSType =
        ShapeFunctionBase<ShapeFunction<LHSDerived, LHSFESType, TrialSpace>>;

      using RHSType =
        ShapeFunctionBase<ShapeFunction<RHSDerived, RHSFESType, TestSpace>>;

      using InnerIntegrandType = Dot<LHSType, RHSType>;
      using IntegrandType = Mult<CoefficientType, InnerIntegrandType>;

      using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;
      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;
      using ScalarType   = typename FormLanguage::Traits<InnerIntegrandType>::ScalarType;

      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);

      constexpr
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getRHS().getLHS().getLeaf(), integrand.getRHS().getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_distortion(Math::nan<Real>()),
          m_weight(Math::nan<Real>()),
          m_set(false)
      {}

      constexpr
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_distortion(Math::nan<Real>()),
          m_weight(Math::nan<Real>()),
          m_set(false)
      {}

      constexpr
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_polytope(std::move(other.m_polytope)),
          m_qf(std::move(other.m_qf)),
          m_p(std::move(other.m_p)),
          m_distortion(std::move(other.m_distortion)),
          m_weight(std::move(other.m_weight)),
          m_trialBasis(std::move(other.m_trialBasis)),
          m_testBasis(std::move(other.m_testBasis)),
          m_matrix(std::move(other.m_matrix)),
          m_set(std::move(other.m_set)),
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
        const bool recompute = !m_set || m_geometry != geometry;

        const auto& integrand = getIntegrand();
        const auto& coeff = integrand.getLHS();
        const auto& inner = integrand.getRHS();
        const auto& lhs = inner.getLHS();
        const auto& rhs = inner.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        if (recompute)
        {
          m_set = true;
          m_geometry = geometry;

          m_qf.emplace(geometry);
          assert(m_qf->getSize() == 1);

          m_p.emplace(polytope, m_qf->getPoint(0));
          m_weight = m_qf->getWeight(0);

          const auto& rc = m_qf->getPoint(0);

          if constexpr (std::is_same_v<LHSRangeType, ScalarType>)
          {
            const P1Element<LHSRangeType> trialfe(geometry);
            m_trialBasis.resize(trialfe.getCount());
            for (size_t i = 0; i < trialfe.getCount(); ++i)
              m_trialBasis[i] = trialfe.getBasis(i)(rc);
          }
          else if constexpr (std::is_same_v<LHSRangeType, Math::Vector<ScalarType>>)
          {
            const P1Element<LHSRangeType> trialfe(geometry, trialfes.getVectorDimension());
            m_trialBasis.resize(trialfe.getCount());
            for (size_t i = 0; i < trialfe.getCount(); ++i)
              m_trialBasis[i] = trialfe.getBasis(i)(rc);
          }
          else
          {
            assert(false);
          }

          if (trialfes == testfes)
          {
            m_testBasis = m_trialBasis;
          }
          else
          {
            if constexpr (std::is_same_v<RHSRangeType, ScalarType>)
            {
              const P1Element<RHSRangeType> testfe(geometry);
              m_testBasis.resize(testfe.getCount());
              for (size_t i = 0; i < testfe.getCount(); ++i)
                m_testBasis[i] = testfe.getBasis(i)(rc);
            }
            else if constexpr (std::is_same_v<RHSRangeType, Math::Vector<ScalarType>>)
            {
              const P1Element<RHSRangeType> testfe(geometry, testfes.getVectorDimension());
              m_testBasis.resize(testfe.getCount());
              for (size_t i = 0; i < testfe.getCount(); ++i)
                m_testBasis[i] = testfe.getBasis(i)(rc);
            }
            else
            {
              assert(false);
            }
          }

          m_matrix.resize(m_testBasis.size(), m_trialBasis.size());
        }

        assert(m_p);
        auto& p = *m_p;
        p.setPolytope(polytope);
        m_distortion = p.getDistortion();

        const ScalarType csv = coeff(p);

        const size_t ntr = m_trialBasis.size();
        const size_t nte = m_testBasis.size();

        for (size_t te = 0; te < nte; ++te)
          for (size_t tr = 0; tr < ntr; ++tr)
          {
            if constexpr (std::is_same_v<LHSRangeType, ScalarType>)
              m_matrix(te, tr) = csv * m_testBasis[te] * m_trialBasis[tr];
            else
              m_matrix(te, tr) = csv * Math::dot(m_testBasis[te], m_trialBasis[tr]);
          }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_weight * m_distortion * m_matrix(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Optional<QF::Centroid> m_qf;
      Optional<Geometry::Point> m_p;

      Real m_distortion;
      Real m_weight;

      std::vector<LHSRangeType> m_trialBasis;
      std::vector<RHSRangeType> m_testBasis;
      Math::Matrix<ScalarType> m_matrix;

      bool m_set;
      Optional<Geometry::Polytope::Type> m_geometry;
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
        Div<ShapeFunction<LHSDerived, P1<Math::Vector<Real>, LHSMesh>, TrialSpace>>,
        P1<Math::Vector<Real>, LHSMesh>, TrialSpace>,
      ShapeFunctionBase<
        ShapeFunction<RHSDerived, P1<Real, RHSMesh>, TestSpace>,
        P1<Real, RHSMesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<typename FormLanguage::Traits<P1<Real, LHSMesh>>::ScalarType>
  {
    public:
      using ScalarType = typename FormLanguage::Traits<P1<Real, LHSMesh>>::ScalarType;
      using TrialFESType = P1<Math::Vector<Real>, LHSMesh>;
      using TestFESType  = P1<Real, RHSMesh>;

      using LHSType = ShapeFunctionBase<
        Div<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>>,
        TrialFESType, TrialSpace>;

      using RHSType = ShapeFunctionBase<
        ShapeFunction<RHSDerived, TestFESType, TestSpace>,
        TestFESType, TestSpace>;

      using IntegrandType = Dot<LHSType, RHSType>;
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      constexpr
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_distortion(Math::nan<Real>()),
          m_weight(Math::nan<Real>()),
          m_set(false)
      {}

      constexpr
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_distortion(Math::nan<Real>()),
          m_weight(Math::nan<Real>()),
          m_set(false)
      {}

      constexpr
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_polytope(std::move(other.m_polytope)),
          m_qf(std::move(other.m_qf)),
          m_p(std::move(other.m_p)),
          m_distortion(std::move(other.m_distortion)),
          m_weight(std::move(other.m_weight)),
          m_trialDiv(std::move(other.m_trialDiv)),
          m_testBasis(std::move(other.m_testBasis)),
          m_matrix(std::move(other.m_matrix)),
          m_set(std::move(other.m_set)),
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
        const bool recompute = !m_set || m_geometry != geometry;

        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes  = rhs.getFiniteElementSpace();

        if (recompute)
        {
          m_set = true;
          m_geometry = geometry;

          m_qf.emplace(geometry);
          assert(m_qf->getSize() == 1);

          m_p.emplace(polytope, m_qf->getPoint(0));
          m_weight = m_qf->getWeight(0);

          const auto& rc = m_qf->getPoint(0);

          const P1Element<Math::Vector<Real>> trialfe(geometry, trialfes.getVectorDimension());
          const P1Element<Real> testfe(geometry);

          m_trialDiv.resize(trialfe.getCount());
          m_testBasis.resize(testfe.getCount());
          for (size_t i = 0; i < testfe.getCount(); ++i)
            m_testBasis[i] = testfe.getBasis(i)(rc);

          m_matrix.resize(m_testBasis.size(), m_trialDiv.size());

          // precompute reference gradients for divergence
          m_refGrad.resize(trialfe.getCount());
          for (size_t i = 0; i < trialfe.getCount(); ++i)
          {
            auto& refg = m_refGrad[i];
            refg.resize(trialfes.getVectorDimension());
            for (size_t comp = 0; comp < trialfes.getVectorDimension(); ++comp)
            {
              refg[comp].resize(polytope.getDimension());
              const auto& basis = trialfe.getBasis(i);
              for (size_t j = 0; j < polytope.getDimension(); ++j)
                refg[comp](j) = basis.template getDerivative<1>(comp, j)(rc);
            }
          }
        }

        assert(m_p);
        auto& p = *m_p;
        p.setPolytope(polytope);
        m_distortion = p.getDistortion();

        const auto& Jinv = p.getJacobianInverse();
        const size_t d = polytope.getDimension();
        const size_t vdim = m_refGrad.empty() ? 0 : m_refGrad.front().size();

        for (size_t i = 0; i < m_refGrad.size(); ++i)
        {
          ScalarType div = 0;
          for (size_t comp = 0; comp < std::min(vdim, d); ++comp)
          {
            ScalarType physComp = 0;
            for (size_t j = 0; j < d; ++j)
              physComp += Jinv(comp, j) * m_refGrad[i][comp](j);
            div += physComp;
          }
          m_trialDiv[i] = div;
        }

        const size_t ntr = m_trialDiv.size();
        const size_t nte = m_testBasis.size();
        for (size_t te = 0; te < nte; ++te)
          for (size_t tr = 0; tr < ntr; ++tr)
            m_matrix(te, tr) = m_testBasis[te] * m_trialDiv[tr];

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_weight * m_distortion * m_matrix(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Optional<QF::Centroid> m_qf;
      Optional<Geometry::Point> m_p;

      ScalarType m_distortion;
      ScalarType m_weight;

      std::vector<ScalarType> m_trialDiv;
      std::vector<ScalarType> m_testBasis;
      Math::Matrix<ScalarType> m_matrix;

      std::vector<std::vector<Math::SpatialVector<ScalarType>>> m_refGrad;

      bool m_set;
      Optional<Geometry::Polytope::Type> m_geometry;
  };

  template <class LHSDerived, class RHSDerived, class LHSMesh, class RHSMesh>
  using P1DivTrialIntegrand =
    Dot<
      ShapeFunctionBase<
        Div<ShapeFunction<LHSDerived, P1<Math::Vector<Real>, LHSMesh>, TrialSpace>>,
        P1<Math::Vector<Real>, LHSMesh>, TrialSpace>,
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
        Div<ShapeFunction<RHSDerived, P1<Math::Vector<Real>, RHSMesh>, TestSpace>>,
        P1<Math::Vector<Real>, RHSMesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<typename FormLanguage::Traits<P1<Real, LHSMesh>>::ScalarType>
  {
    public:
      using ScalarType = typename FormLanguage::Traits<P1<Real, LHSMesh>>::ScalarType;
      using TrialFESType = P1<Real, LHSMesh>;
      using TestFESType  = P1<Math::Vector<Real>, RHSMesh>;

      using LHSType = ShapeFunctionBase<
        ShapeFunction<LHSDerived, TrialFESType, TrialSpace>,
        TrialFESType, TrialSpace>;

      using RHSType = ShapeFunctionBase<
        Div<ShapeFunction<RHSDerived, TestFESType, TestSpace>>,
        TestFESType, TestSpace>;

      using IntegrandType = Dot<LHSType, RHSType>;
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      constexpr
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_distortion(Math::nan<Real>()),
          m_weight(Math::nan<Real>()),
          m_set(false)
      {}

      constexpr
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_distortion(Math::nan<Real>()),
          m_weight(Math::nan<Real>()),
          m_set(false)
      {}

      constexpr
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_polytope(std::move(other.m_polytope)),
          m_qf(std::move(other.m_qf)),
          m_p(std::move(other.m_p)),
          m_distortion(std::move(other.m_distortion)),
          m_weight(std::move(other.m_weight)),
          m_trialBasis(std::move(other.m_trialBasis)),
          m_testDiv(std::move(other.m_testDiv)),
          m_matrix(std::move(other.m_matrix)),
          m_set(std::move(other.m_set)),
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
        const bool recompute = !m_set || m_geometry != geometry;

        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes  = rhs.getFiniteElementSpace();

        if (recompute)
        {
          m_set = true;
          m_geometry = geometry;

          m_qf.emplace(geometry);
          assert(m_qf->getSize() == 1);

          m_p.emplace(polytope, m_qf->getPoint(0));
          m_weight = m_qf->getWeight(0);

          const auto& rc = m_qf->getPoint(0);

          const P1Element<Real> trialfe(geometry);
          const P1Element<Math::Vector<Real>> testfe(geometry, testfes.getVectorDimension());

          m_trialBasis.resize(trialfe.getCount());
          for (size_t i = 0; i < trialfe.getCount(); ++i)
            m_trialBasis[i] = trialfe.getBasis(i)(rc);

          m_testDiv.resize(testfe.getCount());
          m_matrix.resize(m_testDiv.size(), m_trialBasis.size());

          m_refGrad.resize(testfe.getCount());
          for (size_t i = 0; i < testfe.getCount(); ++i)
          {
            auto& refg = m_refGrad[i];
            refg.resize(testfes.getVectorDimension());
            for (size_t comp = 0; comp < testfes.getVectorDimension(); ++comp)
            {
              refg[comp].resize(polytope.getDimension());
              const auto& basis = testfe.getBasis(i);
              for (size_t j = 0; j < polytope.getDimension(); ++j)
                refg[comp](j) = basis.template getDerivative<1>(comp, j)(rc);
            }
          }
        }

        assert(m_p);
        auto& p = *m_p;
        p.setPolytope(polytope);
        m_distortion = p.getDistortion();

        const auto& Jinv = p.getJacobianInverse();
        const size_t d = polytope.getDimension();
        const size_t vdim = m_refGrad.empty() ? 0 : m_refGrad.front().size();

        for (size_t i = 0; i < m_refGrad.size(); ++i)
        {
          ScalarType div = 0;
          for (size_t comp = 0; comp < std::min(vdim, d); ++comp)
          {
            ScalarType physComp = 0;
            for (size_t j = 0; j < d; ++j)
              physComp += Jinv(comp, j) * m_refGrad[i][comp](j);
            div += physComp;
          }
          m_testDiv[i] = div;
        }

        const size_t ntr = m_trialBasis.size();
        const size_t nte = m_testDiv.size();
        for (size_t te = 0; te < nte; ++te)
          for (size_t tr = 0; tr < ntr; ++tr)
            m_matrix(te, tr) = m_trialBasis[tr] * m_testDiv[te];

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_weight * m_distortion * m_matrix(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Optional<QF::Centroid> m_qf;
      Optional<Geometry::Point> m_p;

      ScalarType m_distortion;
      ScalarType m_weight;

      std::vector<ScalarType> m_trialBasis;
      std::vector<ScalarType> m_testDiv;
      Math::Matrix<ScalarType> m_matrix;

      std::vector<std::vector<Math::SpatialVector<ScalarType>>> m_refGrad;

      bool m_set;
      Optional<Geometry::Polytope::Type> m_geometry;
  };

  template <class LHSDerived, class RHSDerived, class LHSMesh, class RHSMesh>
  using P1DivTestIntegrand =
    Dot<
      ShapeFunctionBase<
        ShapeFunction<LHSDerived, P1<Real, LHSMesh>, TrialSpace>,
        P1<Real, LHSMesh>, TrialSpace>,
      ShapeFunctionBase<
        Div<ShapeFunction<RHSDerived, P1<Math::Vector<Real>, RHSMesh>, TestSpace>>,
        P1<Math::Vector<Real>, RHSMesh>, TestSpace>>;

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
      using LHSFESType = P1<LHSRange, LHSMesh>;

      using RHSFESType = P1<RHSRange, RHSMesh>;

      using LHSType =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<LHSDerived, LHSFESType, TrialSpace>>>;

      using LHSOperandType =
        ShapeFunction<LHSDerived, LHSFESType, TrialSpace>;

      using LHSOperandRangeType =
        typename FormLanguage::Traits<LHSOperandType>::RangeType;

      using RHSType =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<RHSDerived, RHSFESType, TestSpace>>>;

      using RHSOperandType =
        ShapeFunction<RHSDerived, RHSFESType, TestSpace>;

      using RHSOperandRangeType =
        typename FormLanguage::Traits<RHSOperandType>::RangeType;

      using IntegrandType = Dot<LHSType, RHSType>;

      using IntegrandRangeType =
        typename FormLanguage::Traits<IntegrandType>::RangeType;

      using ScalarType =
        typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent =
        LocalBilinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSOperandRangeType, Math::Vector<ScalarType>>);
      static_assert(std::is_same_v<RHSOperandRangeType, Math::Vector<ScalarType>>);

      constexpr
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_distortion(Math::nan<Real>()),
          m_weight(Math::nan<Real>()),
          m_set(false)
      {}

      constexpr
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_distortion(Math::nan<Real>()),
          m_weight(Math::nan<Real>()),
          m_set(false)
      {}

      constexpr
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_polytope(std::move(other.m_polytope)),
          m_qf(std::move(other.m_qf)),
          m_p(std::move(other.m_p)),
          m_distortion(std::move(other.m_distortion)),
          m_weight(std::move(other.m_weight)),
          m_trialJac(std::move(other.m_trialJac)),
          m_testJac(std::move(other.m_testJac)),
          m_matrix(std::move(other.m_matrix)),
          m_set(std::move(other.m_set)),
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
        const bool recompute = !m_set || m_geometry != geometry;

        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        if (recompute)
        {
          m_set = true;
          m_geometry = geometry;

          const size_t d = polytope.getDimension();

          const P1Element<Math::Vector<ScalarType>> trialfe(geometry, trialfes.getVectorDimension());
          const P1Element<Math::Vector<ScalarType>> testfe(geometry, testfes.getVectorDimension());

          m_qf.emplace(geometry);
          assert(m_qf->getSize() == 1);

          m_p.emplace(polytope, m_qf->getPoint(0));
          m_weight = m_qf->getWeight(0);

          const auto& rc = m_qf->getPoint(0);

          m_trialJac.resize(trialfe.getCount());
          for (size_t local = 0; local < trialfe.getCount(); ++local)
          {
            auto& J = m_trialJac[local];
            J.resize(trialfes.getVectorDimension(), d);
            const auto& basis = trialfe.getBasis(local);
            for (size_t i = 0; i < trialfes.getVectorDimension(); ++i)
            {
              for (size_t j = 0; j < d; ++j)
                J(i, j) = basis.template getDerivative<1>(i, j)(rc);
            }
          }

          if (trialfes == testfes)
          {
            m_testJac = m_trialJac;
          }
          else
          {
            m_testJac.resize(testfe.getCount());
            for (size_t local = 0; local < testfe.getCount(); ++local)
            {
              auto& J = m_testJac[local];
              J.resize(testfes.getVectorDimension(), d);
              const auto& basis = testfe.getBasis(local);
              for (size_t i = 0; i < testfes.getVectorDimension(); ++i)
              {
                for (size_t j = 0; j < d; ++j)
                  J(i, j) = basis.template getDerivative<1>(i, j)(rc);
              }
            }
          }

          m_matrix.resize(m_testJac.size(), m_trialJac.size());
        }

        assert(m_p);
        auto& p = *m_p;
        p.setPolytope(polytope);
        m_distortion = p.getDistortion();

        const auto& Jinv = p.getJacobianInverse();

        if (trialfes == testfes)
        {
          const size_t n = m_trialJac.size();
          for (size_t i = 0; i < n; ++i)
          {
            const auto Ji = m_trialJac[i] * Jinv;
            m_matrix(i, i) = Ji.squaredNorm();

            for (size_t j = 0; j < i; ++j)
              m_matrix(i, j) = Math::dot(m_trialJac[j] * Jinv, Ji);
          }

          m_matrix.template triangularView<Eigen::Upper>() =
            m_matrix.template triangularView<Eigen::Lower>().transpose();
        }
        else
        {
          const size_t ntr = m_trialJac.size();
          const size_t nte = m_testJac.size();

          for (size_t te = 0; te < nte; ++te)
          {
            const auto Jte = m_testJac[te] * Jinv;
            for (size_t tr = 0; tr < ntr; ++tr)
              m_matrix(te, tr) = Math::dot(m_trialJac[tr] * Jinv, Jte);
          }
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_weight * m_distortion * m_matrix(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Optional<QF::Centroid> m_qf;
      Optional<Geometry::Point> m_p;

      Real m_distortion;
      Real m_weight;

      std::vector<Math::SpatialMatrix<ScalarType>> m_trialJac, m_testJac;

      Math::Matrix<ScalarType> m_matrix;

      bool m_set;
      Optional<Geometry::Polytope::Type> m_geometry;
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
      using LHSFESType = P1<LHSRange, LHSMesh>;

      using RHSFESType = P1<RHSRange, RHSMesh>;

      using CoefficientType = FunctionBase<CoefficientDerived>;

      using CoefficientRangeType =
        typename FormLanguage::Traits<CoefficientType>::RangeType;

      using MultiplicandType =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<LHSDerived, LHSFESType, TrialSpace>>>;

      using LHSType =
        ShapeFunctionBase<Mult<CoefficientType, MultiplicandType>>;

      using LHSOperandType =
        ShapeFunction<LHSDerived, LHSFESType, TrialSpace>;

      using LHSOperandRangeType =
        typename FormLanguage::Traits<LHSOperandType>::RangeType;

      using RHSType =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>>>;

      using RHSOperandType =
        ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>;

      using RHSOperandRangeType =
        typename FormLanguage::Traits<RHSOperandType>::RangeType;

      using IntegrandType = Dot<LHSType, RHSType>;

      using IntegrandRangeType =
        typename FormLanguage::Traits<IntegrandType>::RangeType;

      using ScalarType =
        typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent =
        LocalBilinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSOperandRangeType, Math::Vector<ScalarType>>);
      static_assert(std::is_same_v<RHSOperandRangeType, Math::Vector<ScalarType>>);

      constexpr
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_weight(Math::nan<Real>()),
          m_distortion(Math::nan<Real>()),
          m_set(false)
      {}

      constexpr
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_weight(Math::nan<Real>()),
          m_distortion(Math::nan<Real>()),
          m_set(false)
      {}

      constexpr
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_polytope(std::move(other.m_polytope)),
          m_qf(std::move(other.m_qf)),
          m_p(std::move(other.m_p)),
          m_weight(std::move(other.m_weight)),
          m_distortion(std::move(other.m_distortion)),
          m_cmv(std::move(other.m_cmv)),
          m_trialJac(std::move(other.m_trialJac)),
          m_testJac(std::move(other.m_testJac)),
          m_matrix(std::move(other.m_matrix)),
          m_set(std::move(other.m_set)),
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
        const bool recompute = !m_set || m_geometry != geometry;

        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& coeff = lhs.getDerived().getLHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        if (recompute)
        {
          m_set = true;
          m_geometry = geometry;

          const size_t d = polytope.getDimension();

          const P1Element<Math::Vector<ScalarType>> trialfe(geometry, trialfes.getVectorDimension());
          const P1Element<Math::Vector<ScalarType>> testfe(geometry, testfes.getVectorDimension());

          m_qf.emplace(geometry);
          assert(m_qf->getSize() == 1);

          m_p.emplace(polytope, m_qf->getPoint(0));
          m_weight = m_qf->getWeight(0);

          const auto& rc = m_qf->getPoint(0);

          m_trialJac.resize(trialfe.getCount());
          for (size_t local = 0; local < trialfe.getCount(); ++local)
          {
            auto& J = m_trialJac[local];
            J.resize(trialfes.getVectorDimension(), d);
            const auto& basis = trialfe.getBasis(local);
            for (size_t i = 0; i < trialfes.getVectorDimension(); ++i)
            {
              for (size_t j = 0; j < d; ++j)
                J(i, j) = basis.template getDerivative<1>(i, j)(rc);
            }
          }

          if (trialfes == testfes)
          {
            m_testJac = m_trialJac;
          }
          else
          {
            m_testJac.resize(testfe.getCount());
            for (size_t local = 0; local < testfe.getCount(); ++local)
            {
              auto& J = m_testJac[local];
              J.resize(testfes.getVectorDimension(), d);
              const auto& basis = testfe.getBasis(local);
              for (size_t i = 0; i < testfes.getVectorDimension(); ++i)
              {
                for (size_t j = 0; j < d; ++j)
                  J(i, j) = basis.template getDerivative<1>(i, j)(rc);
              }
            }
          }

          m_matrix.resize(m_testJac.size(), m_trialJac.size());
        }

        assert(m_p);
        auto& p = *m_p;
        p.setPolytope(polytope);
        m_distortion = p.getDistortion();

        const auto& Jinv = p.getJacobianInverse();

        if constexpr (std::is_same_v<CoefficientRangeType, ScalarType>)
        {
          const ScalarType csv = coeff.getValue(p);

          if (trialfes == testfes)
          {
            const size_t n = m_trialJac.size();
            for (size_t i = 0; i < n; ++i)
            {
              const auto Ji = m_trialJac[i] * Jinv;
              m_matrix(i, i) = csv * Ji.squaredNorm();

              for (size_t j = 0; j < i; ++j)
                m_matrix(i, j) = csv * Math::dot(m_trialJac[j] * Jinv, Ji);
            }

            m_matrix.template triangularView<Eigen::Upper>() =
              m_matrix.template triangularView<Eigen::Lower>().transpose();
          }
          else
          {
            const size_t ntr = m_trialJac.size();
            const size_t nte = m_testJac.size();

            for (size_t te = 0; te < nte; ++te)
            {
              const auto Jte = m_testJac[te] * Jinv;
              for (size_t tr = 0; tr < ntr; ++tr)
                m_matrix(te, tr) = csv * Math::dot(m_trialJac[tr] * Jinv, Jte);
            }
          }
        }
        else if constexpr (std::is_same_v<CoefficientRangeType, Math::Matrix<ScalarType>>)
        {
          coeff.getValue(m_cmv, p);

          if (trialfes == testfes)
          {
            const size_t n = m_trialJac.size();
            for (size_t i = 0; i < n; ++i)
            {
              const auto Ji = m_trialJac[i] * Jinv;
              m_matrix(i, i) = Math::dot(m_cmv * Ji, Ji);

              for (size_t j = 0; j < i; ++j)
                m_matrix(i, j) = Math::dot(m_cmv * (m_trialJac[j] * Jinv), Ji);
            }

            for (size_t i = 0; i < n; ++i)
              for (size_t j = i + 1; j < n; ++j)
                m_matrix(i, j) = Math::dot(m_cmv * (m_trialJac[j] * Jinv), m_trialJac[i] * Jinv);
          }
          else
          {
            const size_t ntr = m_trialJac.size();
            const size_t nte = m_testJac.size();

            for (size_t te = 0; te < nte; ++te)
            {
              const auto Jte = m_testJac[te] * Jinv;
              for (size_t tr = 0; tr < ntr; ++tr)
                m_matrix(te, tr) = Math::dot(m_cmv * (m_trialJac[tr] * Jinv), Jte);
            }
          }
        }
        else
        {
          assert(false);
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_weight * m_distortion * m_matrix(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Optional<QF::Centroid> m_qf;
      Optional<Geometry::Point> m_p;

      Real m_weight;
      Real m_distortion;

      Math::Matrix<ScalarType> m_cmv;

      std::vector<Math::SpatialMatrix<ScalarType>> m_trialJac, m_testJac;

      Math::Matrix<ScalarType> m_matrix;

      bool m_set;
      Optional<Geometry::Polytope::Type> m_geometry;
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
          ShapeFunctionBase<Jacobian<ShapeFunction<LHSDerived, P1<Math::Vector<Real>, Mesh>, TrialSpace>>>>>,
      ShapeFunctionBase<
        Jacobian<ShapeFunction<RHSDerived, P1<Math::Vector<Real>, Mesh>, TestSpace>>>>&)
  ->
  QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Mult<
          FunctionBase<LHSFunctionDerived>,
          ShapeFunctionBase<
            Jacobian<ShapeFunction<LHSDerived, P1<Math::Vector<Real>, Mesh>, TrialSpace>>>>>,
      ShapeFunctionBase<
        Jacobian<ShapeFunction<RHSDerived, P1<Math::Vector<Real>, Mesh>, TestSpace>>>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   * @brief Specialization for ∫ (∇u · f)·v in the case of P1 shape functions
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
      using TrialFESType = P1<LHSRange, LHSMesh>;
      using TestFESType  = P1<RHSRange, RHSMesh>;

      using TrialSFType =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>>,
          TrialFESType, TrialSpace>;

      using CoefficientType = FunctionBase<CoefficientDerived>;

      using LHSType =
        ShapeFunctionBase<
          Mult<TrialSFType, CoefficientType>,
          TrialFESType, TrialSpace>;

      using RHSType =
        ShapeFunctionBase<
          ShapeFunction<RHSDerived, TestFESType, TestSpace>,
          TestFESType, TestSpace>;

      using IntegrandType = Dot<LHSType, RHSType>;

      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      constexpr
      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_distortion(Math::nan<Real>()),
          m_weight(Math::nan<Real>()),
          m_set(false)
      {}

      constexpr
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_distortion(Math::nan<Real>()),
          m_weight(Math::nan<Real>()),
          m_set(false)
      {}

      constexpr
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_polytope(std::move(other.m_polytope)),
          m_qf(std::move(other.m_qf)),
          m_p(std::move(other.m_p)),
          m_distortion(std::move(other.m_distortion)),
          m_weight(std::move(other.m_weight)),
          m_grad(std::move(other.m_grad)),
          m_basis(std::move(other.m_basis)),
          m_matrix(std::move(other.m_matrix)),
          m_set(std::move(other.m_set)),
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
        const bool recompute = !m_set || m_geometry != geometry;

        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();

        // The coefficient is the RHS of the Mult node
        const auto& coeff = lhs.getDerived().getRHS();

        if (recompute)
        {
          m_set = true;
          m_geometry = geometry;

          const size_t d = polytope.getDimension();
          const P1Element<ScalarType> scalarFe(geometry);
          const size_t n = scalarFe.getCount();

          m_qf.emplace(geometry);
          assert(m_qf->getSize() == 1);
          m_p.emplace(polytope, m_qf->getPoint(0));
          m_weight = m_qf->getWeight(0);

          const auto& rc = m_qf->getPoint(0);

          // Pre-compute reference gradients of scalar P1 basis
          m_grad.resize(n);
          for (size_t a = 0; a < n; ++a)
          {
            m_grad[a].resize(d);
            const auto& basis = scalarFe.getBasis(a);
            for (size_t j = 0; j < d; ++j)
              m_grad[a](j) = basis.template getDerivative<1>(j)(rc);
          }

          // Pre-compute scalar P1 basis values at centroid
          m_basis.resize(n);
          for (size_t b = 0; b < n; ++b)
            m_basis[b] = scalarFe.getBasis(b)(rc);
        }

        assert(m_p);
        auto& p = *m_p;
        p.setPolytope(polytope);
        m_distortion = p.getDistortion();

        const auto& Jinv = p.getJacobianInverse();
        const size_t n = m_grad.size();

        // Evaluate the coefficient f at the centroid
        const auto fval = coeff.getValue(p);

        // Determine vdim
        const size_t vdim = trialfes.getVectorDimension();

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);
        m_matrix.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_matrix.setZero();

        // Compute physical gradients dotted with the coefficient
        // physGrad[a] = Jinv^T * refGrad[a],  gradDotF[a] = physGrad[a] · f
        for (size_t a = 0; a < n; ++a)
        {
          const auto physGrad = Jinv.transpose() * m_grad[a];
          const ScalarType gradDotF = Math::dot(physGrad, fval);

          for (size_t b = 0; b < n; ++b)
          {
            const ScalarType val = gradDotF * m_basis[b];

            // Fill block-diagonal entries: only (c,c) blocks are non-zero
            for (size_t c = 0; c < vdim; ++c)
            {
              const size_t row = b * vdim + c;
              const size_t col = a * vdim + c;
              m_matrix(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col)) = val;
            }
          }
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_weight * m_distortion * m_matrix(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Optional<QF::Centroid> m_qf;
      Optional<Geometry::Point> m_p;

      Real m_distortion;
      Real m_weight;

      std::vector<Math::SpatialVector<ScalarType>> m_grad;
      std::vector<ScalarType> m_basis;

      Math::Matrix<ScalarType> m_matrix;

      bool m_set;
      Optional<Geometry::Polytope::Type> m_geometry;
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
      using ScalarType = typename FormLanguage::Traits<Range>::ScalarType;

      using KernelType = Kernel;

      using TrialFESType = P1<Range, Mesh>;

      using TestFESType = P1<Range, Mesh>;

      using LHSType =
        Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFESType, TrialSpace>>;

      using RHSType =
        ShapeFunctionBase<RHSDerived, TestFESType, TestSpace>;

      using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;

      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      using IntegrandType = Dot<LHSType, RHSType>;

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

      QuadratureRule& setPolytope(const Geometry::Polytope& trp, const Geometry::Polytope& tep)
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

        if constexpr (std::is_same_v<Range, ScalarType>)
        {
          const P1Element<Range> trialfe(trg);
          const P1Element<Range> testfe(teg);

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

                assert(m_qfs.getSize() == 1);
                const auto r = m_qfs.getPoint(0).value();
                const auto w = m_qfs.getWeight(0);

                const ScalarType xi = 1 - r;
                const ScalarType eta1 = r;
                const ScalarType eta2 = r;
                const ScalarType eta3 = r;

                rx0[0] = xi;
                rx0[1] = xi * (1 - eta1 + eta1 * eta2);

                rz1[0] = xi;
                rz1[1] = xi * (1 - eta1 + eta1 * eta2);

                rz2[0] = xi * (1 - eta1 * eta2);
                rz2[1] = xi * eta1 * (1 - eta2);

                rx3[0] = xi * (1 - eta1 * eta2);
                rx3[1] = xi * eta1 * (1 - eta2);

                rz4[0] = xi;
                rz4[1] = xi * eta1 * (1 - eta2);

                rx5[0] = xi;
                rx5[1] = xi * eta1 * (1 - eta2);

                rz0[0] = xi * (1 - eta1 * eta2 * eta3);
                rz0[1] = xi * (1 - eta1);

                rx1[0] = xi * (1 - eta1 * eta2 * eta3);
                rx1[1] = xi * (1 - eta1);

                rx2[0] = xi;
                rx2[1] = xi * eta1 * (1 - eta2 + eta2 * eta3);

                rz3[0] = xi;
                rz3[1] = xi * eta1 * (1 - eta2 + eta2 * eta3);

                rx4[0] = xi * (1 - eta1 * eta2 * eta3);
                rx4[1] = xi * eta1 * (1 - eta2 * eta3);

                rz5[0] = xi * (1 - eta1 * eta2 * eta3);
                rz5[1] = xi * eta1 * (1 - eta2 * eta3);

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

                m_distortion = xi * xi * xi * eta1 * eta1 * eta2;
                m_weight = w * w * w * w;
                m_matrix.resize(testfe.getCount(), trialfe.getCount());

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
                    m_matrix(l, m) = s0 * kernel(x0, z0) * trb(rx0) * teb(rz0);
                    m_matrix(l, m) += s1 * kernel(x1, z1) * trb(rx1) * teb(rz1);
                    m_matrix(l, m) += s2 * kernel(x2, z2) * trb(rx2) * teb(rz2);
                    m_matrix(l, m) += s3 * kernel(x3, z3) * trb(rx3) * teb(rz3);
                    m_matrix(l, m) += s4 * kernel(x4, z4) * trb(rx4) * teb(rz4);
                    m_matrix(l, m) += s5 * kernel(x5, z5) * trb(rx5) * teb(rz5);
                  }
                }
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
        else if constexpr (std::is_same_v<Range, Math::Vector<ScalarType>>)
        {
          const P1Element<Range> trialfe(trg, trialfes.getVectorDimension());
          const P1Element<Range> testfe(teg, testfes.getVectorDimension());

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

                assert(m_qfs.getSize() == 1);
                const auto r = m_qfs.getPoint(0).value();
                const auto w = m_qfs.getWeight(0);

                const ScalarType xi = 1 - r;
                const ScalarType eta1 = r;
                const ScalarType eta2 = r;
                const ScalarType eta3 = r;

                rx0[0] = xi;
                rx0[1] = xi * (1 - eta1 + eta1 * eta2);

                rz1[0] = xi;
                rz1[1] = xi * (1 - eta1 + eta1 * eta2);

                rz2[0] = xi * (1 - eta1 * eta2);
                rz2[1] = xi * eta1 * (1 - eta2);

                rx3[0] = xi * (1 - eta1 * eta2);
                rx3[1] = xi * eta1 * (1 - eta2);

                rz4[0] = xi;
                rz4[1] = xi * eta1 * (1 - eta2);

                rx5[0] = xi;
                rx5[1] = xi * eta1 * (1 - eta2);

                rz0[0] = xi * (1 - eta1 * eta2 * eta3);
                rz0[1] = xi * (1 - eta1);

                rx1[0] = xi * (1 - eta1 * eta2 * eta3);
                rx1[1] = xi * (1 - eta1);

                rx2[0] = xi;
                rx2[1] = xi * eta1 * (1 - eta2 + eta2 * eta3);

                rz3[0] = xi;
                rz3[1] = xi * eta1 * (1 - eta2 + eta2 * eta3);

                rx4[0] = xi * (1 - eta1 * eta2 * eta3);
                rx4[1] = xi * eta1 * (1 - eta2 * eta3);

                rz5[0] = xi * (1 - eta1 * eta2 * eta3);
                rz5[1] = xi * eta1 * (1 - eta2 * eta3);

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

                m_distortion = xi * xi * xi * eta1 * eta1 * eta2;
                m_weight = w * w * w * w;
                m_matrix.resize(testfe.getCount(), trialfe.getCount());

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
                    m_matrix(l, m) = s0 * (m_k0 * trb(rx0)).dot(teb(rz0));
                    m_matrix(l, m) += s1 * (m_k1 * trb(rx1)).dot(teb(rz1));
                    m_matrix(l, m) += s2 * (m_k2 * trb(rx2)).dot(teb(rz2));
                    m_matrix(l, m) += s3 * (m_k3 * trb(rx3)).dot(teb(rz3));
                    m_matrix(l, m) += s4 * (m_k4 * trb(rx4)).dot(teb(rz4));
                    m_matrix(l, m) += s5 * (m_k5 * trb(rx5)).dot(teb(rz5));
                  }
                }
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

      const QF::Centroid m_qfs;
      Optional<QF::Centroid> m_qftr;
      Optional<QF::Centroid> m_qfte;

      Real m_weight;
      Real m_distortion;

      ScalarType m_sk;
      Math::Matrix<ScalarType> m_mk;

      Math::Vector<ScalarType> m_trv, m_tev;
      Math::Matrix<ScalarType> m_k0, m_k1, m_k2, m_k3, m_k4, m_k5;

      Math::Matrix<ScalarType> m_matrix;
  };
}

#endif
