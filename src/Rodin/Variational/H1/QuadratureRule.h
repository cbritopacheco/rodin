#ifndef RODIN_VARIATIONAL_H1_QUADRATURERULE_H
#define RODIN_VARIATIONAL_H1_QUADRATURERULE_H

#include "Rodin/QF/GenericPolytopeQuadrature.h"

#include "H1.h"
#include "H1Element.h"

namespace Rodin::Variational
{
  template <
    size_t K,
    class LHSDerived, class RHSDerived,
    class Scalar, class Mesh>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Grad<ShapeFunction<LHSDerived, H1<K, Scalar, Mesh>, TrialSpace>>,
        H1<K, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        Grad<ShapeFunction<RHSDerived, H1<K, Scalar, Mesh>, TestSpace>>,
        H1<K, Scalar, Mesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<
            ShapeFunctionBase<
              Grad<ShapeFunction<LHSDerived, H1<K, Scalar, Mesh>, TrialSpace>>,
              H1<K, Scalar, Mesh>, TrialSpace>,
            ShapeFunctionBase<
              Grad<ShapeFunction<RHSDerived, H1<K, Scalar, Mesh>, TestSpace>>,
              H1<K, Scalar, Mesh>, TestSpace>>>::ScalarType>
  {
    public:
      using FESType = H1<K, Scalar, Mesh>;

      using LHSType =
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, FESType, TrialSpace>>, FESType, TrialSpace>;

      using RHSType =
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, FESType, TestSpace>>,  FESType, TestSpace>;

      using IntegrandType = Dot<LHSType, RHSType>;

      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point),
          m_dim(0),
          m_trCount(0),
          m_teCount(0)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point),
          m_dim(0),
          m_trCount(0),
          m_teCount(0)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_ps(std::move(other.m_ps)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_dim(std::exchange(other.m_dim, 0)),
          m_trCount(std::exchange(other.m_trCount, 0)),
          m_teCount(std::exchange(other.m_teCount, 0)),
          m_trRefGrad(std::move(other.m_trRefGrad)),
          m_teRefGrad(std::move(other.m_teRefGrad)),
          m_mat(std::move(other.m_mat))
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

        const auto& integrand = *m_integrand;
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();

        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes  = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes.getFiniteElement(d, idx);

        // For grad-grad, the integrand is degree 2*(K-1) in reference coordinates.
        // Keep the "get(order, geometry)" strategy; just choose the right order.
        // If your FE order is K, trialfe.getOrder()==K etc.
        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe .getOrder();
        const size_t k    = std::max(k_tr, k_te);
        const size_t order = (k > 0) ? 2 * (k - 1) : 0;

        const auto geometry = polytope.getGeometry();

        const bool recompute_qf = (!m_set || m_order != order || m_geometry != geometry);

        if (recompute_qf)
        {
          m_set = true;
          m_order = order;
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

        m_mat.resize(nte, ntr);
        m_mat.setZero();

        // Resize ref-grad caches if needed (qp-major)
        const bool recompute_sizes =
          (m_dim != d) || (m_trCount != ntr) || (m_teCount != nte) || (m_qps != m_ps.size());

        if (recompute_sizes)
        {
          m_dim = d;
          m_qps = m_ps.size();
          m_trCount = ntr;
          m_teCount = nte;

          m_trRefGrad.resize(m_qps * m_trCount);
          m_teRefGrad.resize(m_qps * m_teCount);
          for (auto& v : m_trRefGrad) v.resize(d);
          for (auto& v : m_teRefGrad) v.resize(d);
        }

        // Tabulations (reference gradients) tied to this QF instance
        const auto& trTab = trialfe.getTabulation(*m_qf);
        const auto& teTab = testfe .getTabulation(*m_qf);

        static thread_local Math::SpatialVector<ScalarType> gtr_phys;
        static thread_local Math::SpatialVector<ScalarType> gte_phys;
        gtr_phys.resize(d);
        gte_phys.resize(d);

        // Assemble local matrix
        for (size_t qp = 0; qp < m_qps; ++qp)
        {
          const auto& p = m_ps[qp];
          const auto JinvT = p.getJacobianInverse().transpose();

          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          // cache reference gradients for this qp (optional; keeps API uniform)
          for (size_t a = 0; a < m_trCount; ++a)
          {
            const auto gref = trTab.getGradient(qp, a); // span/ptr length d
            auto& dst = m_trRefGrad[qp * m_trCount + a];
            for (size_t j = 0; j < d; ++j) dst(j) = gref[j];
          }
          for (size_t b = 0; b < m_teCount; ++b)
          {
            const auto gref = teTab.getGradient(qp, b);
            auto& dst = m_teRefGrad[qp * m_teCount + b];
            for (size_t j = 0; j < d; ++j) dst(j) = gref[j];
          }

          // A(b,a) += wdet * (J^{-T} ghat_a)·(J^{-T} ghat_b)
          for (size_t b = 0; b < m_teCount; ++b)
          {
            gte_phys = JinvT * m_teRefGrad[qp * m_teCount + b];

            for (size_t a = 0; a < m_trCount; ++a)
            {
              gtr_phys = JinvT * m_trRefGrad[qp * m_trCount + a];
              m_mat(b, a) += wdet * Math::dot(gtr_phys, gte_phys);
            }
          }
        }

        return *this;
      }

      inline ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_mat(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::vector<Geometry::Point> m_ps;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      // sizes
      size_t m_dim;
      size_t m_qps;
      size_t m_trCount;
      size_t m_teCount;

      // reference gradients (qp-major)
      std::vector<Math::SpatialVector<ScalarType>> m_trRefGrad;
      std::vector<Math::SpatialVector<ScalarType>> m_teRefGrad;

      // local matrix (rows=test, cols=trial)
      Math::Matrix<ScalarType> m_mat;
    };

  // CTAD helper
  template <size_t K, class LHSDerived, class RHSDerived, class Scalar, class Mesh>
  QuadratureRule(
    const Dot<
      ShapeFunctionBase<
        Grad<ShapeFunction<LHSDerived, H1<K, Scalar, Mesh>, TrialSpace>>,
        H1<K, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        Grad<ShapeFunction<RHSDerived, H1<K, Scalar, Mesh>, TestSpace>>,
        H1<K, Scalar, Mesh>, TestSpace>>&)
  -> QuadratureRule<
       Dot<
         ShapeFunctionBase<
           Grad<ShapeFunction<LHSDerived, H1<K, Scalar, Mesh>, TrialSpace>>,
           H1<K, Scalar, Mesh>, TrialSpace>,
         ShapeFunctionBase<
           Grad<ShapeFunction<RHSDerived, H1<K, Scalar, Mesh>, TestSpace>>,
           H1<K, Scalar, Mesh>, TestSpace>>>;
}

#endif

