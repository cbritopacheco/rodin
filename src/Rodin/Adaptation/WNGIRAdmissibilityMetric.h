/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_WNGIRADMISSIBILITYMETRIC_H
#define RODIN_ADAPTATION_WNGIRADMISSIBILITYMETRIC_H

#include "WNGIRDetail.h"
#include "WNGIRParameters.h"

namespace Rodin::Adaptation::Detail
{
    template <class TrialFunction, class TestFunction, class Displacement>
    class WNGIRAdmissibilityMetric final
      : public Variational::LocalBilinearFormIntegratorBase<
          typename TrialFunction::ScalarType>
    {
      public:
        using ScalarType = typename TrialFunction::ScalarType;
        using Parent =
          Variational::LocalBilinearFormIntegratorBase<ScalarType>;

        WNGIRAdmissibilityMetric(
            const TrialFunction& du,
            const TestFunction& v,
            const Displacement& current,
            const WNGIRParameters& parameters)
          : Parent(du.getLeaf(), v.getLeaf()),
            m_du(du),
            m_v(v),
            m_current(current),
            m_parameters(parameters)
        {}

        WNGIRAdmissibilityMetric(
            const WNGIRAdmissibilityMetric&) = default;

        const Geometry::Polytope& getPolytope() const final override
        {
          assert(m_polytope);
          return *m_polytope;
        }

        WNGIRAdmissibilityMetric& setPolytope(
            const Geometry::Polytope& polytope) final override
        {
          m_polytope = &polytope;

          const std::size_t dim = polytope.getDimension();
          const Real d = static_cast<Real>(dim);
          const auto geometry = polytope.getGeometry();
          const auto idx = polytope.getIndex();

          const auto& trialFES = m_du.getFiniteElementSpace();
          const auto& testFES = m_v.getFiniteElementSpace();
          const auto& trialFE = trialFES.getFiniteElement(dim, idx);
          const auto& testFE = testFES.getFiniteElement(dim, idx);
          const auto& params = m_parameters.get();
          const std::size_t qOrder = params.quadratureOrder > 0
            ? params.quadratureOrder
            : std::max<std::size_t>(2, 2 * trialFE.getOrder());
          const auto& qf =
            QF::PolytopeQuadratureFormula::get(qOrder, geometry);
          const auto& quad = polytope.getQuadrature(qf);

          const std::size_t ntr = trialFE.getCount();
          const std::size_t nte = testFE.getCount();
          m_matrix.resize(
              static_cast<Eigen::Index>(nte),
              static_cast<Eigen::Index>(ntr));
          m_matrix.setZero();

          for (std::size_t qp = 0; qp < quad.getSize(); ++qp)
          {
            const auto& pt = quad.getPoint(qp);
            const Variational::IntegrationPoint ip(pt, &qf, qp);
            const Real w = qf.getWeight(qp) * pt.getDistortion();
            const auto F =
              deformationGradient(m_current.get(), polytope, ip, dim);
            const Real jK = F.determinant();
            // Near-singular: cofactor via inverse is unreliable; skip.
            if (std::abs(jK) < Real(1e-14))
              continue;
            // Shape terms (Q_rel and its barriers/hinge) need j > 0; the
            // size hinge needs only j and a_j, so it runs for ALL j —
            // including inverted cells, which it pulls back to validity.
            const bool shapeOK = jK > Real(0);
            Real frob2 = 0;
            Real qK = 0;
            Real sJ0 = 0, sQ0 = 0;
            bool jActive = false, qActive = false, qualActive = false;
            if (shapeOK)
            {
              frob2 = F.squaredNorm();
              qK = frob2 / (d * std::pow(jK, Real(2) / d));
              sJ0 = jK - params.jSafe;
              sQ0 = params.qMax - qK;
              jActive = params.gammaJ > Real(0)
                && sJ0 > Real(0) && sJ0 < params.s0J;
              qActive = params.gammaQ > Real(0)
                && sQ0 > Real(0) && sQ0 < params.s0Q;
              qualActive = params.gammaQual > Real(0) && qK > params.qStar;
            }
            const bool sizeActive =
              params.gammaSize > Real(0) && jK < params.jStar;
            if (!jActive && !qActive && !qualActive && !sizeActive)
              continue;

            const auto Jinv = pt.getJacobianInverse();
            const auto& rc = pt.getReferenceCoordinates();
            const auto FinvT = F.inverse().transpose();
            const auto dQdF = shapeOK
              ? Math::SpatialMatrix<Real>(
                  (Real(2) / d) * std::pow(jK, -Real(2) / d)
                  * (F - (frob2 / d) * FinvT))
              : makeZeroMatrix(dim);

            std::vector<Real> aJTrial(ntr), aQTrial(ntr);
            std::vector<Real> aJTest(nte), aQTest(nte);
            auto fillActions =
              [&](const auto& fe, std::vector<Real>& aJ, std::vector<Real>& aQ)
            {
              for (std::size_t l = 0; l < fe.getCount(); ++l)
              {
                const auto jp = physicalJacobian(fe, l, rc, Jinv, dim);
                Real accJ = 0;
                Real accQ = 0;
                for (std::size_t r = 0; r < dim; ++r)
                  for (std::size_t c = 0; c < dim; ++c)
                  {
                    const auto rr = static_cast<Eigen::Index>(r);
                    const auto cc = static_cast<Eigen::Index>(c);
                    accJ += jK * FinvT(rr, cc) * jp(rr, cc);
                    accQ += dQdF(rr, cc) * jp(rr, cc);
                  }
                aJ[l] = accJ;
                aQ[l] = accQ;
              }
            };
            fillActions(trialFE, aJTrial, aQTrial);
            fillActions(testFE, aJTest, aQTest);

            if (jActive)
            {
              const Real bpp = Real(1) / (sJ0 * sJ0);
              for (std::size_t te = 0; te < nte; ++te)
                for (std::size_t tr = 0; tr < ntr; ++tr)
                  m_matrix(
                      static_cast<Eigen::Index>(te),
                      static_cast<Eigen::Index>(tr))
                    += w * params.gammaJ * bpp
                    * aJTrial[tr] * aJTest[te];
            }
            if (qActive)
            {
              const Real bpp = Real(1) / (sQ0 * sQ0);
              for (std::size_t te = 0; te < nte; ++te)
                for (std::size_t tr = 0; tr < ntr; ++tr)
                  m_matrix(
                      static_cast<Eigen::Index>(te),
                      static_cast<Eigen::Index>(tr))
                    += w * params.gammaQ * bpp
                    * aQTrial[tr] * aQTest[te];
            }
            if (qualActive)
            {
              // Gauss–Newton Hessian of the quality hinge: ρ''(Q)=1.
              for (std::size_t te = 0; te < nte; ++te)
                for (std::size_t tr = 0; tr < ntr; ++tr)
                  m_matrix(
                      static_cast<Eigen::Index>(te),
                      static_cast<Eigen::Index>(tr))
                    += w * params.gammaQual
                    * aQTrial[tr] * aQTest[te];
            }
            if (sizeActive)
            {
              // GN Hessian of the size hinge: ρ''(j)=1 for j < jStar.
              for (std::size_t te = 0; te < nte; ++te)
                for (std::size_t tr = 0; tr < ntr; ++tr)
                  m_matrix(
                      static_cast<Eigen::Index>(te),
                      static_cast<Eigen::Index>(tr))
                    += w * params.gammaSize
                    * aJTrial[tr] * aJTest[te];
            }
          }
          return *this;
        }

        ScalarType integrate(std::size_t tr, std::size_t te) final override
        {
          return m_matrix(
              static_cast<Eigen::Index>(te),
              static_cast<Eigen::Index>(tr));
        }

        Geometry::Region getRegion() const final override
        {
          return Geometry::Region::Cells;
        }

        WNGIRAdmissibilityMetric* copy() const noexcept final override
        {
          return new WNGIRAdmissibilityMetric(*this);
        }

      private:
        TrialFunction m_du;
        TestFunction m_v;
        std::reference_wrapper<const Displacement> m_current;
        std::reference_wrapper<const WNGIRParameters> m_parameters;
        const Geometry::Polytope* m_polytope = nullptr;
        Math::Matrix<Real> m_matrix;
    };
}

#endif
