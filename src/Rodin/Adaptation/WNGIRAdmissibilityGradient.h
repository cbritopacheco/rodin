/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file WNGIRAdmissibilityGradient.h
 * @brief Linear-form integrator assembling the admissibility gradient
 * contribution to the WNGIR step system.
 */
#ifndef RODIN_ADAPTATION_WNGIRADMISSIBILITYGRADIENT_H
#define RODIN_ADAPTATION_WNGIRADMISSIBILITYGRADIENT_H

#include "WNGIRDetail.h"
#include "WNGIRParameters.h"

namespace Rodin::Adaptation::Detail
{
    template <class TestFunction, class Displacement>
    class WNGIRAdmissibilityGradient final
      : public Variational::LinearFormIntegratorBase<
          typename TestFunction::ScalarType>
    {
      public:
        /// @brief Scalar value type.
        using ScalarType = typename TestFunction::ScalarType;
        /// @brief Parent class type.
        using Parent = Variational::LinearFormIntegratorBase<ScalarType>;

        WNGIRAdmissibilityGradient(
            const TestFunction& v,
            const Displacement& current,
            const WNGIRParameters& parameters)
          : Parent(v.getLeaf()),
            m_v(v),
            m_current(current),
            m_parameters(parameters)
        {}

        WNGIRAdmissibilityGradient(
            const WNGIRAdmissibilityGradient&) = default;

        const Geometry::Polytope& getPolytope() const final override
        {
          assert(m_polytope);
          return *m_polytope;
        }

        WNGIRAdmissibilityGradient& setPolytope(
            const Geometry::Polytope& polytope) final override
        {
          m_polytope = &polytope;

          const std::size_t dim = polytope.getDimension();
          const Real d = static_cast<Real>(dim);
          const auto geometry = polytope.getGeometry();
          const auto idx = polytope.getIndex();

          const auto& testFES = m_v.getFiniteElementSpace();
          const auto& testFE = testFES.getFiniteElement(dim, idx);
          const std::size_t qOrder = m_parameters.get().quadratureOrder > 0
            ? m_parameters.get().quadratureOrder
            : std::max<std::size_t>(2, 2 * testFE.getOrder());
          const auto& qf =
            QF::PolytopeQuadratureFormula::get(qOrder, geometry);
          const auto& quad = polytope.getQuadrature(qf);

          const std::size_t nte = testFE.getCount();
          m_vector.resize(static_cast<Eigen::Index>(nte));
          m_vector.setZero();

          for (std::size_t qp = 0; qp < quad.getSize(); ++qp)
          {
            const auto& pt = quad.getPoint(qp);
            const Variational::IntegrationPoint ip(pt, &qf, qp);
            const Real w = qf.getWeight(qp) * pt.getDistortion();
            const auto F =
              deformationGradient(m_current.get(), polytope, ip, dim);
            const Real jK = F.determinant();
            if (std::abs(jK) < Real(1e-14))
              continue;
            const auto& params = m_parameters.get();
            // Shape/Q terms need j > 0; size terms run for ALL j.
            const bool shapeOK = jK > Real(0);
            Real frob2 = 0, qK = 0, sJ0 = 0, sQ0 = 0;
            bool jBarrierActive = false, qBarrierActive = false;
            bool qualHingeActive = false;
            if (shapeOK)
            {
              frob2 = F.squaredNorm();
              qK = frob2 / (d * std::pow(jK, Real(2) / d));
              sJ0 = jK - params.jSafe;
              sQ0 = params.qMax - qK;
              jBarrierActive =
                (params.includeAdmissibilityGradient
                 || params.includeQualityGradient)
                && params.gammaJ > Real(0)
                && sJ0 > Real(0) && sJ0 < params.s0J;
              qBarrierActive = params.includeAdmissibilityGradient
                && params.gammaQ > Real(0)
                && sQ0 > Real(0) && sQ0 < params.s0Q;
              qualHingeActive =
                params.includeQualityGradient
                && params.gammaQual > Real(0) && qK > params.qStar;
            }
            const bool sizeHingeActive =
              params.includeQualityGradient
              && params.gammaSize > Real(0) && jK < params.jStar;
            if (!jBarrierActive && !qBarrierActive
                && !qualHingeActive && !sizeHingeActive)
              continue;

            const auto Jinv = pt.getJacobianInverse();
            const auto& rc = pt.getReferenceCoordinates();
            const auto FinvT = F.inverse().transpose();
            const auto dQdF = shapeOK
              ? Math::SpatialMatrix<Real>(
                  (Real(2) / d) * std::pow(jK, -Real(2) / d)
                  * (F - (frob2 / d) * FinvT))
              : makeZeroMatrix(dim);

            for (std::size_t te = 0; te < nte; ++te)
            {
              const auto jp = physicalJacobian(testFE, te, rc, Jinv, dim);
              Real aJ = 0;
              Real aQ = 0;
              for (std::size_t r = 0; r < dim; ++r)
                for (std::size_t c = 0; c < dim; ++c)
                {
                  const auto rr = static_cast<Eigen::Index>(r);
                  const auto cc = static_cast<Eigen::Index>(c);
                  aJ += jK * FinvT(rr, cc) * jp(rr, cc);
                  aQ += dQdF(rr, cc) * jp(rr, cc);
                }

              Real val = 0;
              if (jBarrierActive)
              {
                const Real bp = -Real(1) / sJ0 + Real(1) / params.s0J;
                val += -params.gammaJ * bp * aJ;
              }
              if (qBarrierActive)
              {
                const Real bp = -Real(1) / sQ0 + Real(1) / params.s0Q;
                val += params.gammaQ * bp * aQ;
              }
              if (qualHingeActive)
              {
                const Real excess = qK - params.qStar;
                val += -params.gammaQual * excess * aQ;
              }
              if (sizeHingeActive)
              {
                const Real deficit = params.jStar - jK;
                val += params.gammaSize * deficit * aJ;
              }
              m_vector(static_cast<Eigen::Index>(te)) += w * val;
            }
          }
          return *this;
        }

        ScalarType integrate(std::size_t local) final override
        {
          return m_vector(static_cast<Eigen::Index>(local));
        }

        Geometry::Region getRegion() const final override
        {
          return Geometry::Region::Cells;
        }

        WNGIRAdmissibilityGradient* copy() const noexcept final override
        {
          return new WNGIRAdmissibilityGradient(*this);
        }

      private:
        TestFunction m_v;
        std::reference_wrapper<const Displacement> m_current;
        std::reference_wrapper<const WNGIRParameters> m_parameters;
        const Geometry::Polytope* m_polytope = nullptr;
        Math::Vector<Real> m_vector;
    };
}

#endif
