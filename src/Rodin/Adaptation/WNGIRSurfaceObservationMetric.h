/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_WNGIRSURFACEOBSERVATIONMETRIC_H
#define RODIN_ADAPTATION_WNGIRSURFACEOBSERVATIONMETRIC_H

#include "WNGIRDetail.h"
#include "WNGIRParameters.h"

namespace Rodin::Adaptation::Detail
{
    /**
     * @brief Bilinear-form integrator for the WNGIR surface observation metric.
     */
  template <class PhiDerived, class GradDerived, class TrialFunction, class TestFunction,
    class Displacement>
  class WNGIRSurfaceObservationMetric final
    : public Variational::LocalBilinearFormIntegratorBase<
        typename TrialFunction::ScalarType>
  {
    public:
      using ScalarType = typename TrialFunction::ScalarType;
      using Parent = Variational::LocalBilinearFormIntegratorBase<ScalarType>;
        /// @brief Level-set function type.
      using PhiType = Variational::RealFunctionBase<PhiDerived>;
        /// @brief Level-set gradient function type.
      using GradType = Variational::VectorFunctionBase<Real, GradDerived>;

        /// @brief Constructs the surface observation metric integrator.
      WNGIRSurfaceObservationMetric(const PhiType& phi, const GradType& grad,
        const TrialFunction& du, const TestFunction& v, const Displacement& current,
        const WNGIRParameters& parameters, Real sigma2)
        : Parent(du.getLeaf(), v.getLeaf()),
          m_phi(phi.copy()),
          m_grad(grad.copy()),
          m_du(std::cref(du)),
          m_v(std::cref(v)),
          m_current(current),
          m_parameters(parameters),
          m_sigma2(sigma2)
      {}

        /// @brief Copy constructor.
      WNGIRSurfaceObservationMetric(const WNGIRSurfaceObservationMetric& other)
        : Parent(other),
          m_phi(other.m_phi ? other.m_phi->copy() : nullptr),
          m_grad(other.m_grad ? other.m_grad->copy() : nullptr),
          m_du(other.m_du),
          m_v(other.m_v),
          m_current(other.m_current),
          m_parameters(other.m_parameters),
          m_sigma2(other.m_sigma2),
          m_polytope(other.m_polytope)
      {}

        /// @brief Returns the current polytope.
      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

        /// @brief Sets the current face and assembles the local matrix.
      WNGIRSurfaceObservationMetric& setPolytope(
        const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const std::size_t faceDim = polytope.getDimension();
        const Index faceIdx = polytope.getIndex();
        const auto& trialFES = m_du.get().getFiniteElementSpace();
        const auto& testFES = m_v.get().getFiniteElementSpace();
        const auto& trialFE = trialFES.getFiniteElement(faceDim, faceIdx);
        const auto& testFE = testFES.getFiniteElement(faceDim, faceIdx);
        const auto& params = m_parameters.get();
        const std::size_t qOrder = params.quadratureOrder > 0
          ? params.quadratureOrder
          : std::max<std::size_t>(2, 2 * std::max(trialFE.getOrder(), testFE.getOrder()));
        const auto& qf =
          QF::PolytopeQuadratureFormula::get(qOrder, polytope.getGeometry());
        const auto& quad = polytope.getQuadrature(qf);

        const std::size_t ntr = trialFE.getCount();
        const std::size_t nte = testFE.getCount();
        m_matrix.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_matrix.setZero();

        constexpr Real epsG = Real(1e-12);
        for (std::size_t qp = 0; qp < quad.getSize(); ++qp)
        {
          const auto& pt = quad.getPoint(qp);
          const Variational::IntegrationPoint ip(pt, &qf, qp);
          const Real w = qf.getWeight(qp) * pt.getDistortion();
          const auto uq = m_current.get().getValue(ip);
          Math::SpatialVector<Real> displacement(polytope.getMesh().getDimension());
          displacement.setZero();
          for (std::size_t r = 0; r < static_cast<std::size_t>(displacement.size()); ++r)
            displacement(static_cast<Eigen::Index>(r)) = uq(r);

          const auto moved = makeTranslatedPoint(pt,
            pt.getPhysicalCoordinates() + displacement, params.pointLocationTolerance);
          const Real r = m_phi->getValue(moved);
          const auto g = m_grad->getValue(moved);

          const auto& rc = pt.getReferenceCoordinates();
          if (params.observationMetric == WNGIRObservationMetric::Isotropic)
          {
            const Real obsWeight = g.dot(g) + epsG +
              (params.residualStabilizedObservationMetric ? (r * r) / m_sigma2 : Real(0));
            for (std::size_t te = 0; te < nte; ++te)
            {
              const auto testValue = testFE.getBasis(te)(rc);
              for (std::size_t tr = 0; tr < ntr; ++tr)
              {
                const auto trialValue = trialFE.getBasis(tr)(rc);
                m_matrix(static_cast<Eigen::Index>(te), static_cast<Eigen::Index>(tr)) +=
                  w * params.gammaObs * obsWeight * trialValue.dot(testValue);
              }
            }
          }
          else
          {
            const Real obsWeight =
              params.observationMetric == WNGIRObservationMetric::RankOneIRLS ||
                params.observationMetric == WNGIRObservationMetric::HybridRankOneIRLS
              ? std::exp(-r * r / m_sigma2)
              : Real(1);
            const Real gradientNorm2 = g.dot(g);
            for (std::size_t te = 0; te < nte; ++te)
            {
              const auto testValue = testFE.getBasis(te)(rc);
              const Real testNormal = g.dot(testValue);
              for (std::size_t tr = 0; tr < ntr; ++tr)
              {
                const auto trialValue = trialFE.getBasis(tr)(rc);
                const Real trialNormal = g.dot(trialValue);
                Real metricValue = trialNormal * testNormal;
                if (params.observationMetric == WNGIRObservationMetric::HybridRankOneIRLS)
                {
                  const Real tangentialProduct = trialValue.dot(testValue) -
                    trialNormal * testNormal / (gradientNorm2 + epsG);
                  metricValue +=
                    params.observationTangentialFloor * gradientNorm2 * tangentialProduct;
                }
                m_matrix(static_cast<Eigen::Index>(te), static_cast<Eigen::Index>(tr)) +=
                  w * params.gammaObs * obsWeight * metricValue;
              }
            }
          }
        }
        return *this;
      }

        /// @brief Returns an entry of the current face matrix.
      ScalarType integrate(std::size_t tr, std::size_t te) final override
      {
        return m_matrix(static_cast<Eigen::Index>(te), static_cast<Eigen::Index>(tr));
      }

        /// @brief Returns the integration region.
      Geometry::Region getRegion() const final override
      {
        return Geometry::Region::Faces;
      }

      WNGIRSurfaceObservationMetric* copy() const noexcept final override
      {
        return new WNGIRSurfaceObservationMetric(*this);
      }

    private:
      std::unique_ptr<PhiType> m_phi;
      std::unique_ptr<GradType> m_grad;
      std::reference_wrapper<const TrialFunction> m_du;
      std::reference_wrapper<const TestFunction> m_v;
      std::reference_wrapper<const Displacement> m_current;
      std::reference_wrapper<const WNGIRParameters> m_parameters;
      Real m_sigma2;
      const Geometry::Polytope* m_polytope = nullptr;
      Math::Matrix<Real> m_matrix;
  };
}

#endif
