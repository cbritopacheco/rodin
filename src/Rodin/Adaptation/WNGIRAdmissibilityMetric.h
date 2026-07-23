/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_WNGIRADMISSIBILITYMETRIC_H
#define RODIN_ADAPTATION_WNGIRADMISSIBILITYMETRIC_H

#include "CellDeformation.h"
#include "WNGIRConstraintState.h"
#include "WNGIRParameters.h"
#include "WNGIRValidationWeights.h"

namespace Rodin::Adaptation::Detail
{
    /**
     * @brief Bilinear-form integrator for WNGIR admissibility and quality metrics.
     */
  template <class TrialFunction, class TestFunction, class Displacement>
  class WNGIRAdmissibilityMetric final
    : public Variational::LocalBilinearFormIntegratorBase<
        typename TrialFunction::ScalarType>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = typename TrialFunction::ScalarType;
      /// @brief Parent class type.
      using Parent = Variational::LocalBilinearFormIntegratorBase<ScalarType>;

        /// @brief Constructs the admissibility-metric integrator.
      WNGIRAdmissibilityMetric(const TrialFunction& du, const TestFunction& v,
        const Displacement& current, const WNGIRParameters& parameters,
        const Displacement* predictor = nullptr, Real constraintMultiplier = Real(1))
        : Parent(du.getLeaf(), v.getLeaf()),
          m_du(std::cref(du)),
          m_v(std::cref(v)),
          m_current(current),
          m_parameters(parameters),
          m_predictor(predictor),
          m_constraintMultiplier(constraintMultiplier)
      {}

        /// @brief Copy constructor.
      WNGIRAdmissibilityMetric(const WNGIRAdmissibilityMetric&) = default;

        /// @brief Returns the current polytope.
      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

        /// @brief Sets the current polytope and assembles the local matrix.
      WNGIRAdmissibilityMetric& setPolytope(
        const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const std::size_t dim = polytope.getDimension();
        const Real d = static_cast<Real>(dim);
        const auto geometry = polytope.getGeometry();
        const auto idx = polytope.getIndex();

        const auto& trialFES = m_du.get().getFiniteElementSpace();
        const auto& testFES = m_v.get().getFiniteElementSpace();
        const auto& trialFE = trialFES.getFiniteElement(dim, idx);
        const auto& testFE = testFES.getFiniteElement(dim, idx);
        const auto& params = m_parameters.get();
        const std::size_t qOrder = params.quadratureOrder > 0
          ? params.quadratureOrder
          : std::max<std::size_t>(2, 2 * trialFE.getOrder());
        const auto& qf = QF::PolytopeQuadratureFormula::get(qOrder, geometry);
        const auto& quad = polytope.getQuadrature(qf);
        const WNGIRValidationWeights validationWeights(qf);

        const std::size_t ntr = trialFE.getCount();
        const std::size_t nte = testFE.getCount();
        m_matrix.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_matrix.setZero();
        m_aJTrial.resize(ntr);
        m_aQTrial.resize(ntr);
        m_aJTest.resize(nte);
        m_aQTest.resize(nte);

        auto trialJacobian = Variational::Jacobian(m_du.get());
        auto testJacobian = Variational::Jacobian(m_v.get());
        auto currentJacobian = Variational::Jacobian(m_current.get());
        auto predictorJacobian =
          Variational::Jacobian(m_predictor ? *m_predictor : m_current.get());
        const bool sameSpace = trialFES == testFES;

        CellDeformation def(dim);
        for (std::size_t qp = 0; qp < quad.getSize(); ++qp)
        {
          const auto& pt = quad.getPoint(qp);
          const Variational::IntegrationPoint ip(pt, &qf, qp);
          const Real w = validationWeights.getWeight(qp) * pt.getDistortion();
          def.setDisplacementGradient(currentJacobian.getValue(ip));
            // Near-singular: cofactor via inverse is unreliable; skip.
          if (!def.isInvertible())
            continue;
          Optional<Math::SpatialMatrix<Real>> predictorGradient;
          if (m_predictor)
            predictorGradient = predictorJacobian.getValue(ip);
          const WNGIRConstraintState state(
            def, params, predictorGradient ? &*predictorGradient : nullptr);
          const Real jWeight = state.getJacobianWeight();
          const Real qWeight = state.getDistortionWeight();
          const Real qualWeight = state.getQualityWeight();
          const Real sizeWeight = state.getSizeWeight();
          if (jWeight <= Real(0) && qWeight <= Real(0) && qualWeight <= Real(0) &&
            sizeWeight <= Real(0))
            continue;

          const bool shapeOK = def.isAdmissible();

          trialJacobian.setIntegrationPoint(ip);
          for (std::size_t l = 0; l < ntr; ++l)
          {
            const auto jp = trialJacobian.getBasis(l);
            m_aJTrial[l] = def.getJacobianAction(jp);
            m_aQTrial[l] = shapeOK ? def.getRelativeDistortionAction(jp) : Real(0);
          }
          if (sameSpace)
          {
            m_aJTest = m_aJTrial;
            m_aQTest = m_aQTrial;
          }
          else
          {
            testJacobian.setIntegrationPoint(ip);
            for (std::size_t l = 0; l < nte; ++l)
            {
              const auto jp = testJacobian.getBasis(l);
              m_aJTest[l] = def.getJacobianAction(jp);
              m_aQTest[l] = shapeOK ? def.getRelativeDistortionAction(jp) : Real(0);
            }
          }

          if (jWeight > Real(0))
          {
            for (std::size_t te = 0; te < nte; ++te)
              for (std::size_t tr = 0; tr < ntr; ++tr)
                m_matrix(static_cast<Eigen::Index>(te), static_cast<Eigen::Index>(tr)) +=
                  w * m_constraintMultiplier * params.gammaJ * jWeight * m_aJTrial[tr] *
                  m_aJTest[te];
          }
          if (qWeight > Real(0))
          {
            for (std::size_t te = 0; te < nte; ++te)
              for (std::size_t tr = 0; tr < ntr; ++tr)
                m_matrix(static_cast<Eigen::Index>(te), static_cast<Eigen::Index>(tr)) +=
                  w * m_constraintMultiplier * params.gammaQ * qWeight * m_aQTrial[tr] *
                  m_aQTest[te];
          }
          if (qualWeight > Real(0))
          {
              // Gauss–Newton Hessian of the quality hinge: ρ''(Q)=1.
            for (std::size_t te = 0; te < nte; ++te)
              for (std::size_t tr = 0; tr < ntr; ++tr)
                m_matrix(static_cast<Eigen::Index>(te), static_cast<Eigen::Index>(tr)) +=
                  w * params.gammaQual * qualWeight * m_aQTrial[tr] * m_aQTest[te];
          }
          if (sizeWeight > Real(0))
          {
              // GN Hessian of the size hinge: ρ''(j)=1 for j < jStar.
            for (std::size_t te = 0; te < nte; ++te)
              for (std::size_t tr = 0; tr < ntr; ++tr)
                m_matrix(static_cast<Eigen::Index>(te), static_cast<Eigen::Index>(tr)) +=
                  w * params.gammaSize * sizeWeight * m_aJTrial[tr] * m_aJTest[te];
          }
        }
        return *this;
      }

        /// @brief Returns an entry of the current element matrix.
      ScalarType integrate(std::size_t tr, std::size_t te) final override
      {
        return m_matrix(static_cast<Eigen::Index>(te), static_cast<Eigen::Index>(tr));
      }

        /// @brief Returns the integration region.
      Geometry::Region getRegion() const final override
      {
        return Geometry::Region::Cells;
      }

      /// @brief Clones this integrator.
      WNGIRAdmissibilityMetric* copy() const noexcept final override
      {
        return new WNGIRAdmissibilityMetric(*this);
      }

    private:
      std::reference_wrapper<const TrialFunction> m_du;
      std::reference_wrapper<const TestFunction> m_v;
      std::reference_wrapper<const Displacement> m_current;
      std::reference_wrapper<const WNGIRParameters> m_parameters;
      const Displacement* m_predictor;
      Real m_constraintMultiplier;
      const Geometry::Polytope* m_polytope = nullptr;
      std::vector<Real> m_aJTrial;
      std::vector<Real> m_aQTrial;
      std::vector<Real> m_aJTest;
      std::vector<Real> m_aQTest;
      Math::Matrix<Real> m_matrix;
  };
}

#endif
