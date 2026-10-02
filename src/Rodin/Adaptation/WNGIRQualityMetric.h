/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_WNGIRQUALITYMETRIC_H
#define RODIN_ADAPTATION_WNGIRQUALITYMETRIC_H

#include <Eigen/Eigenvalues>
#include "Rodin/Alert.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/Variational/Jacobian.h"
#include "CellDeformation.h"
#include "WNGIRParameters.h"

namespace Rodin::Adaptation::Detail
{
  /**
   * @brief Positive-semidefinite frozen shape curvature for WNGIR increments.
   *
   * Assembles h*kappaS times the Hessian of (d/4)(Q-1).
   * The tensor is projected onto its positive-semidefinite part, retaining its
   * nonnegative eigenvalues. This adds no quality force.
   * It is frozen at the outer displacement and reused by the inner QP.
   */
  template <class TrialFunction, class TestFunction, class Displacement>
  class WNGIRQualityMetric final : public Variational::LocalBilinearFormIntegratorBase<
                                     typename TrialFunction::ScalarType>
  {
    public:
      using ScalarType = typename TrialFunction::ScalarType;
      using Parent = Variational::LocalBilinearFormIntegratorBase<ScalarType>;

      WNGIRQualityMetric(const TrialFunction& trial, const TestFunction& test,
        const Displacement& current, const WNGIRParameters& parameters)
        : Parent(trial.getLeaf(), test.getLeaf()),
          m_trial(trial),
          m_test(test),
          m_current(current),
          m_parameters(parameters)
      {}

      WNGIRQualityMetric(const WNGIRQualityMetric&) = default;

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      WNGIRQualityMetric& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;
        const auto d = polytope.getDimension();
        const auto index = polytope.getIndex();
        const auto& trialFE =
          m_trial.get().getFiniteElementSpace().getFiniteElement(d, index);
        const auto& testFE =
          m_test.get().getFiniteElementSpace().getFiniteElement(d, index);
        const auto& parameters = m_parameters.get();
        const auto order = parameters.quadratureOrder > 0
          ? parameters.quadratureOrder
          : std::max<size_t>(2, 2 * std::max(trialFE.getOrder(), testFE.getOrder()));
        const auto& qf =
          QF::PolytopeQuadratureFormula::get(order, polytope.getGeometry());
        const auto& quadrature = polytope.getQuadrature(qf);
        m_matrix = Math::Matrix<Real>::Zero(testFE.getCount(), trialFE.getCount());
        auto currentJacobian = Variational::Jacobian(m_current.get());
        auto trialJacobian = Variational::Jacobian(m_trial.get());
        auto testJacobian = Variational::Jacobian(m_test.get());
        const Real coefficient = parameters.h * parameters.kappaS;
        for (size_t q = 0; q < quadrature.getSize(); ++q)
        {
          const auto& point = quadrature.getPoint(q);
          const Variational::IntegrationPoint ip(point, &qf, q);
          CellDeformation deformation(d);
          deformation.setDisplacementGradient(currentJacobian.getValue(ip));
          trialJacobian.setIntegrationPoint(ip);
          testJacobian.setIntegrationPoint(ip);
          Math::Matrix<Real> curvature(d * d, d * d);
          for (size_t a = 0; a < d * d; ++a)
            for (size_t b = 0; b < d * d; ++b)
            {
              Math::SpatialMatrix<Real> G(d, d);
              G.setZero();
              Math::SpatialMatrix<Real> H(d, d);
              H.setZero();
              G(a / d, a % d) = Real(1);
              H(b / d, b % d) = Real(1);
              curvature(a, b) =
                Real(d) / Real(4) * deformation.getRelativeDistortionSecondAction(G, H);
            }
          {
            curvature = (Real(0.5) * (curvature + curvature.transpose())).eval();
            Eigen::SelfAdjointEigenSolver<Math::Matrix<Real>> eigen(curvature);
            if (eigen.info() != Eigen::Success || !eigen.eigenvalues().allFinite())
              Alert::Exception() << "WNGIR shape-curvature eigensolve failed."
                                 << Alert::Raise;
            curvature = eigen.eigenvectors() *
              eigen.eigenvalues().cwiseMax(Real(0)).asDiagonal() *
              eigen.eigenvectors().transpose();
          }
          const Real weight = coefficient * qf.getWeight(q) * point.getDistortion();
          std::vector<Math::Vector<Real>> trialImages(trialFE.getCount());
          for (size_t trial = 0; trial < trialFE.getCount(); ++trial)
          {
            const auto& G = trialJacobian.getBasis(trial);
            Math::Vector<Real> g(d * d);
            for (size_t a = 0; a < d * d; ++a)
              g(a) = G(a / d, a % d);
            trialImages[trial] = curvature * g;
          }
          for (size_t test = 0; test < testFE.getCount(); ++test)
          {
            const auto& H = testJacobian.getBasis(test);
            Math::Vector<Real> h(d * d);
            for (size_t a = 0; a < d * d; ++a)
              h(a) = H(a / d, a % d);
            for (size_t trial = 0; trial < trialFE.getCount(); ++trial)
              m_matrix(test, trial) += weight * h.dot(trialImages[trial]);
          }
        }
        return *this;
      }

      ScalarType integrate(size_t trial, size_t test) final override
      {
        return m_matrix(test, trial);
      }

      Geometry::Region getRegion() const final override
      {
        return Geometry::Region::Cells;
      }

      WNGIRQualityMetric* copy() const noexcept final override
      {
        return new WNGIRQualityMetric(*this);
      }

    private:
      std::reference_wrapper<const TrialFunction> m_trial;
      std::reference_wrapper<const TestFunction> m_test;
      std::reference_wrapper<const Displacement> m_current;
      std::reference_wrapper<const WNGIRParameters> m_parameters;
      const Geometry::Polytope* m_polytope = nullptr;
      Math::Matrix<Real> m_matrix;
  };
}

#endif
