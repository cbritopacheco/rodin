/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_WNGIRPRIMALBARRIERMETRIC_H
#define RODIN_ADAPTATION_WNGIRPRIMALBARRIERMETRIC_H

#include "CellDeformation.h"
#include "WNGIRPrimalBarrierState.h"
#include "WNGIRValidationWeights.h"

namespace Rodin::Adaptation::Detail
{
  /// @brief Newton Hessian of the primal barrier for the linearized QP.
  template <class TrialFunction, class TestFunction, class Displacement>
  class WNGIRPrimalBarrierMetric final
    : public Variational::LocalBilinearFormIntegratorBase<
        typename TrialFunction::ScalarType>
  {
    public:
      using ScalarType = typename TrialFunction::ScalarType;
      using Parent = Variational::LocalBilinearFormIntegratorBase<ScalarType>;

      WNGIRPrimalBarrierMetric(const TrialFunction& du, const TestFunction& z,
        const Displacement& current, const Displacement& inner,
        const WNGIRParameters& parameters, Real barrierCoefficient)
        : Parent(du.getLeaf(), z.getLeaf()),
          m_du(du), m_z(z), m_current(current), m_inner(inner),
          m_parameters(parameters), m_barrierCoefficient(barrierCoefficient)
      {}

      WNGIRPrimalBarrierMetric(const WNGIRPrimalBarrierMetric&) = default;

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      WNGIRPrimalBarrierMetric& setPolytope(
        const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;
        const std::size_t dim = polytope.getDimension();
        const Index index = polytope.getIndex();
        const auto& trialFES = m_du.get().getFiniteElementSpace();
        const auto& testFES = m_z.get().getFiniteElementSpace();
        const auto& trialFE = trialFES.getFiniteElement(dim, index);
        const auto& testFE = testFES.getFiniteElement(dim, index);
        const auto& parameters = m_parameters.get();
        const std::size_t order = parameters.quadratureOrder > 0
          ? parameters.quadratureOrder
          : std::max<std::size_t>(2, 2 * trialFE.getOrder());
        const auto& qf =
          QF::PolytopeQuadratureFormula::get(order, polytope.getGeometry());
        const auto& quadrature = polytope.getQuadrature(qf);
        const WNGIRValidationWeights validationWeights(qf);
        const std::size_t nTrial = trialFE.getCount();
        const std::size_t nTest = testFE.getCount();

        m_matrix.resize(
          static_cast<Eigen::Index>(nTest), static_cast<Eigen::Index>(nTrial));
        m_matrix.setZero();
        m_jTrial.resize(nTrial);
        m_qTrial.resize(nTrial);
        m_jTest.resize(nTest);
        m_qTest.resize(nTest);

        auto trialJacobian = Variational::Jacobian(m_du.get());
        auto testJacobian = Variational::Jacobian(m_z.get());
        auto currentJacobian = Variational::Jacobian(m_current.get());
        auto innerJacobian = Variational::Jacobian(m_inner.get());
        const bool sameSpace = trialFES == testFES;
        CellDeformation deformation(dim);

        for (std::size_t q = 0; q < quadrature.getSize(); ++q)
        {
          const auto& point = quadrature.getPoint(q);
          const Variational::IntegrationPoint ip(point, &qf, q);
          deformation.setDisplacementGradient(currentJacobian.getValue(ip));
          if (!deformation.isAdmissible())
            continue;
          const WNGIRPrimalBarrierState state(
            deformation, innerJacobian.getValue(ip), parameters,
            m_barrierCoefficient);
          assert(state.isFeasible());
          if (!state.isFeasible())
            continue;

          trialJacobian.setIntegrationPoint(ip);
          for (std::size_t local = 0; local < nTrial; ++local)
          {
            const auto gradient = trialJacobian.getBasis(local);
            m_jTrial[local] = -deformation.getJacobianAction(gradient);
            m_qTrial[local] = deformation.getRelativeDistortionAction(gradient);
          }
          if (sameSpace)
          {
            m_jTest = m_jTrial;
            m_qTest = m_qTrial;
          }
          else
          {
            testJacobian.setIntegrationPoint(ip);
            for (std::size_t local = 0; local < nTest; ++local)
            {
              const auto gradient = testJacobian.getBasis(local);
              m_jTest[local] = -deformation.getJacobianAction(gradient);
              m_qTest[local] = deformation.getRelativeDistortionAction(gradient);
            }
          }

          const Real weight = validationWeights.getWeight(q) * point.getDistortion();
          const Real weightJ = state.getJacobianHessian();
          const Real weightQ = state.getDistortionHessian();
          for (std::size_t test = 0; test < nTest; ++test)
          {
            for (std::size_t trial = 0; trial < nTrial; ++trial)
            {
              m_matrix(static_cast<Eigen::Index>(test),
                static_cast<Eigen::Index>(trial)) += weight *
                (weightJ * m_jTrial[trial] * m_jTest[test] +
                 weightQ * m_qTrial[trial] * m_qTest[test]);
            }
          }
        }
        return *this;
      }

      ScalarType integrate(std::size_t trial, std::size_t test) final override
      {
        return m_matrix(
          static_cast<Eigen::Index>(test), static_cast<Eigen::Index>(trial));
      }

      Geometry::Region getRegion() const final override
      {
        return Geometry::Region::Cells;
      }

      WNGIRPrimalBarrierMetric* copy() const noexcept final override
      {
        return new WNGIRPrimalBarrierMetric(*this);
      }

    private:
      std::reference_wrapper<const TrialFunction> m_du;
      std::reference_wrapper<const TestFunction> m_z;
      std::reference_wrapper<const Displacement> m_current;
      std::reference_wrapper<const Displacement> m_inner;
      std::reference_wrapper<const WNGIRParameters> m_parameters;
      Real m_barrierCoefficient;
      const Geometry::Polytope* m_polytope = nullptr;
      std::vector<Real> m_jTrial;
      std::vector<Real> m_qTrial;
      std::vector<Real> m_jTest;
      std::vector<Real> m_qTest;
      Math::Matrix<Real> m_matrix;
  };
}

#endif
