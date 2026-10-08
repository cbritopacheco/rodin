/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_WNGIR_HINGEMETRIC_H
#define RODIN_ADAPTATION_WNGIR_HINGEMETRIC_H

#include "../CellDeformation.h"
#include "Hinge.h"

namespace Rodin::Adaptation::WNGIR
{
  /// @brief Hessian of the affine quadratic quality hinges.
  template <class TrialFunction, class TestFunction, class Displacement>
  class HingeMetric final : public Variational::LocalBilinearFormIntegratorBase<
                                   typename TrialFunction::ScalarType>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = typename TrialFunction::ScalarType;
      /// @brief Parent class type.
      using Parent = Variational::LocalBilinearFormIntegratorBase<ScalarType>;

      /**
       * @brief Constructs the affine quadratic hinge tangent.
       * @param du Trial function for the displacement correction.
       * @param z Test function for the displacement correction.
       * @param current Frozen outer displacement.
       * @param inner Current inner displacement increment.
       * @param parameters Quality budgets and hinge activation weights.
       * @param hingeCoefficient Effective hinge coefficient for this outer model.
       */
      HingeMetric(const TrialFunction& du, const TestFunction& z,
        const Displacement& current, const Displacement& inner,
        const Parameters& parameters, Real hingeCoefficient)
        : Parent(du.getLeaf(), z.getLeaf()),
          m_du(du),
          m_z(z),
          m_current(current),
          m_inner(inner),
          m_parameters(parameters),
          m_hingeCoefficient(hingeCoefficient)
      {}

      /**
       * @brief Copy constructor.
       * @param other Integrator to copy, retaining its field references.
       */
      HingeMetric(const HingeMetric& other) = default;

      /**
       * @brief Returns the current polytope.
       * @returns The current polytope.
       */
      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      /**
       * @brief Binds to a polytope and assembles the local tangent.
       * @param polytope Cell to integrate.
       * @returns This integrator after binding and assembly.
       */
      HingeMetric& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;
        const std::size_t dim = polytope.getDimension();
        const Index index = polytope.getIndex();
        const auto& trialFES = m_du.get().getFiniteElementSpace();
        const auto& testFES = m_z.get().getFiniteElementSpace();
        const auto& trialFE = trialFES.getFiniteElement(dim, index);
        const auto& testFE = testFES.getFiniteElement(dim, index);
        const auto& parameters = m_parameters.get();
        const std::size_t order = parameters.quadrature.getVolumeOrder(
          trialFE.getOrder(), polytope.getTransformation().getOrder(),
          Geometry::Polytope::Traits(polytope.getGeometry()).getVertexCount() == dim + 1);
        const auto& qf =
          QF::PolytopeQuadratureFormula::get(order, polytope.getGeometry());
        const auto& quadrature = polytope.getQuadrature(qf);
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
          const HingeState state(
            deformation, innerJacobian.getValue(ip), parameters, m_hingeCoefficient);
          trialJacobian.setIntegrationPoint(ip);
          for (std::size_t local = 0; local < nTrial; ++local)
          {
            const auto gradient = trialJacobian.getBasis(local);
            m_jTrial[local] = state.getJacobianRow(gradient);
            m_qTrial[local] = state.getDistortionRow(gradient);
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
              m_jTest[local] = state.getJacobianRow(gradient);
              m_qTest[local] = state.getDistortionRow(gradient);
            }
          }

          const Real weight = qf.getWeight(q) * point.getDistortion();
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

      /**
       * @brief Returns an entry of the assembled local system.
       * @returns Integral computed by the quadrature rule.
       * @param trial Local trial basis index.
       * @param test Local test basis index.
       */
      ScalarType integrate(std::size_t trial, std::size_t test) final override
      {
        return m_matrix(
          static_cast<Eigen::Index>(test), static_cast<Eigen::Index>(trial));
      }

      /**
       * @brief Returns the integration region.
       * @returns The integration region.
       */
      Geometry::Region getRegion() const final override
      {
        return Geometry::Region::Cells;
      }

      /**
       * @brief Clones this integrator.
       * @returns Newly allocated copy owned by the caller.
       */
      HingeMetric* copy() const noexcept final override
      {
        return new HingeMetric(*this);
      }

    private:
      std::reference_wrapper<const TrialFunction> m_du;
      std::reference_wrapper<const TestFunction> m_z;
      std::reference_wrapper<const Displacement> m_current;
      std::reference_wrapper<const Displacement> m_inner;
      std::reference_wrapper<const Parameters> m_parameters;
      Real m_hingeCoefficient;
      const Geometry::Polytope* m_polytope = nullptr;
      std::vector<Real> m_jTrial;
      std::vector<Real> m_qTrial;
      std::vector<Real> m_jTest;
      std::vector<Real> m_qTest;
      Math::Matrix<Real> m_matrix;
  };
}

#endif
