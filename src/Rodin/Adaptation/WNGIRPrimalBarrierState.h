/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_WNGIRPRIMALBARRIERSTATE_H
#define RODIN_ADAPTATION_WNGIRPRIMALBARRIERSTATE_H

#include "CellDeformation.h"
#include "WNGIRParameters.h"

namespace Rodin::Adaptation::Detail
{
  /// @brief Pointwise slacks and Newton coefficients of the primal QP barrier.
  class WNGIRPrimalBarrierState
  {
    public:
      WNGIRPrimalBarrierState(const CellDeformation& deformation,
        const Math::SpatialMatrix<Real>& innerGradient,
        const WNGIRParameters& parameters, Real barrierCoefficient)
      {
        if (!deformation.isAdmissible())
          return;
        m_jAction = -deformation.getJacobianAction(innerGradient);
        m_qAction = deformation.getRelativeDistortionAction(innerGradient);
        m_jSlack = deformation.getJacobian() - parameters.jSafe - m_jAction;
        m_qSlack = parameters.qMax -
          deformation.getRelativeDistortion() - m_qAction;
        if (m_jSlack <= Real(0) || m_qSlack <= Real(0))
          return;

        if (parameters.includeAdmissibilityMetric && parameters.gammaJ > Real(0))
        {
          const Real coefficient = barrierCoefficient * parameters.gammaJ;
          m_jHessian = coefficient / (m_jSlack * m_jSlack);
          m_jForce = m_jHessian * m_jAction - coefficient / m_jSlack;
        }
        if (parameters.includeAdmissibilityMetric && parameters.gammaQ > Real(0))
        {
          const Real coefficient = barrierCoefficient * parameters.gammaQ;
          m_qHessian = coefficient / (m_qSlack * m_qSlack);
          m_qForce = m_qHessian * m_qAction - coefficient / m_qSlack;
        }
        m_feasible = true;
      }

      bool isFeasible() const { return m_feasible; }
      Real getJacobianAction() const { return m_jAction; }
      Real getDistortionAction() const { return m_qAction; }
      Real getJacobianSlack() const { return m_jSlack; }
      Real getDistortionSlack() const { return m_qSlack; }
      Real getJacobianHessian() const { return m_jHessian; }
      Real getDistortionHessian() const { return m_qHessian; }
      Real getJacobianForce() const { return m_jForce; }
      Real getDistortionForce() const { return m_qForce; }

    private:
      bool m_feasible = false;
      Real m_jAction = 0;
      Real m_qAction = 0;
      Real m_jSlack = 0;
      Real m_qSlack = 0;
      Real m_jHessian = 0;
      Real m_qHessian = 0;
      Real m_jForce = 0;
      Real m_qForce = 0;
  };
}

#endif
