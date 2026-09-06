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
      /// @brief Constructs the w n g i r primal barrier state.
      WNGIRPrimalBarrierState(const CellDeformation& deformation,
        const Math::SpatialMatrix<Real>& innerGradient, const WNGIRParameters& parameters,
        Real barrierCoefficient)
      {
        if (!deformation.isAdmissible())
          return;
        m_jAction = -deformation.getJacobianAction(innerGradient);
        m_qAction = deformation.getRelativeDistortionAction(innerGradient);
        m_jSlack = deformation.getJacobian() - parameters.jSafe - m_jAction;
        m_qSlack = parameters.qMax - deformation.getRelativeDistortion() - m_qAction;
        if (m_jSlack <= Real(0) || m_qSlack <= Real(0))
          return;

        if (parameters.kappaJ > Real(0))
        {
          const Real coefficient = barrierCoefficient * parameters.kappaJ;
          m_jHessian = coefficient / (m_jSlack * m_jSlack);
          m_jForce = m_jHessian * m_jAction - coefficient / m_jSlack;
        }
        if (parameters.kappaQ > Real(0))
        {
          const Real coefficient = barrierCoefficient * parameters.kappaQ;
          m_qHessian = coefficient / (m_qSlack * m_qSlack);
          m_qForce = m_qHessian * m_qAction - coefficient / m_qSlack;
        }
        m_feasible = true;
      }

      /// @brief Whether feasible.
      bool isFeasible() const
      {
        return m_feasible;
      }
      /// @brief The jacobian action.
      Real getJacobianAction() const
      {
        return m_jAction;
      }
      /// @brief The distortion action.
      Real getDistortionAction() const
      {
        return m_qAction;
      }
      /// @brief The jacobian slack.
      Real getJacobianSlack() const
      {
        return m_jSlack;
      }
      /// @brief The distortion slack.
      Real getDistortionSlack() const
      {
        return m_qSlack;
      }
      /// @brief The jacobian hessian.
      Real getJacobianHessian() const
      {
        return m_jHessian;
      }
      /// @brief The distortion hessian.
      Real getDistortionHessian() const
      {
        return m_qHessian;
      }
      /// @brief The jacobian force.
      Real getJacobianForce() const
      {
        return m_jForce;
      }
      /// @brief The distortion force.
      Real getDistortionForce() const
      {
        return m_qForce;
      }

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
