/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_WNGIR_HINGE_H
#define RODIN_ADAPTATION_WNGIR_HINGE_H

#include <limits>

#include "../CellDeformation.h"
#include "Parameters.h"

namespace Rodin::Adaptation
{
  /// @brief Pointwise slacks and Newton coefficients of the affine quadratic hinges.
  class WNGIRHingeState
  {
    public:
      /// @brief Constructs the WNGIR hinge state.
      WNGIRHingeState(const CellDeformation& deformation,
        const Math::SpatialMatrix<Real>& innerGradient, const WNGIRParameters& parameters,
        Real hingeCoefficient)
        : m_rowDeformation(deformation)
      {
        if (!deformation.isAdmissible())
          return;
        m_jAction = getJacobianRow(innerGradient);
        m_qAction = getDistortionRow(innerGradient);
        m_jSlack = deformation.getJacobian() - parameters.model.jacobian - m_jAction;
        m_qSlack =
          parameters.model.distortion - deformation.getRelativeDistortion() - m_qAction;
        const auto coefficients = [&](Real action, Real slack, Real delta, Real weight,
                                    Real& hessian, Real& force) {
          if (weight <= Real(0) || slack >= delta)
            return;
          hessian = hingeCoefficient * weight / (delta * delta);
          force = hessian * (action - (delta - slack));
        };
        coefficients(m_jAction, m_jSlack,
          parameters.model.qualityGuard * (Real(1) - parameters.model.jacobian),
          parameters.model.jacobianWeight, m_jHessian, m_jForce);
        coefficients(m_qAction, m_qSlack,
          parameters.model.qualityGuard * (parameters.model.distortion - Real(1)),
          parameters.model.distortionWeight, m_qHessian, m_qForce);
      }

      /// @brief Affine quality energy with the construction parameters, evaluated only on demand.
      Real getEnergy(const WNGIRParameters& parameters, Real hingeCoefficient) const
      {
        if (!isAdmissible())
          return std::numeric_limits<Real>::infinity();
        const auto energy = [&](Real slack, Real delta, Real weight) {
          const Real violation = std::max(Real(0), Real(1) - slack / delta);
          return Real(0.5) * hingeCoefficient * weight * violation * violation;
        };
        return energy(m_jSlack,
                 parameters.model.qualityGuard * (Real(1) - parameters.model.jacobian),
                 parameters.model.jacobianWeight) +
          energy(m_qSlack,
            parameters.model.qualityGuard * (parameters.model.distortion - Real(1)),
            parameters.model.distortionWeight);
      }

      /// @brief Whether the frozen deformation permits evaluation; affine slacks may be negative.
      bool isAdmissible() const
      {
        return m_rowDeformation.isAdmissible();
      }
      /// Negative Jacobian differential at the frozen outer state.
      Real getJacobianRow(const Math::SpatialMatrix<Real>& gradient) const
      {
        return -m_rowDeformation.getJacobianAction(gradient);
      }
      /// Distortion differential at the frozen outer state.
      Real getDistortionRow(const Math::SpatialMatrix<Real>& gradient) const
      {
        return m_rowDeformation.getRelativeDistortionAction(gradient);
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
      CellDeformation m_rowDeformation;
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
