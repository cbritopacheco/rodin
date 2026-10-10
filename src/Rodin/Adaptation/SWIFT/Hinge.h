/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_SWIFT_HINGE_H
#define RODIN_ADAPTATION_SWIFT_HINGE_H

#include <limits>

#include "../CellDeformation.h"
#include "Parameters.h"

namespace Rodin::Adaptation::SWIFT
{
  /// @brief Pointwise slacks and Newton coefficients of the affine quadratic hinges.
  class HingeState
  {
    public:
      /**
       * @brief Constructs the SWIFT hinge state.
       * @param deformation Frozen outer deformation.
       * @param innerGradient Gradient of the proposed inner increment.
       * @param parameters Quality bounds, guard widths and relative hinge weights.
       * @param hingeCoefficient Effective penalty coefficient for this outer model.
       */
      HingeState(const CellDeformation& deformation,
        const Math::SpatialMatrix<Real>& innerGradient, const Parameters& parameters,
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

      /**
       * @brief Affine quality energy with the construction parameters, evaluated only on demand.
       * @param parameters Quality guards and relative hinge weights used at construction.
       * @param hingeCoefficient Effective penalty coefficient used at construction.
       * @param weightJ Frozen Jacobian witness measure.
       * @param weightQ Frozen distortion witness measure.
       * @returns The squared-hinge energy, or infinity for an inadmissible frozen state.
       */
      Real getEnergy(const Parameters& parameters, Real hingeCoefficient,
        Real weightJ = Real(1), Real weightQ = Real(1)) const
      {
        if (!isAdmissible())
          return std::numeric_limits<Real>::infinity();
        const auto energy = [&](Real slack, Real delta, Real weight) {
          const Real violation = std::max(Real(0), Real(1) - slack / delta);
          return Real(0.5) * hingeCoefficient * weight * violation * violation;
        };
        return weightJ * energy(m_jSlack,
                 parameters.model.qualityGuard * (Real(1) - parameters.model.jacobian),
                 parameters.model.jacobianWeight) +
          weightQ * energy(m_qSlack,
            parameters.model.qualityGuard * (parameters.model.distortion - Real(1)),
            parameters.model.distortionWeight);
      }

      /**
       * @brief Whether the frozen deformation permits evaluation; affine slacks may be negative.
       * @returns Whether the frozen deformation has a positive Jacobian.
       */
      bool isAdmissible() const
      {
        return m_rowDeformation.isAdmissible();
      }
      /**
       * @brief Negative Jacobian differential at the frozen outer state.
       * @param gradient Incremental deformation gradient.
       * @returns The negative determinant differential.
       */
      Real getJacobianRow(const Math::SpatialMatrix<Real>& gradient) const
      {
        return -m_rowDeformation.getJacobianAction(gradient);
      }
      /**
       * @brief Distortion differential at the frozen outer state.
       * @param gradient Incremental deformation gradient.
       * @returns The relative-distortion differential.
       */
      Real getDistortionRow(const Math::SpatialMatrix<Real>& gradient) const
      {
        return m_rowDeformation.getRelativeDistortionAction(gradient);
      }
      /**
       * @brief The jacobian action.
       * @returns The negative determinant differential on the inner increment.
       */
      Real getJacobianAction() const
      {
        return m_jAction;
      }
      /**
       * @brief The distortion action.
       * @returns The distortion differential on the inner increment.
       */
      Real getDistortionAction() const
      {
        return m_qAction;
      }
      /**
       * @brief The jacobian slack.
       * @returns The affine distance above the Jacobian floor.
       */
      Real getJacobianSlack() const
      {
        return m_jSlack;
      }
      /**
       * @brief The distortion slack.
       * @returns The affine distance below the distortion ceiling.
       */
      Real getDistortionSlack() const
      {
        return m_qSlack;
      }
      /**
       * @brief The jacobian hessian.
       * @returns The active Jacobian hinge curvature coefficient, or zero.
       */
      Real getJacobianHessian() const
      {
        return m_jHessian;
      }
      /**
       * @brief The distortion hessian.
       * @returns The active distortion hinge curvature coefficient, or zero.
       */
      Real getDistortionHessian() const
      {
        return m_qHessian;
      }
      /**
       * @brief The jacobian force.
       * @returns The coefficient of the Jacobian row in the Newton load.
       */
      Real getJacobianForce() const
      {
        return m_jForce;
      }
      /**
       * @brief The distortion force.
       * @returns The coefficient of the distortion row in the Newton load.
       */
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
