/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_WNGIR_RESIDUAL_H
#define RODIN_ADAPTATION_WNGIR_RESIDUAL_H

#include "Rodin/Math.h"
#include "Rodin/Types.h"
#include "Rodin/Variational/IntegrationPoint.h"

#include "Loss.h"

namespace Rodin::Adaptation::WNGIR
{
  /// @brief Pointwise robust residual state at a deformed interface point.
  class ResidualState
  {
    public:
      template <class PhiType, class GradType, class DeformationType>
      /**
       * @brief Constructs the WNGIR residual state.
       * @param deformation Deformation data.
       * @param ip Integration point at which the expression is evaluated.
       * @param phi Observation field.
       * @param grad Gradient of the observation field.
       * @param loss Loss applied to the observation residual.
       */
      ResidualState(const PhiType& phi, const GradType& grad,
        const DeformationType& deformation, const Variational::IntegrationPoint& ip,
        const Loss& loss)
      {
        const auto& moved = deformation.getMovedPoint(ip);
        m_residual = phi.getValue(moved);
        m_gradient = grad.getValue(moved);
        m_weight = loss.getWeight(m_residual);
      }

      /**
       * @brief The residual.
       * @returns The residual.
       */
      Real getResidual() const
      {
        return m_residual;
      }

      /**
       * @brief The gradient.
       * @returns Derivative evaluated at the supplied point.
       */
      const Math::SpatialVector<Real>& getGradient() const
      {
        return m_gradient;
      }

      /**
       * @brief The weight.
       * @returns The weight.
       */
      Real getWeight() const
      {
        return m_weight;
      }

    private:
      Real m_residual;
      Math::SpatialVector<Real> m_gradient;
      Real m_weight;
  };
}

#endif
