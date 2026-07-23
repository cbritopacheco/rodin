/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_WNGIRRESIDUALSTATE_H
#define RODIN_ADAPTATION_WNGIRRESIDUALSTATE_H

#include <cmath>

#include "Rodin/Math.h"
#include "Rodin/Types.h"
#include "Rodin/Variational/IntegrationPoint.h"

namespace Rodin::Adaptation::Detail
{
  /// @brief Pointwise Welsch residual state at a deformed interface point.
  class WNGIRResidualState
  {
    public:
      template <class PhiType, class GradType, class DeformationType>
      /// @brief Constructs the w n g i r residual state.
      WNGIRResidualState(const PhiType& phi, const GradType& grad,
        const DeformationType& deformation, const Variational::IntegrationPoint& ip,
        Real sigma2, bool weighted)
      {
        const auto& moved = deformation.getMovedPoint(ip);
        m_residual = phi.getValue(moved);
        m_gradient = grad.getValue(moved);
        m_weight = weighted ? std::exp(-m_residual * m_residual / sigma2) : Real(1);
      }

      /// @brief The residual.
      Real getResidual() const
      {
        return m_residual;
      }

      /// @brief The gradient.
      const Math::SpatialVector<Real>& getGradient() const
      {
        return m_gradient;
      }

      /// @brief The weight.
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
