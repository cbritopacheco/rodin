/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file GreenLagrangeStrain.h
 * @brief Green-Lagrange strain tensor computation.
 *
 * Computes the Green-Lagrange strain tensor:
 * @f[
 *   \mathbf{E} = \frac{1}{2}(\mathbf{F}^T \mathbf{F} - \mathbf{I})
 *              = \frac{1}{2}(\mathbf{C} - \mathbf{I})
 * @f]
 */
#ifndef RODIN_SOLID_POSTPROCESSING_GREENLAGRANGESTRAIN_H
#define RODIN_SOLID_POSTPROCESSING_GREENLAGRANGESTRAIN_H

#include "Rodin/Types.h"
#include "Rodin/Math/Matrix.h"

#include "Rodin/Solid/Kinematics/KinematicState.h"

namespace Rodin::Solid
{
  /**
   * @brief Computes the Green-Lagrange strain tensor from a kinematic state.
   *
   * @f[
   *   \mathbf{E} = \frac{1}{2}(\mathbf{C} - \mathbf{I})
   * @f]
   *
   * ## Usage
   * @code
   * KinematicState state(d);
   * state.setDisplacementGradient(gradU).update();
   * GreenLagrangeStrain glStrain;
   * Math::Matrix<Real> E;
   * glStrain.getGreenLagrangeStrain(E, state);
   * @endcode
   */
  class GreenLagrangeStrain
  {
    public:
      GreenLagrangeStrain() = default;

      /**
       * @brief Computes @f$ \mathbf{E} = \tfrac{1}{2}(\mathbf{C} - \mathbf{I}) @f$.
       * @param[out] E Output strain tensor
       * @param[in] state Kinematic state
       */
      void getGreenLagrangeStrain(Math::Matrix<Real>& E, const KinematicState& state) const
      {
        const auto& C = state.getRightCauchyGreenTensor();
        const size_t d = state.getDimension();
        E = 0.5 * (C - Math::Matrix<Real>::Identity(d, d));
      }
  };
}

#endif
