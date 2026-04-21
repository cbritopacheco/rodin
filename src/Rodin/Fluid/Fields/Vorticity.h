/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Vorticity.h
 * @brief Vorticity tensor/vector helpers.
 */
#ifndef RODIN_FLUID_FIELDS_VORTICITY_H
#define RODIN_FLUID_FIELDS_VORTICITY_H

#include "Rodin/Types.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/SpatialMatrix.h"

#include "Rodin/Fluid/Local/FlowPoint.h"

namespace Rodin::Fluid
{
  /**
   * @brief Computes vorticity quantities from the velocity gradient.
   */
  class Vorticity
  {
    public:
      /**
       * @brief Computes the spin tensor
       * @f$\mathbf{W}=\frac12(\nabla\mathbf{u}-\nabla\mathbf{u}^T)@f$.
       */
      void getVorticityTensor(Math::SpatialMatrix<Real>& W, const FlowPoint& fp) const
      {
        const auto& gradU = fp.getVelocityGradient();
        W = 0.5 * (gradU - gradU.transpose());
      }

      /**
       * @brief Computes the vorticity vector/scalar curl.
       *
       * - In 3D: returns @f$\nabla\times\mathbf{u}\in\mathbb{R}^3@f$.
       * - In 2D: returns one-component vector with out-of-plane vorticity
       *   @f$\omega_z = \partial_x u_y - \partial_y u_x@f$.
       * - In 1D: returns zero.
       */
      void getVorticityVector(Math::SpatialVector<Real>& omega, const FlowPoint& fp) const
      {
        const auto& gradU = fp.getVelocityGradient();
        const size_t d = fp.getDimension();

        if (d == 3)
        {
          omega.resize(3);
          omega(0) = gradU(2, 1) - gradU(1, 2);
          omega(1) = gradU(0, 2) - gradU(2, 0);
          omega(2) = gradU(1, 0) - gradU(0, 1);
        }
        else if (d == 2)
        {
          omega.resize(1);
          omega(0) = gradU(1, 0) - gradU(0, 1);
        }
        else
        {
          omega.resize(1);
          omega(0) = 0.0;
        }
      }
  };
}

#endif
