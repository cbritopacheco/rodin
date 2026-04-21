/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file StrainRate.h
 * @brief Symmetric strain-rate tensor helper.
 */
#ifndef RODIN_FLUID_FIELDS_STRAINRATE_H
#define RODIN_FLUID_FIELDS_STRAINRATE_H

#include "Rodin/Types.h"
#include "Rodin/Math/SpatialMatrix.h"

#include "Rodin/Fluid/Local/FlowPoint.h"

namespace Rodin::Fluid
{
  /**
   * @brief Computes the symmetric strain-rate tensor
   * @f$\mathbf{D}(\mathbf{u}) = \frac12(\nabla\mathbf{u} + \nabla\mathbf{u}^T)@f$.
   */
  class StrainRate
  {
    public:
      void getStrainRate(Math::SpatialMatrix<Real>& D, const FlowPoint& fp) const
      {
        const auto& gradU = fp.getVelocityGradient();
        D = 0.5 * (gradU + gradU.transpose());
      }

      void getStrainRate(
          Math::SpatialMatrix<Real>& D,
          const Math::SpatialMatrix<Real>& velocityGradient) const
      {
        D = 0.5 * (velocityGradient + velocityGradient.transpose());
      }
  };
}

#endif
