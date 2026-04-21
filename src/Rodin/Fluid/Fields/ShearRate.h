/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ShearRate.h
 * @brief Scalar shear-rate helper for generalized Newtonian rheology.
 */
#ifndef RODIN_FLUID_FIELDS_SHEARRATE_H
#define RODIN_FLUID_FIELDS_SHEARRATE_H

#include "Rodin/Types.h"
#include "Rodin/Math/Common.h"
#include "Rodin/Math/SpatialMatrix.h"

#include "Rodin/Fluid/Fields/StrainRate.h"
#include "Rodin/Fluid/Local/FlowPoint.h"

namespace Rodin::Fluid
{
  /**
   * @brief Computes the scalar shear-rate invariant.
   *
   * Convention:
   * @f[
   *   \dot{\gamma} = \sqrt{2\,\mathbf{D}:\mathbf{D}}
   * @f]
   * where @f$\mathbf{D}@f$ is the symmetric strain-rate tensor.
   */
  class ShearRate
  {
    public:
      Real getShearRate(const FlowPoint& fp) const
      {
        Math::SpatialMatrix<Real> D;
        m_strainRate.getStrainRate(D, fp);
        return getShearRate(D);
      }

      Real getShearRate(const Math::SpatialMatrix<Real>& D) const
      {
        return Math::sqrt(2.0 * D.dot(D));
      }

    private:
      StrainRate m_strainRate;
  };
}

#endif
