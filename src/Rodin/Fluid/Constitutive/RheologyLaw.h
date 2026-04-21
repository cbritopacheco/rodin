/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file RheologyLaw.h
 * @brief CRTP base class for fluid rheology constitutive closures.
 */
#ifndef RODIN_FLUID_CONSTITUTIVE_RHEOLOGYLAW_H
#define RODIN_FLUID_CONSTITUTIVE_RHEOLOGYLAW_H

#include "Rodin/Types.h"
#include "Rodin/Math/SpatialMatrix.h"

#include "Rodin/Fluid/Local/FlowPoint.h"

namespace Rodin::Fluid
{
  /**
   * @brief CRTP base class for local fluid rheology laws.
   *
   * Derived laws provide:
   * - @c setCache(Cache&, const FlowPoint&)
   * - @c getDeviatoricStress(SpatialMatrix&, const Cache&, const FlowPoint&)
   * - @c getTangent(SpatialMatrix&, const Cache&, const FlowPoint&, const SpatialMatrix&)
   *
   * Tangent is the directional derivative of deviatoric stress with respect to
   * the velocity gradient perturbation @f$\delta\nabla\mathbf{u}@f$.
   */
  template <class Derived>
  class RheologyLaw
  {
    public:
      template <class Cache>
      void setCache(Cache& cache, const FlowPoint& fp) const
      {
        static_cast<const Derived&>(*this).setCache(cache, fp);
      }

      template <class Cache>
      void getDeviatoricStress(
          Math::SpatialMatrix<Real>& tau,
          const Cache& cache,
          const FlowPoint& fp) const
      {
        static_cast<const Derived&>(*this).getDeviatoricStress(tau, cache, fp);
      }

      template <class Cache>
      void getTangent(
          Math::SpatialMatrix<Real>& dtau,
          const Cache& cache,
          const FlowPoint& fp,
          const Math::SpatialMatrix<Real>& dGradU) const
      {
        static_cast<const Derived&>(*this).getTangent(dtau, cache, fp, dGradU);
      }

    protected:
      RheologyLaw() = default;
      RheologyLaw(const RheologyLaw&) = default;
      RheologyLaw(RheologyLaw&&) = default;
      RheologyLaw& operator=(const RheologyLaw&) = default;
      RheologyLaw& operator=(RheologyLaw&&) = default;
  };
}

#endif
