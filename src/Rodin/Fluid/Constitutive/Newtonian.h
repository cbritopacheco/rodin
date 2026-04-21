/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Newtonian.h
 * @brief Constant-viscosity Newtonian rheology.
 */
#ifndef RODIN_FLUID_CONSTITUTIVE_NEWTONIAN_H
#define RODIN_FLUID_CONSTITUTIVE_NEWTONIAN_H

#include "Rodin/Types.h"
#include "Rodin/Math/SpatialMatrix.h"

#include "Rodin/Fluid/Fields/StrainRate.h"
#include "Rodin/Fluid/Constitutive/RheologyLaw.h"

namespace Rodin::Fluid
{
  /**
   * @brief Newtonian rheology with constant dynamic viscosity @f$\mu@f$.
   *
   * @f[
   *   \boldsymbol{\tau} = 2\mu\mathbf{D}
   * @f]
   */
  class Newtonian final : public RheologyLaw<Newtonian>
  {
    public:
      struct Cache
      {
        Math::SpatialMatrix<Real> D;
      };

      explicit Newtonian(Real dynamicViscosity)
        : m_mu(dynamicViscosity)
      {}

      Newtonian(const Newtonian&) = default;
      Newtonian(Newtonian&&) = default;

      Real getDynamicViscosity() const
      {
        return m_mu;
      }

      void setCache(Cache& cache, const FlowPoint& fp) const
      {
        m_strainRate.getStrainRate(cache.D, fp);
      }

      void getDeviatoricStress(
          Math::SpatialMatrix<Real>& tau,
          const Cache& cache,
          const FlowPoint&) const
      {
        tau = 2.0 * m_mu * cache.D;
      }

      void getTangent(
          Math::SpatialMatrix<Real>& dtau,
          const Cache&,
          const FlowPoint&,
          const Math::SpatialMatrix<Real>& dGradU) const
      {
        const Math::SpatialMatrix<Real> dD = 0.5 * (dGradU + dGradU.transpose());
        dtau = 2.0 * m_mu * dD;
      }

    private:
      Real m_mu;
      StrainRate m_strainRate;
  };
}

#endif
