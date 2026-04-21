/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file HerschelBulkley.h
 * @brief Regularized Herschel-Bulkley rheology.
 */
#ifndef RODIN_FLUID_CONSTITUTIVE_HERSCHELBULKLEY_H
#define RODIN_FLUID_CONSTITUTIVE_HERSCHELBULKLEY_H

#include "Rodin/Types.h"
#include "Rodin/Math/Common.h"
#include "Rodin/Math/SpatialMatrix.h"

#include "Rodin/Fluid/Fields/StrainRate.h"
#include "Rodin/Fluid/Fields/ShearRate.h"
#include "Rodin/Fluid/Constitutive/RheologyLaw.h"

namespace Rodin::Fluid
{
  /**
   * @brief Herschel-Bulkley rheology with Papanastasiou-regularized yield term.
   *
   * Apparent viscosity:
   * @f[
   *   \mu_\text{eff}(\dot\gamma)
   *   = K\,(\dot\gamma_r)^{n-1}
   *   + \tau_y\,\frac{1 - e^{-m\dot\gamma_r}}{\dot\gamma_r},
   *   \qquad \dot\gamma_r = \sqrt{\dot\gamma^2 + \varepsilon^2}
   * @f]
   */
  class HerschelBulkley final : public RheologyLaw<HerschelBulkley>
  {
    public:
      struct Cache
      {
        Math::SpatialMatrix<Real> D;
        Real shearRate;
        Real regularizedShearRate;
        Real effectiveViscosity;
        Real dEffectiveViscosityDRegularizedShearRate;
      };

      HerschelBulkley(
          Real consistencyIndex,
          Real flowIndex,
          Real yieldStress,
          Real regularization,
          Real shearRateFloor = 1e-12)
        : m_k(consistencyIndex),
          m_n(flowIndex),
          m_tauY(yieldStress),
          m_m(regularization),
          m_eps(shearRateFloor)
      {}

      HerschelBulkley(const HerschelBulkley&) = default;
      HerschelBulkley(HerschelBulkley&&) = default;

      void setCache(Cache& cache, const FlowPoint& fp) const
      {
        m_strainRate.getStrainRate(cache.D, fp);
        cache.shearRate = m_shearRate.getShearRate(cache.D);
        cache.regularizedShearRate = Math::sqrt(
            cache.shearRate * cache.shearRate + m_eps * m_eps);

        const Real g = cache.regularizedShearRate;
        const Real expTerm = Math::exp(-m_m * g);
        const Real yieldTerm = (1.0 - expTerm) / g;

        cache.effectiveViscosity =
            m_k * Math::pow(g, m_n - 1.0)
          + m_tauY * yieldTerm;

        cache.dEffectiveViscosityDRegularizedShearRate =
            m_k * (m_n - 1.0) * Math::pow(g, m_n - 2.0)
          + m_tauY * ((m_m * expTerm * g - (1.0 - expTerm)) / (g * g));
      }

      void getDeviatoricStress(
          Math::SpatialMatrix<Real>& tau,
          const Cache& cache,
          const FlowPoint&) const
      {
        tau = 2.0 * cache.effectiveViscosity * cache.D;
      }

      void getTangent(
          Math::SpatialMatrix<Real>& dtau,
          const Cache& cache,
          const FlowPoint&,
          const Math::SpatialMatrix<Real>& dGradU) const
      {
        const Math::SpatialMatrix<Real> dD = 0.5 * (dGradU + dGradU.transpose());
        const Real dGammaRegularized =
            2.0 * cache.D.dot(dD) / cache.regularizedShearRate;
        const Real dMu = cache.dEffectiveViscosityDRegularizedShearRate * dGammaRegularized;

        dtau = 2.0 * cache.effectiveViscosity * dD + 2.0 * dMu * cache.D;
      }

    private:
      Real m_k;
      Real m_n;
      Real m_tauY;
      Real m_m;
      Real m_eps;
      StrainRate m_strainRate;
      ShearRate m_shearRate;
  };
}

#endif
