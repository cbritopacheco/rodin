/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Bingham.h
 * @brief Regularized Bingham rheology (Papanastasiou regularization).
 */
#ifndef RODIN_FLUID_CONSTITUTIVE_BINGHAM_H
#define RODIN_FLUID_CONSTITUTIVE_BINGHAM_H

#include "Rodin/Types.h"
#include "Rodin/Math/Common.h"
#include "Rodin/Math/SpatialMatrix.h"

#include "Rodin/Fluid/Fields/StrainRate.h"
#include "Rodin/Fluid/Fields/ShearRate.h"
#include "Rodin/Fluid/Constitutive/RheologyLaw.h"

namespace Rodin::Fluid
{
  /**
   * @brief Bingham rheology using Papanastasiou regularization.
   *
   * Apparent viscosity:
   * @f[
   *   \mu_\text{eff}(\dot\gamma)
   *   = \mu_p + \tau_y\,\frac{1 - e^{-m\dot\gamma_r}}{\dot\gamma_r},
   *   \qquad \dot\gamma_r = \sqrt{\dot\gamma^2 + \varepsilon^2}
   * @f]
   * where @f$\mu_p@f$ is plastic viscosity, @f$\tau_y@f$ is yield stress and
   * @f$m@f$ is the Papanastasiou regularization parameter.
   */
  class Bingham final : public RheologyLaw<Bingham>
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

      Bingham(
          Real plasticViscosity,
          Real yieldStress,
          Real regularization,
          Real shearRateFloor = 1e-12)
        : m_muPlastic(plasticViscosity),
          m_tauY(yieldStress),
          m_m(regularization),
          m_eps(shearRateFloor)
      {}

      Bingham(const Bingham&) = default;
      Bingham(Bingham&&) = default;

      void setCache(Cache& cache, const FlowPoint& fp) const
      {
        m_strainRate.getStrainRate(cache.D, fp);
        cache.shearRate = m_shearRate.getShearRate(cache.D);
        cache.regularizedShearRate = Math::sqrt(
            cache.shearRate * cache.shearRate + m_eps * m_eps);

        const Real g = cache.regularizedShearRate;
        const Real expTerm = Math::exp(-m_m * g);
        const Real f = (1.0 - expTerm) / g;

        cache.effectiveViscosity = m_muPlastic + m_tauY * f;
        cache.dEffectiveViscosityDRegularizedShearRate =
            m_tauY * ((m_m * expTerm * g - (1.0 - expTerm)) / (g * g));
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
      Real m_muPlastic;
      Real m_tauY;
      Real m_m;
      Real m_eps;
      StrainRate m_strainRate;
      ShearRate m_shearRate;
  };
}

#endif
