/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file CarreauYasuda.h
 * @brief Carreau-Yasuda generalized Newtonian rheology.
 */
#ifndef RODIN_FLUID_CONSTITUTIVE_CARREAUYASUDA_H
#define RODIN_FLUID_CONSTITUTIVE_CARREAUYASUDA_H

#include "Rodin/Types.h"
#include "Rodin/Math/Common.h"
#include "Rodin/Math/SpatialMatrix.h"

#include "Rodin/Fluid/Fields/StrainRate.h"
#include "Rodin/Fluid/Fields/ShearRate.h"
#include "Rodin/Fluid/Constitutive/RheologyLaw.h"

namespace Rodin::Fluid
{
  /**
   * @brief Carreau-Yasuda generalized Newtonian rheology.
   *
   * Apparent viscosity:
   * @f[
   *   \mu_\text{eff}(\dot\gamma)
   *   = \mu_\infty + (\mu_0 - \mu_\infty)
   *     \left[1 + (\lambda\dot\gamma_r)^a\right]^{\frac{n-1}{a}},
   *   \qquad \dot\gamma_r = \sqrt{\dot\gamma^2 + \varepsilon^2}
   * @f]
   */
  class CarreauYasuda final : public RheologyLaw<CarreauYasuda>
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

      CarreauYasuda(
          Real zeroShearViscosity,
          Real infiniteShearViscosity,
          Real timeConstant,
          Real flowIndex,
          Real yasudaExponent,
          Real regularization = 1e-12)
        : m_mu0(zeroShearViscosity),
          m_muInf(infiniteShearViscosity),
          m_lambda(timeConstant),
          m_n(flowIndex),
          m_a(yasudaExponent),
          m_eps(regularization)
      {}

      CarreauYasuda(const CarreauYasuda&) = default;
      CarreauYasuda(CarreauYasuda&&) = default;

      void setCache(Cache& cache, const FlowPoint& fp) const
      {
        m_strainRate.getStrainRate(cache.D, fp);
        cache.shearRate = m_shearRate.getShearRate(cache.D);
        cache.regularizedShearRate = Math::sqrt(
            cache.shearRate * cache.shearRate + m_eps * m_eps);

        const Real xi = m_lambda * cache.regularizedShearRate;
        const Real q = 1.0 + Math::pow(xi, m_a);
        const Real exponent = (m_n - 1.0) / m_a;

        cache.effectiveViscosity = m_muInf + (m_mu0 - m_muInf) * Math::pow(q, exponent);

        cache.dEffectiveViscosityDRegularizedShearRate =
            (m_mu0 - m_muInf)
          * (m_n - 1.0)
          * Math::pow(q, exponent - 1.0)
          * Math::pow(xi, m_a - 1.0)
          * m_lambda;
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
      Real m_mu0;
      Real m_muInf;
      Real m_lambda;
      Real m_n;
      Real m_a;
      Real m_eps;
      StrainRate m_strainRate;
      ShearRate m_shearRate;
  };
}

#endif
