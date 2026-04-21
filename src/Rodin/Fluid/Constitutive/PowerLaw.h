/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file PowerLaw.h
 * @brief Generalized Newtonian power-law rheology.
 */
#ifndef RODIN_FLUID_CONSTITUTIVE_POWERLAW_H
#define RODIN_FLUID_CONSTITUTIVE_POWERLAW_H

#include "Rodin/Types.h"
#include "Rodin/Math/Common.h"
#include "Rodin/Math/SpatialMatrix.h"

#include "Rodin/Fluid/Fields/StrainRate.h"
#include "Rodin/Fluid/Fields/ShearRate.h"
#include "Rodin/Fluid/Constitutive/RheologyLaw.h"

namespace Rodin::Fluid
{
  /**
   * @brief Power-law generalized Newtonian rheology.
   *
   * Apparent viscosity:
   * @f[
   *   \mu_\text{eff}(\dot\gamma)
   *    = K\,(\dot\gamma_r)^{n-1},\qquad
   *   \dot\gamma_r = \sqrt{\dot\gamma^2 + \varepsilon^2}
   * @f]
   * with @f$K@f$ consistency index, @f$n@f$ flow index, and @f$\varepsilon@f$
   * a small regularization to avoid singular behavior at vanishing shear.
   */
  class PowerLaw final : public RheologyLaw<PowerLaw>
  {
    public:
      struct Cache
      {
        Math::SpatialMatrix<Real> D;
        Real shearRate;
        Real regularizedShearRate;
        Real effectiveViscosity;
      };

      PowerLaw(Real consistencyIndex, Real flowIndex, Real regularization = 1e-12)
        : m_k(consistencyIndex),
          m_n(flowIndex),
          m_eps(regularization)
      {}

      PowerLaw(const PowerLaw&) = default;
      PowerLaw(PowerLaw&&) = default;

      Real getConsistencyIndex() const
      {
        return m_k;
      }

      Real getFlowIndex() const
      {
        return m_n;
      }

      Real getRegularization() const
      {
        return m_eps;
      }

      void setCache(Cache& cache, const FlowPoint& fp) const
      {
        m_strainRate.getStrainRate(cache.D, fp);
        cache.shearRate = m_shearRate.getShearRate(cache.D);
        cache.regularizedShearRate = Math::sqrt(
            cache.shearRate * cache.shearRate + m_eps * m_eps);
        cache.effectiveViscosity = m_k * Math::pow(cache.regularizedShearRate, m_n - 1.0);
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

        const Real DdD = cache.D.dot(dD);
        const Real g = cache.regularizedShearRate;

        const Real dmu =
            2.0 * m_k * (m_n - 1.0)
          * Math::pow(g, m_n - 3.0)
          * DdD;

        dtau = 2.0 * cache.effectiveViscosity * dD + 2.0 * dmu * cache.D;
      }

    private:
      Real m_k;
      Real m_n;
      Real m_eps;
      StrainRate m_strainRate;
      ShearRate m_shearRate;
  };
}

#endif
