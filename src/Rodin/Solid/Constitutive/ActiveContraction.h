/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ActiveContraction.h
 * @brief Active fiber-contraction decorator for hyperelastic constitutive laws.
 *
 * This law composes a passive hyperelastic law with an additive active
 * contribution aligned with a user-provided fiber direction:
 *
 * @f[
 *   W_{\mathrm{active}} = \frac{1}{2} T_a \alpha (\lambda_f - \lambda_0)^2,
 *   \quad \lambda_f = \sqrt{I_4},
 *   \quad I_4 = \mathbf{a}_0 \cdot \mathbf{C} \mathbf{a}_0
 * @f]
 *
 * where @f$T_a@f$ is the active tension scale, @f$\alpha@f$ is the activation
 * (injected through @c Tags::Activation), and @f$\mathbf{a}_0@f$ is the fiber
 * direction in the reference configuration (injected through
 * @c Tags::FiberDirection).
 *
 * The resulting first Piola-Kirchhoff stress is
 * @f[
 *   \mathbf{P} = \mathbf{P}_{\mathrm{passive}}
 *   + T_a \alpha\left(1 - \frac{\lambda_0}{\lambda_f}\right)
 *     (\mathbf{F}\mathbf{a}_0) \otimes \mathbf{a}_0.
 * @f]
 *
 * No Schur/local elimination machinery is required: this follows Rodin's
 * existing constitutive/integrator architecture and stays fully local at each
 * quadrature point.
 */
#ifndef RODIN_SOLID_CONSTITUTIVE_ACTIVECONTRACTION_H
#define RODIN_SOLID_CONSTITUTIVE_ACTIVECONTRACTION_H

#include <cmath>
#include <algorithm>
#include <cassert>

#include "Rodin/Types.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialVector.h"

#include "Rodin/Solid/Local/ConstitutivePoint.h"

#include "HyperElasticLaw.h"

namespace Rodin::Solid
{
  /**
   * @brief Decorates a passive hyperelastic law with fiber-aligned active contraction.
   *
   * @tparam PassiveLaw The passive law type (must satisfy HyperElasticLaw API)
   */
  template <class PassiveLaw>
  class ActiveContraction final : public HyperElasticLaw<ActiveContraction<PassiveLaw>>
  {
    public:
      struct Cache
      {
        typename PassiveLaw::Cache passive;
        Math::SpatialVector<Real> a0;
        Math::SpatialVector<Real> f;
        Real activation;
        Real lambdaF;
      };

      /**
       * @param passive Passive constitutive law.
       * @param activeTensionScale Active tension scale @f$T_a@f$.
       * @param referenceFiberStretch Active reference fiber stretch @f$\lambda_0@f$.
       */
      ActiveContraction(
          const PassiveLaw& passive,
          Real activeTensionScale,
          Real referenceFiberStretch = 1.0)
        : m_passive(passive),
          m_activeTensionScale(activeTensionScale),
          m_referenceFiberStretch(referenceFiberStretch)
      {}

      void setCache(Cache& cache, const ConstitutivePoint& cp) const
      {
        m_passive.setCache(cache.passive, cp);

        assert(cp.has<Tags::FiberDirection>());
        cache.a0 = cp.get<Tags::FiberDirection>();
        const Real norm = cache.a0.norm();
        assert(norm > 0.0);
        cache.a0 /= norm;

        cache.activation = cp.has<Tags::Activation>()
          ? cp.get<Tags::Activation>()
          : 1.0;

        const auto& F = cp.getKinematicState().getDeformationGradient();
        cache.f = F * cache.a0;

        const Real I4 = std::max<Real>(cache.f.dot(cache.f), 1e-16);
        cache.lambdaF = std::sqrt(I4);
      }

      Real getStrainEnergyDensity(const Cache& cache, const ConstitutivePoint& cp) const
      {
        const Real passiveEnergy = m_passive.getStrainEnergyDensity(cache.passive, cp);
        const Real delta = cache.lambdaF - m_referenceFiberStretch;
        return passiveEnergy + 0.5 * m_activeTensionScale * cache.activation * delta * delta;
      }

      void getFirstPiolaKirchhoffStress(
          Math::SpatialMatrix<Real>& P,
          const Cache& cache,
          const ConstitutivePoint& cp) const
      {
        m_passive.getFirstPiolaKirchhoffStress(P, cache.passive, cp);

        const Real coef = m_activeTensionScale * cache.activation
                        * (1.0 - m_referenceFiberStretch / cache.lambdaF);

        P += coef * (cache.f * cache.a0.transpose());
      }

      void getMaterialTangent(
          Math::SpatialMatrix<Real>& dP,
          const Cache& cache,
          const ConstitutivePoint& cp,
          const Math::SpatialMatrix<Real>& dF) const
      {
        m_passive.getMaterialTangent(dP, cache.passive, cp, dF);

        const Math::SpatialVector<Real> df = dF * cache.a0;
        const Real fdotdf = cache.f.dot(df);

        const Real coef = m_activeTensionScale * cache.activation
                        * (1.0 - m_referenceFiberStretch / cache.lambdaF);

        const Real dCoef = m_activeTensionScale * cache.activation
                         * m_referenceFiberStretch * fdotdf
                         / (cache.lambdaF * cache.lambdaF * cache.lambdaF);

        dP += dCoef * (cache.f * cache.a0.transpose())
            + coef * (df * cache.a0.transpose());
      }

      const PassiveLaw& getPassiveLaw() const { return m_passive; }
      Real getActiveTensionScale() const { return m_activeTensionScale; }
      Real getReferenceFiberStretch() const { return m_referenceFiberStretch; }

    private:
      PassiveLaw m_passive;
      Real m_activeTensionScale;
      Real m_referenceFiberStretch;
  };
}

#endif
