/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ActiveContraction.h
 * @brief Passive hyperelastic law plus active fiber contraction.
 */
#ifndef RODIN_SOLID_CONSTITUTIVE_ACTIVECONTRACTION_H
#define RODIN_SOLID_CONSTITUTIVE_ACTIVECONTRACTION_H

#include <cmath>

#include "Rodin/Types.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Solid/Local/FiberKinematics.h"
#include "Rodin/Solid/Local/ConstitutivePoint.h"

#include "HyperElasticLaw.h"
#include "ActiveFiberLaw.h"

namespace Rodin::Solid
{
  /**
   * @brief Adds active fiber stress to a passive law.
   *
   * The passive law can be any Rodin hyperelastic law. The active contribution is
   * aligned with `Tags::FiberDirection` and uses `Tags::ActiveExtension`; if
   * dynamic tags are present (`TimeStep`, `PreviousActiveExtension`,
   * `PreviousActiveGamma`, `PreviousActiveBeta`, `ElectricalActivation`), the
   * condensed dynamic tangent is used.
   */
  template <class PassiveLaw, class ActiveLaw = ActiveFiberLaw>
  class ActiveContraction final : public HyperElasticLaw<ActiveContraction<PassiveLaw, ActiveLaw>>
  {
    public:
      struct Cache
      {
        typename PassiveLaw::Cache passive;
        FiberKinematics fiber;
        Real strain = 0.0;
        Real activeExtension = 0.0;
        typename ActiveLaw::Response active;
      };

      ActiveContraction(
          const PassiveLaw& passiveLaw,
          const ActiveLaw& activeLaw = ActiveLaw())
        : m_passiveLaw(passiveLaw),
          m_activeLaw(activeLaw)
      {}

      const PassiveLaw& getPassiveLaw() const
      {
        return m_passiveLaw;
      }

      const ActiveLaw& getActiveLaw() const
      {
        return m_activeLaw;
      }

      void setCache(Cache& cache, const ConstitutivePoint& cp) const
      {
        m_passiveLaw.setCache(cache.passive, cp);

        cache.fiber = FiberKinematics(cp);
        cache.strain = cache.fiber.strain();
        cache.activeExtension = cp.has<Tags::ActiveExtension>()
          ? cp.get<Tags::ActiveExtension>()
          : m_activeLaw.getParameters().initial.extension;

        if (hasDynamicData(cp))
        {
          typename ActiveLaw::State oldState;
          oldState.gamma = cp.get<Tags::PreviousActiveGamma>();
          oldState.beta = cp.get<Tags::PreviousActiveBeta>();
          const auto newState = m_activeLaw.update(
              cp.get<Tags::TimeStep>(),
              oldState,
              cp.get<Tags::PreviousActiveExtension>(),
              cache.activeExtension,
              cp.get<Tags::ElectricalActivation>());
          cache.active = m_activeLaw.evaluateDynamic(
              cp.get<Tags::TimeStep>(),
              oldState,
              newState,
              cache.strain,
              cp.get<Tags::PreviousActiveExtension>(),
              cache.activeExtension,
              cp.get<Tags::ElectricalActivation>());
        }
        else
        {
          cache.active =
            m_activeLaw.evaluateStatic(cache.strain, cache.activeExtension);
        }
      }

      Real getStrainEnergyDensity(const Cache& cache, const ConstitutivePoint& cp) const
      {
        const Real passiveEnergy =
          m_passiveLaw.getStrainEnergyDensity(cache.passive, cp);
        const Real denom = 1.0 + 2.0 * cache.activeExtension;
        const Real activeEnergy =
          0.5 * m_activeLaw.getParameters().stiffness
          * (cache.strain - cache.activeExtension)
          * (cache.strain - cache.activeExtension)
          / (denom * denom);
        return passiveEnergy + activeEnergy;
      }

      void getFirstPiolaKirchhoffStress(
          Math::SpatialMatrix<Real>& P,
          const Cache& cache,
          const ConstitutivePoint& cp) const
      {
        m_passiveLaw.getFirstPiolaKirchhoffStress(P, cache.passive, cp);
        P = P + cache.active.stress
              * cp.getKinematicState().getDeformationGradient()
              * cache.fiber.tensor();
      }

      void getMaterialTangent(
          Math::SpatialMatrix<Real>& dP,
          const Cache& cache,
          const ConstitutivePoint& cp,
          const Math::SpatialMatrix<Real>& dF) const
      {
        m_passiveLaw.getMaterialTangent(dP, cache.passive, cp, dF);

        const Real dStrain = cache.fiber.dStrain(dF);
        dP = dP
           + cache.active.stress * dF * cache.fiber.tensor()
           + cache.active.tangent * dStrain
             * cp.getKinematicState().getDeformationGradient()
             * cache.fiber.tensor();
      }

    private:
      static bool hasDynamicData(const ConstitutivePoint& cp)
      {
        return cp.has<Tags::TimeStep>()
            && cp.has<Tags::PreviousActiveExtension>()
            && cp.has<Tags::PreviousActiveGamma>()
            && cp.has<Tags::PreviousActiveBeta>()
            && cp.has<Tags::ElectricalActivation>();
      }

      PassiveLaw m_passiveLaw;
      ActiveLaw m_activeLaw;
  };
}

#endif
