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
        Math::SpatialVector<Real> fiber;
        Math::SpatialMatrix<Real> fiberTensor;
        Math::SpatialVector<Real> deformedFiber;
        Real strain1D = 0.0;
        Real activeExtension = 0.0;
        Real activeStress = 0.0;
        Real condensedStressTangent = 0.0;
        typename ActiveLaw::Response activeResponse;
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

        const auto& state = cp.getKinematicState();
        const size_t d = state.getDimension();
        const auto d8 = static_cast<std::uint8_t>(d);

        cache.fiber = getFiberDirection(cp, d);
        cache.fiberTensor.resize(d8, d8);
        cache.fiberTensor.setZero();
        for (size_t i = 0; i < d; ++i)
          for (size_t j = 0; j < d; ++j)
            cache.fiberTensor(static_cast<std::uint8_t>(i), static_cast<std::uint8_t>(j)) =
              cache.fiber[i] * cache.fiber[j];

        cache.deformedFiber = state.getDeformationGradient() * cache.fiber;
        cache.strain1D =
          0.5 * (cache.fiber.dot(state.getRightCauchyGreenTensor() * cache.fiber) - 1.0);
        cache.activeExtension = cp.has<Tags::ActiveExtension>()
          ? cp.get<Tags::ActiveExtension>()
          : m_activeLaw.getParameters().initialExtension;

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
          cache.activeResponse = m_activeLaw.evaluateDynamic(
              cp.get<Tags::TimeStep>(),
              oldState,
              newState,
              cache.strain1D,
              cp.get<Tags::PreviousActiveExtension>(),
              cache.activeExtension,
              cp.get<Tags::ElectricalActivation>());
        }
        else
        {
          cache.activeResponse =
            m_activeLaw.evaluateStatic(cache.strain1D, cache.activeExtension);
        }

        cache.activeStress = cache.activeResponse.stress;
        cache.condensedStressTangent = cache.activeResponse.condensedStressTangent;
      }

      Real getStrainEnergyDensity(const Cache& cache, const ConstitutivePoint& cp) const
      {
        const Real passiveEnergy =
          m_passiveLaw.getStrainEnergyDensity(cache.passive, cp);
        const Real denom = 1.0 + 2.0 * cache.activeExtension;
        const Real activeEnergy =
          0.5 * m_activeLaw.getParameters().stiffness
          * (cache.strain1D - cache.activeExtension)
          * (cache.strain1D - cache.activeExtension)
          / (denom * denom);
        return passiveEnergy + activeEnergy;
      }

      void getFirstPiolaKirchhoffStress(
          Math::SpatialMatrix<Real>& P,
          const Cache& cache,
          const ConstitutivePoint& cp) const
      {
        m_passiveLaw.getFirstPiolaKirchhoffStress(P, cache.passive, cp);
        P = P + cache.activeStress
              * cp.getKinematicState().getDeformationGradient()
              * cache.fiberTensor;
      }

      void getMaterialTangent(
          Math::SpatialMatrix<Real>& dP,
          const Cache& cache,
          const ConstitutivePoint& cp,
          const Math::SpatialMatrix<Real>& dF) const
      {
        m_passiveLaw.getMaterialTangent(dP, cache.passive, cp, dF);

        const auto dFa = dF * cache.fiber;
        const Real dStrain1D = cache.deformedFiber.dot(dFa);
        dP = dP
           + cache.activeStress * dF * cache.fiberTensor
           + cache.condensedStressTangent * dStrain1D
             * cp.getKinematicState().getDeformationGradient()
             * cache.fiberTensor;
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

      static Math::SpatialVector<Real> getFiberDirection(
          const ConstitutivePoint& cp,
          size_t d)
      {
        Math::SpatialVector<Real> fiber(static_cast<std::uint8_t>(d));
        if (cp.has<Tags::FiberDirection>())
          fiber = cp.get<Tags::FiberDirection>();
        else
        {
          fiber.setZero();
          fiber[0] = 1.0;
        }

        const Real norm = std::sqrt(fiber.dot(fiber));
        if (norm > 0.0)
          fiber = (1.0 / norm) * fiber;
        return fiber;
      }

      PassiveLaw m_passiveLaw;
      ActiveLaw m_activeLaw;
  };
}

#endif
