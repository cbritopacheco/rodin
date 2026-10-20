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
   * active extension, gamma and beta must also be supplied as independent
   * discontinuous fields. The cache exposes their residual and tangent blocks
   * for a coupled global Newton solve; no local solve or condensation occurs.
   * @tparam PassiveLaw Passive constitutive law providing its response cache.
   * @tparam ActiveLaw Three-field active law providing residual and partial tangent blocks.
   */
  template <class PassiveLaw, class ActiveLaw = ActiveFiberLaw>
  class ActiveContraction final
    : public HyperElasticLaw<ActiveContraction<PassiveLaw, ActiveLaw>>
  {
    public:
      /// @brief Cached passive and active quantities at a quadrature point.
      struct Cache
      {
        /// @brief Passive law cache.
        typename PassiveLaw::Cache passive;

        /// @brief Fiber kinematics built from the constitutive point.
        FiberKinematics fiber;

        /// @brief Fiber strain driving the active branch.
        Real strain = 0.0;

        /// @brief Current active extension.
        Real activeExtension = 0.0;

        /// @brief Active branch response.
        typename ActiveLaw::Response active;

        /// @brief Current independent gamma and beta fields.
        typename ActiveLaw::State newState;

        /// @brief Whether the dynamic active equations were evaluated.
        bool dynamic = false;

        /// @brief Uneliminated active-field residual and tangent blocks.
        typename ActiveLaw::Residual residual;
      };

      /**
       * @brief Constructs the coupled active contraction law.
       * @param passiveLaw Passive constitutive law.
       * @param activeLaw Active constitutive law.
       */
      ActiveContraction(
        const PassiveLaw& passiveLaw, const ActiveLaw& activeLaw = ActiveLaw())
        : m_passiveLaw(passiveLaw),
          m_activeLaw(activeLaw)
      {}

      /**
       * @brief Returns the passive hyperelastic law.
       * @returns The passive hyperelastic law.
       */
      const PassiveLaw& getPassiveLaw() const
      {
        return m_passiveLaw;
      }

      /**
       * @brief Returns the active fiber law.
       * @returns The active fiber law.
       */
      const ActiveLaw& getActiveLaw() const
      {
        return m_activeLaw;
      }

      /**
       * @brief Populates cached passive and active quantities.
       * @param cache Storage for the constitutive invariant cache.
       * @param cp Constitutive evaluation point.
       */
      void getCache(Cache& cache, const ConstitutivePoint& cp) const
      {
        m_passiveLaw.getCache(cache.passive, cp);

        cache.fiber = FiberKinematics(cp);
        cache.strain = cache.fiber.strain();

        if (hasDynamicData(cp))
        {
          if (!cp.has<Tags::TimeStep>() || !cp.has<Tags::PreviousActiveExtension>() ||
            !cp.has<Tags::PreviousActiveGamma>() || !cp.has<Tags::PreviousActiveBeta>() ||
            !cp.has<Tags::ElectricalActivation>() || !cp.has<Tags::ActiveExtension>() ||
            !cp.has<Tags::ActiveGamma>() || !cp.has<Tags::ActiveBeta>())
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Dynamic contraction requires independent active fields."
              << Alert::Raise;
          }
          typename ActiveLaw::State oldState;
          oldState.gamma = cp.get<Tags::PreviousActiveGamma>();
          oldState.beta = cp.get<Tags::PreviousActiveBeta>();
          cache.activeExtension = cp.get<Tags::ActiveExtension>();
          cache.newState.gamma = cp.get<Tags::ActiveGamma>();
          cache.newState.beta = cp.get<Tags::ActiveBeta>();
          const bool midpoint = cp.has<Tags::PreviousFiberStrain>();
          const Real strainFactor = midpoint ? Real(0.5) : Real(1);
          const Real strain = midpoint
            ? Real(0.5) * (cache.strain + cp.get<Tags::PreviousFiberStrain>())
            : cache.strain;
          const Real dt = cp.get<Tags::TimeStep>();
          const Real previousExtension = cp.get<Tags::PreviousActiveExtension>();
          const Real activation = cp.get<Tags::ElectricalActivation>();
          m_activeLaw.getResidual(cache.residual, dt, oldState, cache.newState, strain,
            previousExtension, cache.activeExtension, activation, strainFactor);
          cache.active = m_activeLaw.getResponse(cache.residual, strain,
            previousExtension, cache.activeExtension, strainFactor);
          cache.dynamic = true;
        }
        else
        {
          cache.activeExtension = cp.has<Tags::ActiveExtension>()
            ? cp.get<Tags::ActiveExtension>()
            : m_activeLaw.getParameters().initial.extension;
          cache.active = m_activeLaw.evaluateStatic(cache.strain, cache.activeExtension);
          cache.newState = m_activeLaw.initialState();
          cache.dynamic = false;
          cache.residual = typename ActiveLaw::Residual();
        }
      }

      /**
       * @brief Returns the sum of passive and active strain-energy densities.
       * @param cache Constitutive invariant cache populated by getCache().
       * @param cp Constitutive evaluation point.
       * @returns The sum of passive and active strain-energy densities.
       */
      Real getStrainEnergyDensity(const Cache& cache, const ConstitutivePoint& cp) const
      {
        const Real passiveEnergy = m_passiveLaw.getStrainEnergyDensity(cache.passive, cp);
        const Real denom = 1.0 + 2.0 * cache.activeExtension;
        const Real activeEnergy = 0.5 * m_activeLaw.getParameters().stiffness *
          (cache.strain - cache.activeExtension) *
          (cache.strain - cache.activeExtension) / (denom * denom);
        return passiveEnergy + activeEnergy;
      }

      /**
       * @brief Adds the active contribution to the first Piola-Kirchhoff stress.
       * @param P Storage for the first Piola-Kirchhoff stress.
       * @param cache Constitutive invariant cache populated by getCache().
       * @param cp Constitutive evaluation point.
       */
      void getFirstPiolaKirchhoffStress(Math::SpatialMatrix<Real>& P, const Cache& cache,
        const ConstitutivePoint& cp) const
      {
        m_passiveLaw.getFirstPiolaKirchhoffStress(P, cache.passive, cp);
        P = P +
          cache.active.stress * cp.getKinematicState().getDeformationGradient() *
            cache.fiber.tensor();
      }

      /**
       * @brief Adds the active contribution to the material tangent action.
       * @param dP Storage for the resulting stress increment.
       * @param cache Constitutive invariant cache populated by getCache().
       * @param cp Constitutive evaluation point.
       * @param dF Increment of the deformation gradient.
       */
      void getMaterialTangent(Math::SpatialMatrix<Real>& dP, const Cache& cache,
        const ConstitutivePoint& cp, const Math::SpatialMatrix<Real>& dF) const
      {
        m_passiveLaw.getMaterialTangent(dP, cache.passive, cp, dF);

        // Active fields are independent global unknowns: this block takes the
        // displacement derivative at fixed extension, gamma and beta.
        const Real dStrain = cache.fiber.dStrain(dF);
        const Real fiberTangent =
          cache.dynamic ? cache.active.tangent : cache.active.dStressDe;
        dP = dP + cache.active.stress * dF * cache.fiber.tensor() +
          fiberTangent * dStrain * cp.getKinematicState().getDeformationGradient() *
            cache.fiber.tensor();
      }

    private:
      /**
       * @brief Checks whether the constitutive point contains the active-law time data.
       * @param cp Constitutive point containing the local material data.
       * @returns True when any dynamic time, activation or internal-state input is present.
       */

      static bool hasDynamicData(const ConstitutivePoint& cp)
      {
        return cp.has<Tags::TimeStep>() || cp.has<Tags::PreviousActiveExtension>() ||
          cp.has<Tags::PreviousActiveGamma>() || cp.has<Tags::PreviousActiveBeta>() ||
          cp.has<Tags::ElectricalActivation>() || cp.has<Tags::ActiveGamma>() ||
          cp.has<Tags::ActiveBeta>();
      }

      PassiveLaw m_passiveLaw;
      ActiveLaw m_activeLaw;
  };
}

#endif
