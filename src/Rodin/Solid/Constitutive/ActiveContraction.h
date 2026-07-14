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
        /// @brief Updated active branch internal state.
          typename ActiveLaw::State newState;
        /// @brief Whether the dynamic active update was used.
          bool dynamic = false;
        /// @brief Number of local Newton iterations used.
          size_t localIterations = 0;
      };

      /// @brief Constructs the coupled active contraction law.
      ActiveContraction(
          const PassiveLaw& passiveLaw,
          const ActiveLaw& activeLaw = ActiveLaw())
        : m_passiveLaw(passiveLaw),
          m_activeLaw(activeLaw),
          m_localTolerance(1e-12),
          m_localMaxIterations(50)
      {}

      /// @brief Sets the convergence tolerance for the per-quadrature-point local
      ///        Newton solve on @f$e_c@f$.
      ActiveContraction& setLocalTolerance(Real tol)
      {
        m_localTolerance = tol;
        return *this;
      }

      /// @brief Sets the maximum number of local Newton iterations on @f$e_c@f$.
      ActiveContraction& setLocalMaxIterations(size_t n)
      {
        m_localMaxIterations = n;
        return *this;
      }

      /// @brief Returns the passive hyperelastic law.
      const PassiveLaw& getPassiveLaw() const
      {
        return m_passiveLaw;
      }

      /// @brief Returns the active fiber law.
      const ActiveLaw& getActiveLaw() const
      {
        return m_activeLaw;
      }

      /// @brief Populates cached passive and active quantities.
      void setCache(Cache& cache, const ConstitutivePoint& cp) const
      {
        m_passiveLaw.setCache(cache.passive, cp);

        cache.fiber = FiberKinematics(cp);
        cache.strain = cache.fiber.strain();

        if (hasDynamicData(cp))
        {
          // Per-quadrature-point local Newton on the active extension c^{n+1}.
          //
          // The local equation is the discrete Hill-Maxwell residual
          // R(c; e, c^n, gamma^n, beta^n, dt, u_1) = 0 produced by
          // ActiveLaw::evaluateDynamic.  Once converged, the global tangent is
          // the Schur-condensed quantity Response::tangent.
          typename ActiveLaw::State oldState;
          oldState.gamma = cp.get<Tags::PreviousActiveGamma>();
          oldState.beta = cp.get<Tags::PreviousActiveBeta>();
          const Real dt = cp.get<Tags::TimeStep>();
          const Real ecN = cp.get<Tags::PreviousActiveExtension>();
          const Real activation = cp.get<Tags::ElectricalActivation>();

          // Initial guess: previous extension, optionally overridden by the
          // input via Tags::ActiveExtension (warm start).
          Real c =
            cp.has<Tags::ActiveExtension>() ? cp.get<Tags::ActiveExtension>() : ecN;

          typename ActiveLaw::State newState =
            m_activeLaw.update(dt, oldState, ecN, c, activation);
          typename ActiveLaw::Response resp = m_activeLaw.evaluateDynamic(
            dt, oldState, newState, cache.strain, ecN, c, activation);

          size_t it = 0;
          for (; it < m_localMaxIterations; ++it)
          {
            // Response::residual is R/kcc, i.e. the next Newton update is
            // dc = -residual.  Stop when the update is below tolerance.
            const Real dc = -resp.residual;
            c += dc;
            newState = m_activeLaw.update(dt, oldState, ecN, c, activation);
            resp = m_activeLaw.evaluateDynamic(
              dt, oldState, newState, cache.strain, ecN, c, activation);
            if (std::abs(dc) < m_localTolerance)
              break;
          }
          cache.activeExtension = c;
          cache.active = resp;
          cache.newState = newState;
          cache.dynamic = true;
          cache.localIterations = it;
        }
        else
        {
          cache.activeExtension = cp.has<Tags::ActiveExtension>()
            ? cp.get<Tags::ActiveExtension>()
            : m_activeLaw.getParameters().initial.extension;
          cache.active =
            m_activeLaw.evaluateStatic(cache.strain, cache.activeExtension);
          cache.newState = m_activeLaw.initialState();
          cache.dynamic = false;
          cache.localIterations = 0;
        }
      }

      /// @brief Returns the sum of passive and active strain-energy densities.
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

      /// @brief Adds the active contribution to the first Piola-Kirchhoff stress.
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

      /// @brief Adds the active contribution to the material tangent action.
      void getMaterialTangent(
          Math::SpatialMatrix<Real>& dP,
          const Cache& cache,
          const ConstitutivePoint& cp,
          const Math::SpatialMatrix<Real>& dF) const
      {
        m_passiveLaw.getMaterialTangent(dP, cache.passive, cp, dF);

        // Active contribution: P_a = sigma(e, c) * F * (n_f x n_f).
        //
        //  - Dynamic mode: c was solved locally from R(c; e, ...) = 0, so c
        //    is an implicit function of e and the global tangent picks up the
        //    Schur correction
        //        dsigma/de_total = dsigma/de + dsigma/dc * (-kce/kcc)
        //    that ActiveLaw::evaluateDynamic stores in Response::tangent.
        //
        //  - Static mode: c is supplied externally via Tags::ActiveExtension
        //    and is not eliminated.  The correct material tangent is then the
        //    partial derivative dsigma/de at fixed c, i.e. dStressDe.
        const Real dStrain = cache.fiber.dStrain(dF);
        const Real fiberTangent =
          cache.dynamic ? cache.active.tangent : cache.active.dStressDe;
        dP = dP
           + cache.active.stress * dF * cache.fiber.tensor()
           + fiberTangent * dStrain
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
      Real m_localTolerance;
      size_t m_localMaxIterations;
  };
}

#endif
