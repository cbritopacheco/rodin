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
#include <utility>

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
   * discontinuous fields. Const accessors expose their residual and tangent blocks
   * for a coupled global Newton solve; no local solve or condensation occurs.
   * @tparam PassiveLaw Passive constitutive law with non-const binding and const evaluation.
   * @tparam ActiveLaw Three-field active law providing residual and partial tangent blocks.
   */
  template <class PassiveLaw, class ActiveLaw = ActiveFiberLaw>
  class ActiveContraction final
    : public HyperElasticLaw<ActiveContraction<PassiveLaw, ActiveLaw>>
  {
    public:
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
       * @brief Binds passive and active inputs, reusing unchanged active responses.
       * @param cp Constitutive point supplying kinematics and active fields.
       * @returns This law after successful binding.
       */
      ActiveContraction& setConstitutivePoint(const ConstitutivePoint& cp)
      {
        const bool prepared = std::exchange(m_set, false);
        const auto direction = FiberKinematics::direction(cp);
        const Inputs inputs{read<Tags::TimeStep>(cp),
          read<Tags::PreviousActiveExtension>(cp), read<Tags::PreviousActiveGamma>(cp),
          read<Tags::PreviousActiveBeta>(cp), read<Tags::ElectricalActivation>(cp),
          read<Tags::ActiveExtension>(cp), read<Tags::ActiveGamma>(cp),
          read<Tags::ActiveBeta>(cp), read<Tags::PreviousFiberStrain>(cp)};
        const bool reuse = prepared &&
          this->hasSameKinematics(
            m_passiveLaw.getKinematicState(), cp.getKinematicState()) &&
          this->hasSameDirection(m_direction, direction) && m_inputs == inputs;
        m_passiveLaw.setConstitutivePoint(cp);
        if (!reuse)
        {
          computeCache(m_cache, cp, direction, inputs);
          m_inputs = inputs;
          m_direction = direction;
        }
        m_set = true;
        return *this;
      }

      /**
       * @brief Returns owned kinematics shared with the passive law.
       * @returns Prepared kinematic state.
       */
      const KinematicState& getKinematicState() const
      {
        this->require(m_set, __func__);
        return m_passiveLaw.getKinematicState();
      }

      /**
       * @brief Returns the prepared strain-energy density.
       * @returns Passive and active strain-energy density.
       */
      Real getStrainEnergyDensity() const
      {
        const auto& cache = m_cache;
        this->require(m_set, __func__);
        const Real passiveEnergy = m_passiveLaw.getStrainEnergyDensity();
        const Real denom = 1.0 + 2.0 * cache.activeExtension;
        const Real activeEnergy = 0.5 * m_activeLaw.getParameters().stiffness *
          (cache.strain - cache.activeExtension) *
          (cache.strain - cache.activeExtension) / (denom * denom);
        return passiveEnergy + activeEnergy;
      }

      /**
       * @brief Computes prepared passive and active stress.
       * @param P Output first Piola-Kirchhoff stress.
       */
      void getFirstPiolaKirchhoffStress(Math::SpatialMatrix<Real>& P) const
      {
        const auto& cache = m_cache;
        const auto& state = getKinematicState();
        m_passiveLaw.getFirstPiolaKirchhoffStress(P);
        P =
          P + cache.active.stress * state.getDeformationGradient() * cache.fiber.tensor();
      }

      /**
       * @brief Computes the prepared material tangent action.
       * @param dP Output stress increment.
       * @param dF Deformation-gradient increment.
       */
      void getMaterialTangent(
        Math::SpatialMatrix<Real>& dP, const Math::SpatialMatrix<Real>& dF) const
      {
        const auto& cache = m_cache;
        const auto& state = getKinematicState();
        m_passiveLaw.getMaterialTangent(dP, dF);

        // Active fields are independent global unknowns: this block takes the
        // displacement derivative at fixed extension, gamma and beta.
        const Real dStrain = cache.fiber.dStrain(dF);
        const Real fiberTangent =
          cache.dynamic ? cache.active.tangent : cache.active.dStressDe;
        dP = dP + cache.active.stress * dF * cache.fiber.tensor() +
          fiberTangent * dStrain * state.getDeformationGradient() * cache.fiber.tensor();
      }

      /**
       * @brief Returns fiber kinematics.
       * @returns Prepared fiber kinematics.
       */
      const FiberKinematics& getFiberKinematics() const
      {
        this->require(m_set, __func__);
        return m_cache.fiber;
      }

      /**
       * @brief Returns current fiber strain.
       * @returns Prepared current fiber strain.
       */
      Real getFiberStrain() const
      {
        this->require(m_set, __func__);
        return m_cache.strain;
      }

      /**
       * @brief Returns current active extension.
       * @returns Prepared current active extension.
       */
      Real getActiveExtension() const
      {
        this->require(m_set, __func__);
        return m_cache.activeExtension;
      }

      /**
       * @brief Returns active stress response.
       * @returns Prepared active stress response.
       */
      const typename ActiveLaw::Response& getActiveResponse() const
      {
        this->require(m_set, __func__);
        return m_cache.active;
      }

      /**
       * @brief Returns current independent active state.
       * @returns Prepared current independent active state.
       */
      const typename ActiveLaw::State& getActiveState() const
      {
        this->require(m_set, __func__);
        return m_cache.newState;
      }

      /**
       * @brief Returns whether dynamic active equations were evaluated.
       * @returns Prepared whether dynamic active equations were evaluated.
       */
      bool isDynamic() const
      {
        this->require(m_set, __func__);
        return m_cache.dynamic;
      }

      /**
       * @brief Returns uneliminated active residual and derivatives.
       * @returns Prepared uneliminated active residual and derivatives.
       */
      const typename ActiveLaw::Residual& getActiveResidual() const
      {
        this->require(m_set, __func__);
        return m_cache.residual;
      }

    private:
      /// @brief Optional values preserve both tag presence and consumed inputs.
      struct Inputs
      {
        /// @brief Time step.
        Optional<Real> timeStep;

        /// @brief Previous active extension.
        Optional<Real> previousActiveExtension;

        /// @brief Previous active gamma.
        Optional<Real> previousActiveGamma;

        /// @brief Previous active beta.
        Optional<Real> previousActiveBeta;

        /// @brief Electrical activation.
        Optional<Real> electricalActivation;

        /// @brief Active extension.
        Optional<Real> activeExtension;

        /// @brief Active gamma.
        Optional<Real> activeGamma;

        /// @brief Active beta.
        Optional<Real> activeBeta;

        /// @brief Previous fiber strain for midpoint evaluation.
        Optional<Real> previousFiberStrain;

        /**
         * @brief Compares captured values and their presence exactly.
         * @param other Captured inputs to compare.
         * @return Whether every captured value and its presence agree.
         */
        bool operator==(const Inputs& other) const = default;
      };

      /// @brief Cached active quantities at a quadrature point.
      struct Cache
      {
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
       * @brief Populates cached active quantities.
       * @param cache Storage for the constitutive invariant cache.
       * @param cp Constitutive evaluation point.
       * @param direction Captured normalized reference fiber direction.
       * @param inputs Captured scalar active inputs and tag presence.
       */
      void computeCache(Cache& cache, const ConstitutivePoint& cp,
        const Math::SpatialVector<Real>& direction, const Inputs& inputs) const
      {
        cache.fiber = FiberKinematics(cp.getKinematicState(), direction);
        cache.strain = cache.fiber.strain();

        if (hasDynamicData(inputs))
        {
          if (!inputs.timeStep || !inputs.previousActiveExtension ||
            !inputs.previousActiveGamma || !inputs.previousActiveBeta ||
            !inputs.electricalActivation || !inputs.activeExtension ||
            !inputs.activeGamma || !inputs.activeBeta)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Dynamic contraction requires independent active fields."
              << Alert::Raise;
          }
          typename ActiveLaw::State oldState;
          oldState.gamma = (*inputs.previousActiveGamma);
          oldState.beta = (*inputs.previousActiveBeta);
          cache.activeExtension = (*inputs.activeExtension);
          cache.newState.gamma = (*inputs.activeGamma);
          cache.newState.beta = (*inputs.activeBeta);
          const bool midpoint = inputs.previousFiberStrain.has_value();
          const Real strainFactor = midpoint ? Real(0.5) : Real(1);
          const Real strain = midpoint
            ? Real(0.5) * (cache.strain + (*inputs.previousFiberStrain))
            : cache.strain;
          const Real dt = (*inputs.timeStep);
          const Real previousExtension = (*inputs.previousActiveExtension);
          const Real activation = (*inputs.electricalActivation);
          m_activeLaw.getResidual(cache.residual, dt, oldState, cache.newState, strain,
            previousExtension, cache.activeExtension, activation, strainFactor);
          cache.active = m_activeLaw.getResponse(cache.residual, strain,
            previousExtension, cache.activeExtension, strainFactor);
          cache.dynamic = true;
        }
        else
        {
          cache.activeExtension = inputs.activeExtension.has_value()
            ? (*inputs.activeExtension)
            : m_activeLaw.getParameters().initial.extension;
          cache.active = m_activeLaw.evaluateStatic(cache.strain, cache.activeExtension);
          cache.newState = m_activeLaw.initialState();
          cache.dynamic = false;
          cache.residual = typename ActiveLaw::Residual();
        }
      }

      /**
       * @brief Captures one typed active input without retaining the point.
       * @tparam Tag Active scalar input tag.
       * @param cp Incoming constitutive point.
       * @returns The supplied value, or empty when absent.
       */
      template <class Tag>
      static Optional<Real> read(const ConstitutivePoint& cp)
      {
        if (cp.has<Tag>())
        {
          return cp.get<Tag>();
        }
        return {};
      }

      /// @brief Active response for the captured inputs.
      Cache m_cache{};

      /// @brief Captured scalar inputs and their presence.
      Inputs m_inputs{};

      /// @brief Captured preferred direction before fiber construction.
      Math::SpatialVector<Real> m_direction;

      /// @brief Whether all passive and active preparation succeeded.
      bool m_set = false;

      /**
       * @brief Checks whether the constitutive point contains the active-law time data.
       * @param inputs Captured active inputs and their presence.
       * @returns True when any dynamic time, activation or internal-state input is present.
       */
      static bool hasDynamicData(const Inputs& inputs)
      {
        return inputs.timeStep.has_value() ||
          inputs.previousActiveExtension.has_value() ||
          inputs.previousActiveGamma.has_value() ||
          inputs.previousActiveBeta.has_value() ||
          inputs.electricalActivation.has_value() || inputs.activeGamma.has_value() ||
          inputs.activeBeta.has_value();
      }

      PassiveLaw m_passiveLaw;
      ActiveLaw m_activeLaw;
  };
}

#endif
