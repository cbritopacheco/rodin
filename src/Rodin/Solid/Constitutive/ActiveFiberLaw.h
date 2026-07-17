/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ActiveFiberLaw.h
 * @brief Local active fiber contraction law.
 */
#ifndef RODIN_SOLID_CONSTITUTIVE_ACTIVEFIBERLAW_H
#define RODIN_SOLID_CONSTITUTIVE_ACTIVEFIBERLAW_H

#include <algorithm>
#include <cmath>

#include "Rodin/Types.h"

namespace Rodin::Solid
{
  /**
   * @brief Local active fiber contraction law.
   *
   * Evaluates a one-dimensional fiber stress, the local contraction residual,
   * and the condensed tangent obtained after eliminating the local extension
   * variable @f$e_c@f$.
   */
  class ActiveFiberLaw
  {
    public:
      /// @brief Parameters controlling the local active-fiber law.
      struct Parameters
      {
        /// @brief Initial active extension, stiffness, and stress data.
          struct Initial
          {
          /// @brief Initial active extension.
              Real extension = 0.0;
          /// @brief Initial active branch stiffness.
              Real stiffness = 0.0;
          /// @brief Initial active branch stress.
              Real stress = 0.0;
          };

        /// @brief Passive serial stiffness in the active branch.
          Real stiffness = 1.0;
        /// @brief Viscous damping coefficient for dynamic active extension.
          Real damping = 0.0;
        /// @brief Cross-bridge destruction rate.
          Real destructionRate = 0.0;
        /// @brief Cross-bridge stiffness coefficient.
          Real crossBridgeStiffness = 0.0;
        /// @brief Activation-dependent contractility coefficient.
          Real contractility = 0.0;

        /// @brief Initial active branch state parameters.
          Initial initial;
      };

      /// @brief Internal active-fiber state.
      struct State
      {
        /// @brief Square-root stiffness-like internal variable.
          Real gamma = 0.0;
        /// @brief Stress-like internal variable.
          Real beta = 0.0;

        /// @brief Returns the active stress @f$\gamma\beta@f$.
          Real activeStress() const
          {
            return gamma * beta;
          }
      };

      /// @brief Local response and tangent data for a fiber update.
      struct Response
      {
        /// @brief Active fiber stress.
          Real stress = 0.0;
        /// @brief Partial derivative of stress with respect to fiber strain.
          Real dStressDe = 0.0;
          /// @brief Partial derivative of stress with respect to active extension.
          Real dStressDc = 0.0;
          /// @brief Scaled local residual used as the Newton update.
          Real residual = 0.0;
          /// @brief Mixed derivative of the local residual with respect to strain.
          Real kce = 0.0;
          /// @brief Derivative of the local residual with respect to active extension.
          Real kcc = 1.0;
          /// @brief Condensed tangent after eliminating the active extension.
          Real tangent = 0.0;
      };

      /// @brief Constructs the active law with default parameters.
      ActiveFiberLaw()
        : m_parameters()
      {}

      /// @brief Constructs the active law from parameters.
      explicit ActiveFiberLaw(const Parameters& parameters)
        : m_parameters(parameters)
      {}

      /// @brief Returns the material parameters.
      const Parameters& getParameters() const
      {
        return m_parameters;
      }

      /// @brief Builds the initial internal active-fiber state.
      State initialState() const
      {
        State state;
        state.gamma =
          std::sqrt(std::max<Real>(m_parameters.initial.stiffness, 0.0));
        state.beta = state.gamma > 0.0
          ? m_parameters.initial.stress / state.gamma
          : 0.0;
        return state;
      }

      /// @brief Evaluates active stress at fiber strain @p e and extension @p c.
      Real stress(Real e, Real c) const
      {
        const Real denom = 1.0 + 2.0 * c;
        return m_parameters.stiffness
             * (e - c) / (denom * denom);
      }

      /// @brief Evaluates @f$\partial\sigma/\partial e@f$ at fixed extension.
      Real dStressDe(Real c) const
      {
        const Real denom = 1.0 + 2.0 * c;
        return m_parameters.stiffness / (denom * denom);
      }

      /// @brief Evaluates @f$\partial\sigma/\partial c@f$ at fixed fiber strain.
      Real dStressDc(Real e, Real c) const
      {
        const Real denom = 1.0 + 2.0 * c;
        return m_parameters.stiffness
             * (2.0 * c - 4.0 * e - 1.0)
             / (denom * denom * denom);
      }

      /// @brief Evaluates the static active response.
      Response evaluateStatic(Real e, Real c) const
      {
        Response response;
        response.stress = stress(e, c);
        response.dStressDe = dStressDe(c);
        response.dStressDc = dStressDc(e, c);
        response.kcc = m_parameters.stiffness * (1.0 + 2.0 * e);
        response.kce =
          -m_parameters.stiffness * (1.0 + 4.0 * e - 2.0 * c);
        response.residual =
          (-m_parameters.stiffness
           * (e - c) * (1.0 + 2.0 * e))
          / response.kcc;
        response.tangent =
          response.dStressDe - response.dStressDc * response.kce / response.kcc;
        return response;
      }

      /// @brief Advances the internal active-fiber state.
      State update(
          Real dt,
          const State& oldState,
          Real previousActiveExtension,
          Real activeExtension,
          Real activation) const
      {
        const Real alpha = m_parameters.destructionRate;
        const Real k0 = m_parameters.crossBridgeStiffness;
        const Real sigma0 = m_parameters.contractility;

        const Real delta = activeExtension - previousActiveExtension;
        const Real activationPlus = std::max<Real>(activation, 0.0);
        const Real n0 = starling(previousActiveExtension);

        const Real denominatorGamma =
          1.0
          + dt * std::abs(activation)
          + alpha * std::abs(delta);

        const Real gammaSquare =
          std::max<Real>(1.e-16,
              (oldState.gamma * oldState.gamma + dt * n0 * k0 * activationPlus)
              / denominatorGamma);

        State state;
        state.gamma = std::sqrt(gammaSquare);

        const Real denominatorBeta =
          1.0
          + 0.5 * dt * n0 * k0 * activationPlus / gammaSquare
          + 0.5 * dt * std::abs(activation)
          + 0.5 * alpha * std::abs(delta);

        state.beta =
          (oldState.beta + state.gamma * delta
           + dt * n0 * sigma0 * activationPlus / state.gamma)
          / denominatorBeta;

        return state;
      }

      /// @brief Evaluates the dynamic active response and condensed tangent.
      ///
      /// @param dt Time step @f$\Delta t@f$.
      /// @param oldState Previous internal state @f$(\gamma^n, \beta^n)@f$.
      /// @param newState Updated internal state @f$(\gamma^{n+1}, \beta^{n+1})@f$.
      /// @param e Fiber strain at which the series law is evaluated. For the
      ///   compatible discretization this is the midpoint strain
      ///   @f$e_{1D}^{n+\frac{1}{2}}@f$.
      /// @param previousActiveExtension Previous active extension @f$e_c^n@f$.
      /// @param activeExtension Current active extension @f$e_c^{n+1}@f$.
      /// @param activation Electrical activation @f$u_1@f$.
      /// @param strainFactor Derivative @f$\partial e/\partial e_{1D}^{n+1}@f$
      ///   of the evaluation strain with respect to the current fiber strain.
      ///   It is @f$\frac{1}{2}@f$ for the midpoint strain and @f$1@f$ when
      ///   @p e is the current strain. Only the condensed tangent depends on
      ///   it, since the global tangent differentiates with respect to
      ///   @f$e_{1D}^{n+1}@f$.
      Response evaluateDynamic(Real dt, const State& oldState, const State& newState,
        Real e, Real previousActiveExtension, Real activeExtension, Real activation,
        Real strainFactor = 1.0) const
      {
        const Real midpointExtension =
          0.5 * (activeExtension + previousActiveExtension);
        const Real onePlus2MidpointExtension = 1.0 + 2.0 * midpointExtension;
        const Real activeBranchStress =
          newState.activeStress()
          + m_parameters.damping
          * (activeExtension - previousActiveExtension) / dt;

        Response response;
        response.stress = stress(e, midpointExtension);
        response.dStressDe = dStressDe(midpointExtension);
        response.dStressDc = dStressDc(e, midpointExtension);
        response.kce =
          -m_parameters.stiffness
          * (1.0 + 4.0 * e - 2.0 * midpointExtension);
        response.kcc =
          3.0 * onePlus2MidpointExtension * onePlus2MidpointExtension
            * activeBranchStress
          + onePlus2MidpointExtension * onePlus2MidpointExtension
            * onePlus2MidpointExtension
            * (dActiveStressDc(
                dt, oldState, activation, previousActiveExtension, activeExtension)
              + m_parameters.damping / dt)
          + 0.5 * m_parameters.stiffness * (1.0 + 2.0 * e);
        response.residual =
          (activeBranchStress
             * onePlus2MidpointExtension * onePlus2MidpointExtension
             * onePlus2MidpointExtension
           - m_parameters.stiffness
             * (e - midpointExtension) * (1.0 + 2.0 * e))
          / response.kcc;
        // Chain rule to the current fiber strain: sigma and the local residual
        // are functions of the evaluation strain e, whose derivative with
        // respect to e_1D^{n+1} is strainFactor.
        response.tangent = strainFactor *
          (response.dStressDe - 0.5 * response.dStressDc * response.kce / response.kcc);
        return response;
      }

      /// @brief Evaluates the derivative of active stress with respect to extension.
      Real dActiveStressDc(
          Real dt,
          const State& oldState,
          Real activation,
          Real previousActiveExtension,
          Real activeExtension) const
      {
        const Real alpha = m_parameters.destructionRate;
        const Real k0 = m_parameters.crossBridgeStiffness;
        const Real sigma0 = m_parameters.contractility;

        const Real delta = activeExtension - previousActiveExtension;
        const Real absDelta = std::abs(delta);
        const Real sign = delta > 0.0 ? 1.0 : (delta < 0.0 ? -1.0 : 0.0);
        const Real activationPlus = std::max<Real>(activation, 0.0);
        const Real n0 = starling(previousActiveExtension);

        const Real Dg =
          1.0
          + dt * std::abs(activation)
          + alpha * absDelta;

        const Real Ng =
          oldState.gamma * oldState.gamma
          + dt * n0 * k0 * activationPlus;

        const Real gammaSquare = std::max<Real>(1.e-16, Ng / Dg);
        const Real gamma = std::sqrt(gammaSquare);

        const Real dGammaSquare = -Ng * alpha * sign / (Dg * Dg);
        const Real dGamma = 0.5 * dGammaSquare / gamma;

        const Real Nb =
          oldState.beta
          + gamma * delta
          + dt * n0 * sigma0 * activationPlus / gamma;

        const Real Db =
          1.0
          + 0.5 * dt * n0 * k0 * activationPlus / gammaSquare
          + 0.5 * dt * std::abs(activation)
          + 0.5 * alpha * absDelta;

        const Real dNb =
          dGamma * delta
          + gamma
          - dt * n0 * sigma0 * activationPlus * dGamma / (gamma * gamma);

        const Real dDb =
          -0.5 * dt * n0 * k0 * activationPlus
            * dGammaSquare / (gammaSquare * gammaSquare)
          + 0.5 * alpha * sign;

        const Real dBeta = (dNb * Db - Nb * dDb) / (Db * Db);

        return dGamma * (Nb / Db) + gamma * dBeta;
      }

      /// @brief Evaluates the length-dependent Starling activation factor.
      static Real starling(Real activeExtension)
      {
        const Real x1 = -0.4;
        const Real y1 = 0.0;
        const Real x2 = 0.3;
        const Real y2 = 0.38;
        const Real x3 = 0.73;
        const Real y3 = 0.74;
        const Real x4 = 1.0;
        const Real y4 = 1.0;
        const Real x5 = 1.3;
        const Real y5 = 1.0;
        const Real x6 = 2.4;
        const Real y6 = 0.0;

        Real n0 = 0.0;
        if (activeExtension < x2)
          n0 = ((y2 - y1) / (x2 - x1)) * (activeExtension - x2) + y2;
        else if (activeExtension < x3)
          n0 = ((y3 - y2) / (x3 - x2)) * (activeExtension - x3) + y3;
        else if (activeExtension < x4)
          n0 = ((y4 - y3) / (x4 - x3)) * (activeExtension - x4) + y4;
        else if (activeExtension < x5)
          n0 = y4;
        else if (activeExtension < x6)
          n0 = ((y6 - y5) / (x6 - x5)) * (activeExtension - x6) + y6;
        return std::max<Real>(n0, 0.0);
      }

    private:
      Parameters m_parameters;
  };
}

#endif
