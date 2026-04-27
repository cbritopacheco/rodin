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
      struct Parameters
      {
        Real stiffness = 1.0;
        Real initialExtension = 0.0;
        Real initialActiveStiffness = 0.0;
        Real initialActiveStress = 0.0;
        Real damping = 0.0;
        Real destructionRate = 0.0;
        Real crossBridgeStiffness = 0.0;
        Real contractility = 0.0;
      };

      struct State
      {
        Real gamma = 0.0;
        Real beta = 0.0;

        Real activeStress() const
        {
          return gamma * beta;
        }
      };

      struct Response
      {
        Real stress = 0.0;
        Real dStressWrtStrain = 0.0;
        Real dStressWrtActiveExtension = 0.0;
        Real contractionResidualOverK22 = 0.0;
        Real k21 = 0.0;
        Real k22 = 1.0;
        Real condensedStressTangent = 0.0;
      };

      ActiveFiberLaw()
        : m_parameters()
      {}

      explicit ActiveFiberLaw(const Parameters& parameters)
        : m_parameters(parameters)
      {}

      const Parameters& getParameters() const
      {
        return m_parameters;
      }

      State initialState() const
      {
        State state;
        state.gamma =
          std::sqrt(std::max<Real>(m_parameters.initialActiveStiffness, 0.0));
        state.beta = state.gamma > 0.0
          ? m_parameters.initialActiveStress / state.gamma
          : 0.0;
        return state;
      }

      Real getStress(Real strain1D, Real activeExtension) const
      {
        const Real denom = 1.0 + 2.0 * activeExtension;
        return m_parameters.stiffness
             * (strain1D - activeExtension) / (denom * denom);
      }

      Real getPartialDerivativeStressWrtStrain1D(Real activeExtension) const
      {
        const Real denom = 1.0 + 2.0 * activeExtension;
        return m_parameters.stiffness / (denom * denom);
      }

      Real getPartialDerivativeStressWrtActiveExtension(
          Real strain1D,
          Real activeExtension) const
      {
        const Real denom = 1.0 + 2.0 * activeExtension;
        return m_parameters.stiffness
             * (2.0 * activeExtension - 4.0 * strain1D - 1.0)
             / (denom * denom * denom);
      }

      Response evaluateStatic(Real strain1D, Real activeExtension) const
      {
        Response response;
        response.stress = getStress(strain1D, activeExtension);
        response.dStressWrtStrain =
          getPartialDerivativeStressWrtStrain1D(activeExtension);
        response.dStressWrtActiveExtension =
          getPartialDerivativeStressWrtActiveExtension(strain1D, activeExtension);
        response.k22 = m_parameters.stiffness * (1.0 + 2.0 * strain1D);
        response.k21 =
          -m_parameters.stiffness * (1.0 + 4.0 * strain1D - 2.0 * activeExtension);
        response.contractionResidualOverK22 =
          (-m_parameters.stiffness
           * (strain1D - activeExtension) * (1.0 + 2.0 * strain1D))
          / response.k22;
        response.condensedStressTangent =
          response.dStressWrtStrain
          - response.dStressWrtActiveExtension * response.k21 / response.k22;
        return response;
      }

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
          1.0 + dt * std::abs(activation) + alpha * std::abs(delta);
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

      Response evaluateDynamic(
          Real dt,
          const State& oldState,
          const State& newState,
          Real strain1D,
          Real previousActiveExtension,
          Real activeExtension,
          Real activation) const
      {
        const Real midpointExtension =
          0.5 * (activeExtension + previousActiveExtension);
        const Real onePlus2MidpointExtension = 1.0 + 2.0 * midpointExtension;
        const Real activeBranchStress =
          newState.activeStress()
          + m_parameters.damping
          * (activeExtension - previousActiveExtension) / dt;

        Response response;
        response.stress = getStress(strain1D, midpointExtension);
        response.dStressWrtStrain =
          getPartialDerivativeStressWrtStrain1D(midpointExtension);
        response.dStressWrtActiveExtension =
          getPartialDerivativeStressWrtActiveExtension(strain1D, midpointExtension);
        response.k21 =
          -m_parameters.stiffness
          * (1.0 + 4.0 * strain1D - 2.0 * midpointExtension);
        response.k22 =
          3.0 * onePlus2MidpointExtension * onePlus2MidpointExtension
            * activeBranchStress
          + onePlus2MidpointExtension * onePlus2MidpointExtension
            * onePlus2MidpointExtension
            * (getPartialDerivativeActiveStressWrtActiveExtension(
                dt, oldState, activation, previousActiveExtension, activeExtension)
              + m_parameters.damping / dt)
          + 0.5 * m_parameters.stiffness * (1.0 + 2.0 * strain1D);
        response.contractionResidualOverK22 =
          (activeBranchStress
             * onePlus2MidpointExtension * onePlus2MidpointExtension
             * onePlus2MidpointExtension
           - m_parameters.stiffness
             * (strain1D - midpointExtension) * (1.0 + 2.0 * strain1D))
          / response.k22;
        response.condensedStressTangent =
          response.dStressWrtStrain
          - 0.5 * response.dStressWrtActiveExtension * response.k21 / response.k22;
        return response;
      }

      Real getPartialDerivativeActiveStressWrtActiveExtension(
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

        const Real Dg = 1.0 + dt * std::abs(activation) + alpha * absDelta;
        const Real Ng = oldState.gamma * oldState.gamma
                      + dt * n0 * k0 * activationPlus;
        const Real gammaSquare = std::max<Real>(1.e-16, Ng / Dg);
        const Real gamma = std::sqrt(gammaSquare);

        const Real dGammaSquare = -Ng * alpha * sign / (Dg * Dg);
        const Real dGamma = 0.5 * dGammaSquare / gamma;

        const Real Nb = oldState.beta + gamma * delta
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
