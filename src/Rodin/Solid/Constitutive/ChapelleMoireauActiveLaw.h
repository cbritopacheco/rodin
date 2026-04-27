/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ChapelleMoireauActiveLaw.h
 * @brief Chapelle-Moireau active contraction rheology.
 */
#ifndef RODIN_SOLID_CONSTITUTIVE_CHAPELLEMOIREAUACTIVELAW_H
#define RODIN_SOLID_CONSTITUTIVE_CHAPELLEMOIREAUACTIVELAW_H

#include <algorithm>
#include <cmath>

#include "Rodin/Types.h"

namespace Rodin::Solid
{
  /**
   * @brief Local Chapelle-Moireau/Hill-Maxwell active contraction law.
   *
   * This is the Rodin-side counterpart of Felisce's
   * `ChapelleMoireauActiveLaw`: it evaluates the one-dimensional active fiber
   * stress, the local contraction residual, and the Schur-complement
   * derivative used when the active extension @f$e_c@f$ is eliminated locally.
   */
  class ChapelleMoireauActiveLaw
  {
    public:
      struct Input
      {
        Real Es = 1.0;
        Real initFibDef = 0.0;
        Real initActiveStiffness = 0.0;
        Real initActiveStress = 0.0;
        Real DampingParallel = 0.0;
        Real DestructionRate = 0.0;
        Real CrossBridgeStiffness = 0.0;
        Real Contractility = 0.0;
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

      ChapelleMoireauActiveLaw()
        : m_input()
      {}

      explicit ChapelleMoireauActiveLaw(const Input& input)
        : m_input(input)
      {}

      const Input& getInput() const
      {
        return m_input;
      }

      State initialState() const
      {
        State state;
        state.gamma = std::sqrt(std::max<Real>(m_input.initActiveStiffness, 0.0));
        state.beta = state.gamma > 0.0 ? m_input.initActiveStress / state.gamma : 0.0;
        return state;
      }

      Real getStress(Real strain1D, Real activeExtension) const
      {
        const Real denom = 1.0 + 2.0 * activeExtension;
        return m_input.Es * (strain1D - activeExtension) / (denom * denom);
      }

      Real getPartialDerivativeStressWrtStrain1D(Real activeExtension) const
      {
        const Real denom = 1.0 + 2.0 * activeExtension;
        return m_input.Es / (denom * denom);
      }

      Real getPartialDerivativeStressWrtActiveExtension(
          Real strain1D,
          Real activeExtension) const
      {
        const Real denom = 1.0 + 2.0 * activeExtension;
        return m_input.Es * (2.0 * activeExtension - 4.0 * strain1D - 1.0)
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
        response.k22 = m_input.Es * (1.0 + 2.0 * strain1D);
        response.k21 = -m_input.Es * (1.0 + 4.0 * strain1D - 2.0 * activeExtension);
        response.contractionResidualOverK22 =
          (-m_input.Es * (strain1D - activeExtension) * (1.0 + 2.0 * strain1D))
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
        const Real alpha = m_input.DestructionRate;
        const Real k0 = m_input.CrossBridgeStiffness;
        const Real sigma0 = m_input.Contractility;
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
        const Real fib12 = 0.5 * (activeExtension + previousActiveExtension);
        const Real onePlus2Fib12 = 1.0 + 2.0 * fib12;
        const Real activeBranchStress =
          newState.activeStress()
          + m_input.DampingParallel
          * (activeExtension - previousActiveExtension) / dt;

        Response response;
        response.stress = getStress(strain1D, fib12);
        response.dStressWrtStrain =
          getPartialDerivativeStressWrtStrain1D(fib12);
        response.dStressWrtActiveExtension =
          getPartialDerivativeStressWrtActiveExtension(strain1D, fib12);
        response.k21 = -m_input.Es * (1.0 + 4.0 * strain1D - 2.0 * fib12);
        response.k22 =
          3.0 * onePlus2Fib12 * onePlus2Fib12 * activeBranchStress
          + onePlus2Fib12 * onePlus2Fib12 * onePlus2Fib12
            * (getPartialDerivativeActiveStressWrtActiveExtension(
                dt, oldState, activation, previousActiveExtension, activeExtension)
              + m_input.DampingParallel / dt)
          + 0.5 * m_input.Es * (1.0 + 2.0 * strain1D);
        response.contractionResidualOverK22 =
          (activeBranchStress * onePlus2Fib12 * onePlus2Fib12 * onePlus2Fib12
           - m_input.Es * (strain1D - fib12) * (1.0 + 2.0 * strain1D))
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
        const Real alpha = m_input.DestructionRate;
        const Real k0 = m_input.CrossBridgeStiffness;
        const Real sigma0 = m_input.Contractility;
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
      Input m_input;
  };
}

#endif
