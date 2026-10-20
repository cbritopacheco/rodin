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
#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Alert/Raise.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"

namespace Rodin::Solid
{
  /**
   * @brief Local active fiber contraction law.
   *
   * Evaluates a one-dimensional fiber stress, the local contraction residual,
   * and partial tangent blocks for the independent fields @f$(e_c,\gamma,\beta)@f$.
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

          /**
           * @brief Returns the active stress @f$\gamma\beta@f$.
           * @returns The active stress @f$\gamma\beta@f$.
           */
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

        /// @brief Uneliminated extension residual.
          Real residual = 0.0;
          /// @brief Mixed derivative of the local residual with respect to strain.
          Real kce = 0.0;
          /// @brief Derivative of the local residual with respect to active extension.
          Real kcc = 1.0;

        /// @brief Strain tangent at fixed independent active fields.
          Real tangent = 0.0;
      };

      /**
       * @brief Residual and tangent for the three independent active fields.
       *
       * Rows and columns are ordered as extension, gamma and beta. The strain
       * derivative supplies the coupling to the displacement field. No internal
       * variable is solved for or eliminated by this numerical evaluation.
       */
      struct Residual
      {
        Math::FixedSizeVector<Real, 3> value; ///< Active field residuals.
        Math::FixedSizeMatrix<Real, 3, 3>
          derivative; ///< Derivatives with respect to active fields.
        Math::FixedSizeVector<Real, 3>
          strainDerivative; ///< Derivatives with respect to current fiber strain.

        /// @brief Constructs zero residual and tangent storage.
        Residual()
          : value(),
            derivative(),
            strainDerivative()
        {
          value.setZero();
          derivative.setZero();
          strainDerivative.setZero();
        }
      };

      /**
       * @brief Evaluates the coupled active-field equations and all first derivatives.
       *
       * @f$R_c=d^3(\gamma\beta+\eta(c-c_n)/\Delta t)
       * -k(e-(c+c_n)/2)(1+2e)@f$, where @f$d=1+c+c_n@f$.
       * Gamma and beta satisfy the uneliminated discrete Hill-Maxwell equations.
       * @param output Residual and tangent output.
       * @param dt Positive time-step size.
       * @param oldState Previous gamma and beta fields.
       * @param state Current independent gamma and beta fields.
       * @param e Evaluation fiber strain.
       * @param previousExtension Previous extension field.
       * @param extension Current independent extension field.
       * @param activation Electrical activation.
       * @param strainFactor Derivative of evaluation strain with respect to current strain.
       */
      void getResidual(Residual& output, Real dt, const State& oldState,
        const State& state, Real e, Real previousExtension, Real extension,
        Real activation, Real strainFactor = 1) const
      {
        if (!std::isfinite(dt) || dt <= 0 || !std::isfinite(state.gamma) ||
          state.gamma <= 0 || !std::isfinite(state.beta) ||
          !std::isfinite(oldState.gamma) || !std::isfinite(oldState.beta) ||
          !std::isfinite(e) || !std::isfinite(previousExtension) ||
          !std::isfinite(extension) || !std::isfinite(activation) ||
          !std::isfinite(strainFactor) || 1 + extension + previousExtension <= 0)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Active fields require finite data, positive gamma and time step, "
            << "and positive midpoint stretch." << Alert::Raise;
        }
        const Real delta = extension - previousExtension;
        const Real sign = delta > 0 ? Real(1) : (delta < 0 ? Real(-1) : Real(0));
        const Real alpha = m_parameters.destructionRate;
        const Real source =
          dt * starling(previousExtension) * std::max(activation, Real(0));
        const Real stiffness = source * m_parameters.crossBridgeStiffness;
        const Real contractility = source * m_parameters.contractility;
        const Real gamma = state.gamma;
        const Real beta = state.beta;
        const Real gammaSquared = gamma * gamma;
        const Real gammaDenominator =
          1 + dt * std::abs(activation) + alpha * std::abs(delta);
        const Real betaDenominator = 1 + Real(0.5) * stiffness / gammaSquared +
          Real(0.5) * dt * std::abs(activation) + Real(0.5) * alpha * std::abs(delta);
        const Real midpoint = Real(0.5) * (extension + previousExtension);
        const Real d = 1 + 2 * midpoint;
        const Real branchStress = gamma * beta + m_parameters.damping * delta / dt;
        output.value[0] = d * d * d * branchStress -
          m_parameters.stiffness * (e - midpoint) * (1 + 2 * e);
        const Real gammaNumerator = oldState.gamma * oldState.gamma + stiffness;
        const bool gammaFloor = gammaNumerator / gammaDenominator <= MinimumGammaSquared;
        output.value[1] =
          gammaSquared - std::max(MinimumGammaSquared, gammaNumerator / gammaDenominator);
        output.value[2] =
          betaDenominator * beta - oldState.beta - gamma * delta - contractility / gamma;
        output.derivative.setZero();
        output.derivative(0, 0) = 3 * d * d * branchStress +
          d * d * d * m_parameters.damping / dt +
          Real(0.5) * m_parameters.stiffness * (1 + 2 * e);
        output.derivative(0, 1) = d * d * d * beta;
        output.derivative(0, 2) = d * d * d * gamma;
        output.derivative(1, 0) = gammaFloor
          ? Real(0)
          : gammaNumerator * alpha * sign / (gammaDenominator * gammaDenominator);
        output.derivative(1, 1) = 2 * gamma;
        output.derivative(2, 0) = Real(0.5) * alpha * sign * beta - gamma;
        output.derivative(2, 1) = -stiffness * beta / (gammaSquared * gamma) - delta +
          contractility / gammaSquared;
        output.derivative(2, 2) = betaDenominator;
        output.strainDerivative.setZero();
        output.strainDerivative[0] =
          -strainFactor * m_parameters.stiffness * (1 + 4 * e - 2 * midpoint);
      }

      /**
       * @brief Gets the positive gamma floor used by the discrete state equation.
       * @returns Minimum gamma compatible with the regularized update.
       */
      static Real getMinimumGamma() noexcept
      {
        return std::sqrt(MinimumGammaSquared);
      }

      /// @brief Constructs the active law with default parameters.
      ActiveFiberLaw()
        : m_parameters()
      {}

      /**
       * @brief Constructs the active law from parameters.
       * @param parameters Parameters configuring the operation.
       */
      explicit ActiveFiberLaw(const Parameters& parameters)
        : m_parameters(parameters)
      {}

      /**
       * @brief Returns the material parameters.
       * @returns The material parameters.
       */
      const Parameters& getParameters() const
      {
        return m_parameters;
      }

      /**
       * @brief Builds the initial internal active-fiber state.
       * @returns Initial active-fiber internal state.
       */
      State initialState() const
      {
        State state;
        state.gamma = std::sqrt(std::max<Real>(m_parameters.initial.stiffness, 0.0));
        state.beta = state.gamma > 0.0 ? m_parameters.initial.stress / state.gamma : 0.0;
        return state;
      }

      /**
       * @brief Evaluates active stress at fiber strain @p e and extension @p c.
       * @param e Muscle strain.
       * @param c Contractile state.
       * @returns Active stress at the supplied strain and extension.
       */
      Real stress(Real e, Real c) const
      {
        const Real denom = 1.0 + 2.0 * c;
        return m_parameters.stiffness * (e - c) / (denom * denom);
      }

      /**
       * @brief Evaluates @f$\partial\sigma/\partial e@f$ at fixed extension.
       * @param c Contractile state.
       * @returns Derivative of active stress with respect to strain at fixed extension.
       */
      Real dStressDe(Real c) const
      {
        const Real denom = 1.0 + 2.0 * c;
        return m_parameters.stiffness / (denom * denom);
      }

      /**
       * @brief Evaluates @f$\partial\sigma/\partial c@f$ at fixed fiber strain.
       * @param e Muscle strain.
       * @param c Contractile state.
       * @returns Derivative of active stress with respect to extension at fixed strain.
       */
      Real dStressDc(Real e, Real c) const
      {
        const Real denom = 1.0 + 2.0 * c;
        return m_parameters.stiffness * (2.0 * c - 4.0 * e - 1.0) /
          (denom * denom * denom);
      }

      /**
       * @brief Evaluates the static active response.
       * @param e Muscle strain.
       * @param c Contractile state.
       * @returns Static active stress and tangent response.
       */
      Response evaluateStatic(Real e, Real c) const
      {
        Response response;
        response.stress = stress(e, c);
        response.dStressDe = dStressDe(c);
        response.dStressDc = dStressDc(e, c);
        response.kcc = m_parameters.stiffness * (1.0 + 2.0 * e);
        response.kce = -m_parameters.stiffness * (1.0 + 4.0 * e - 2.0 * c);
        response.residual =
          (-m_parameters.stiffness * (e - c) * (1.0 + 2.0 * e)) / response.kcc;
        response.tangent =
          response.dStressDe - response.dStressDc * response.kce / response.kcc;
        return response;
      }

      /**
       * @brief Advances the internal active-fiber state.
       * @param dt Time-step size.
       * @param oldState State at the previous time step.
       * @param previousActiveExtension Active extension at the previous time step.
       * @param activeExtension Active extension field.
       * @param activation Activation at the current time step.
       * @returns Updated active-fiber internal state.
       */
      State update(Real dt, const State& oldState, Real previousActiveExtension,
        Real activeExtension, Real activation) const
      {
        const Real alpha = m_parameters.destructionRate;
        const Real k0 = m_parameters.crossBridgeStiffness;
        const Real sigma0 = m_parameters.contractility;

        const Real delta = activeExtension - previousActiveExtension;
        const Real activationPlus = std::max<Real>(activation, 0.0);
        const Real n0 = starling(previousActiveExtension);

        const Real denominatorGamma =
          1.0 + dt * std::abs(activation) + alpha * std::abs(delta);

        const Real gammaSquare = std::max<Real>(MinimumGammaSquared,
          (oldState.gamma * oldState.gamma + dt * n0 * k0 * activationPlus) /
            denominatorGamma);

        State state;
        state.gamma = std::sqrt(gammaSquare);

        const Real denominatorBeta = 1.0 +
          0.5 * dt * n0 * k0 * activationPlus / gammaSquare +
          0.5 * dt * std::abs(activation) + 0.5 * alpha * std::abs(delta);

        state.beta = (oldState.beta + state.gamma * delta +
                       dt * n0 * sigma0 * activationPlus / state.gamma) /
          denominatorBeta;

        return state;
      }

      /**
       * @brief Gets fiber stress and partial tangents from prepared active equations.
       * @param residual Active residual and tangent evaluated at the same inputs.
       * @param e Evaluation fiber strain.
       * @param previousActiveExtension Previous independent extension.
       * @param activeExtension Current independent extension.
       * @param strainFactor Derivative of evaluation strain with respect to current strain.
       * @returns Fiber response with no internal-variable elimination.
       */
      Response getResponse(const Residual& residual, Real e, Real previousActiveExtension,
        Real activeExtension, Real strainFactor = 1) const
      {
        const Real midpointExtension =
          Real(0.5) * (activeExtension + previousActiveExtension);
        Response response;
        response.stress = stress(e, midpointExtension);
        response.dStressDe = dStressDe(midpointExtension);
        response.dStressDc = dStressDc(e, midpointExtension);
        response.residual = residual.value[0];
        response.kce = residual.strainDerivative[0];
        response.kcc = residual.derivative(0, 0);
        response.tangent = strainFactor * response.dStressDe;
        return response;
      }

      /**
       * @brief Evaluates the length-dependent Starling activation factor.
       * @param activeExtension Active extension field.
       * @returns Length-dependent activation factor.
       */
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
      /// @brief Absolute floor on gamma squared before division; legacy numerical clamp in state units.
      static constexpr Real MinimumGammaSquared = 1.e-16;

      Parameters m_parameters;
  };
}

#endif
