/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Windkessel.h
 * @brief Windkessel outflow evaluator for the arterial branch.
 *
 * Implements the three-element Windkessel outflow term
 * @f$ Q_{ar} = K_{ar} (p_v - p_{ar}) @f$ when the aortic valve is open
 * (Caruel et al. 2014, §4, arterial outflow model).
 */
#ifndef RODIN_HEART_CCMLC2014_PHYSICS_WINDKESSEL_H
#define RODIN_HEART_CCMLC2014_PHYSICS_WINDKESSEL_H

namespace Rodin::Heart::CCMLC2014::Physics
{

  namespace Rheology {
    struct Newtonian {
    template <typename Scalar, typename EvalData, typename Input>
      static void computeFlow(EvalData& data, const Input& input) {
        const Scalar aorticFlowRate = input.Kar * (data.pvMid - data.parMid);

        data.windkesselOutflow = (aorticFlowRate > Scalar(0)) ? aorticFlowRate : Scalar(0);

        const Scalar midpointDerivative = (aorticFlowRate > Scalar(0))
                                                  ? (Scalar(0.5) * input.Kar)
                                                  : Scalar(0);

        data.dWindkesselOutflow_dPv = midpointDerivative;
        data.dWindkesselOutflow_dPar = -midpointDerivative;

        data.windkesselflowP = (data.parMid - data.pdMid) / input.Rp
        - data.windkesselOutflow;

        data.windkesselflowD = (data.pdMid - data.parMid) / input.Rp
        - (data.pSvMid - data.pdMid) / input.Rd;

        data.dWindkesselflowP_dPar = Scalar(1) / (Scalar(2) * input.Rp)
            - data.dWindkesselOutflow_dPar;

        data.dWindkesselflowP_dPd = -Scalar(1) / (Scalar(2) * input.Rp);

        data.dWindkesselflowD_dPar = -Scalar(1) / (Scalar(2) * input.Rp);

        data.dWindkesselflowD_dPd = Scalar(1) / (Scalar(2) * input.Rp)
            + Scalar(1) / (Scalar(2) * input.Rd);
      }
    };

    struct NonNewtonian {
      struct CY {
        /// @brief Pressure-drop threshold for the Poiseuille fallback.
        static constexpr Real pressureDropTolerance = 1.0e-12;
        /// @brief Minimum shear-rate bracket.
        static constexpr Real minShearRate = 1.0e-8;
        /// @brief Number of RK4 substeps for the WRMS flow integral.
        static constexpr int integralSteps = 100;
        /// @brief Maximum bracketing expansions for outlet scalar solves.
        static constexpr int maxBracketIterations = 100;
        /// @brief Wall shear root solver absolute tolerance.
        static constexpr Real shearAbsoluteTolerance = 1.0e-12;
        /// @brief Wall shear root solver relative tolerance.
        static constexpr Real shearRelativeTolerance = 1.0e-10;
        /// @brief Wall shear root solver step tolerance.
        static constexpr Real shearStepTolerance = 1.0e-12;
        /// @brief Wall shear root solver maximum iterations.
        static constexpr int shearMaxIterations = 50;
        /// @brief Flow inversion root solver absolute tolerance.
        static constexpr Real flowAbsoluteTolerance = 1.0e-10;
        /// @brief Flow inversion root solver relative tolerance.
        static constexpr Real flowRelativeTolerance = 1.0e-9;
        /// @brief Flow inversion root solver step tolerance.
        static constexpr Real flowStepTolerance = 1.0e-12;
        /// @brief Flow inversion root solver maximum iterations.
        static constexpr int flowMaxIterations = 50;
        /// @brief Flow magnitude treated as zero in pressure-drop inversion.
        static constexpr Real zeroFlowTolerance = 1.0e-16;
        /// @brief Minimum pressure-drop bracket.
        static constexpr Real pressureDropBracketMin = 1.0;
        /// @brief Distal capacitor bracket pressure pad.
        static constexpr Real distalPressureBracketPad = 1000.0;

        template <typename Scalar, typename EvalData, typename Input>
        static void computeFlow(EvalData& data, const Input& input) {
            const Scalar aorticFlowRate = input.Kar * (data.pvMid - data.parMid);

            data.windkesselOutflow = (aorticFlowRate > Scalar(0)) ? aorticFlowRate : Scalar(0);

            const Scalar midpointDerivative = (aorticFlowRate > Scalar(0))
                                                      ? (Scalar(0.5) * input.Kar)
                                                      : Scalar(0);

            data.dWindkesselOutflow_dPv = midpointDerivative;
            data.dWindkesselOutflow_dPar = -midpointDerivative;


            auto flowLaw = [&](Real dp, Real L, Real radius) -> std::pair<Real, Real>
            {
              const Real mu0    = input.mu_0;
              const Real muInf  = input.mu_Inf;
              const Real lambda = input.lambda;
              const Real n      = input.n;
              const Real yasuda = input.yasuda;
              const Real delta  = mu0 - muInf;

              const Real sgn = (dp >= 0.0) ? 1.0 : -1.0;
              const Real adp = std::abs(dp);

              const Real R0 =
                8.0 * mu0 * L /
                (std::numbers::pi_v<Real> * std::pow(radius, 4.0));

              if (adp < pressureDropTolerance)
                return {dp / R0, 1.0 / R0};

              const Real tauW = radius * adp / (2.0 * L);

              auto mu = [&](Real g) -> Real
              {
                return muInf
                  + delta * std::pow(
                      1.0 + std::pow(lambda * g, yasuda),
                      (n - 1.0) / yasuda);
              };

              auto dmu = [&](Real g) -> Real
              {
                const Real base = 1.0 + std::pow(lambda * g, yasuda);

                return delta * (n - 1.0)
                  * std::pow(base, (n - 1.0 - yasuda) / yasuda)
                  * std::pow(lambda, yasuda)
                  * std::pow(g, yasuda - 1.0);
              };

              auto tauMinusTauW = [&](Real g) -> std::pair<Real, Real>
              {
                const Real m  = mu(g);
                const Real dm = dmu(g);
                return {g * m - tauW, m + g * dm};
              };

              Math::RootFinding::NewtonRaphson<Real> rootFinder(
                shearAbsoluteTolerance,
                shearRelativeTolerance,
                shearStepTolerance,
                shearMaxIterations);

              Real gHi = std::max<Real>(tauW / muInf, minShearRate);

              for (int k = 0;
                   k < maxBracketIterations && tauMinusTauW(gHi).first < 0.0;
                   ++k)
                gHi *= 2.0;

              if (tauMinusTauW(gHi).first < 0.0)
              {
                std::cerr << "Warning: failed to bracket wall shear rate. "
                          << "Using Poiseuille fallback.\n";
                return {dp / R0, 1.0 / R0};
              }

              const auto gammaRoot =
                rootFinder.solve(tauMinusTauW, 0.5 * gHi, shearStepTolerance, gHi);

              if (!gammaRoot)
              {
                std::cerr << "Warning: failed to solve wall shear rate. "
                          << "Using Poiseuille fallback.\n";
                return {dp / R0, 1.0 / R0};
              }

              const Real gammaW = *gammaRoot;

              auto integrand = [&](Real g) -> Real
              {
                if (g <= 0.0)
                  return 0.0;

                const Real m     = mu(g);
                const Real dm    = dmu(g);
                const Real dtau  = m + g * dm;

                return std::pow(g, 3.0) * m * m * dtau;
              };

              Math::RungeKutta::RK4 integrator;

              const int steps = integralSteps;
              const Real h = gammaW / static_cast<Real>(steps);

              Real I = 0.0;

              auto rhs = [&](Real g, Real y) -> Real
              {
                (void) y;
                return integrand(g);
              };

              for (int i = 0; i < steps; ++i)
              {
                const Real g = static_cast<Real>(i) * h;
                integrator.step(I, g, h, I, rhs);
              }

              if (I <= 0.0 || !std::isfinite(I))
              {
                std::cerr << "Warning: invalid WRMS integral. "
                          << "Using Poiseuille fallback.\n";
                return {dp / R0, 1.0 / R0};
              }

              const Real qAbs =
                std::numbers::pi_v<Real> * std::pow(radius, 3.0) * I /
                std::pow(tauW, 3.0);

              const Real dqAbs =
                (std::numbers::pi_v<Real> * std::pow(radius, 3.0) * gammaW
                 - 3.0 * qAbs) / adp;

              if (!std::isfinite(qAbs) || !std::isfinite(dqAbs) || dqAbs <= 0.0)
              {
                std::cerr << "Warning: invalid WRMS flow derivative. "
                          << "Using Poiseuille fallback.\n";
                return {dp / R0, 1.0 / R0};
              }

              return {sgn * qAbs, dqAbs};
            };

            const Real radiusP = input.proximalRadius;
            const Real lengthP = input.proximalLength;
            const Real radiusD = input.distalRadius;
            const Real lengthD = input.distalLength;

            const Scalar dp_ar_d = data.parMid - data.pdMid;
            const auto [qp, dqp] = flowLaw(dp_ar_d, lengthP, radiusP);

            const Scalar dp_d_sv = data.pSvMid - data.pdMid;
            const auto [qd, dqd] = flowLaw(dp_d_sv, lengthD, radiusD);

            data.windkesselflowP = qp - data.windkesselOutflow;

            data.windkesselflowD = - qp - qd;

            data.dWindkesselflowP_dPar = dqp
                - data.dWindkesselOutflow_dPar;

            data.dWindkesselflowP_dPd = -dqp;

            data.dWindkesselflowD_dPar = -dqp;

            data.dWindkesselflowD_dPd = dqp + dqd;
        }
      };
    };

    }
  /**
   * @brief Evaluates the Windkessel outflow and its pressure derivatives.
   *
   * The outflow is the positive part of the aortic flow rate
   * @f$ Q_{ar}^+ = \max(0, K_{ar}(p_{v,\text{mid}} - p_{ar,\text{mid}})) @f$.
   *
   * @tparam Input Model input parameter type.
   */
  template <class Input, typename RheologyModel = Rheology::Newtonian>
  class WindkesselOutflowEvaluator
  {
    public:
      /**
       * @brief Construct with model input parameters.
       * @param[in] input Model parameters (aortic valve conductance).
       */
      explicit WindkesselOutflowEvaluator(const Input& input)
        : m_input(input)
      {}

      /**
       * @brief Evaluate outflow and pressure derivatives.
       *
       * Reads midpoint pressures from @p data and writes
       * @p data.windkesselOutflow and its pressure derivatives.
       *
       * @tparam EvalData Evaluation data structure type.
       * @param[in,out] data Evaluation data with midpoint pressures set.
       */
      template <class EvalData>
      void evaluate(EvalData& data) const
      {
        using Scalar = decltype(data.y);
        RheologyModel::template computeFlow<Scalar>(data, m_input);
      }

    private:
      const Input& m_input;
  };
}

#endif
