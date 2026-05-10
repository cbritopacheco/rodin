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

#include <cmath>
#include <iostream>
#include <numbers>
#include <utility>

#include "Rodin/Math/RootFinding/NewtonRaphson.h"
#include "Rodin/Math/RungeKutta/RK4.h"

namespace Rodin::Heart::CCMLC2014::Physics
{
  struct SolverConfig
  {
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
  };

  namespace Rheology
  {
    struct Newtonian
    {
      template <typename Input>
      static std::pair<Real, Real> flowLaw(
          Real dp, Real L, Real radius, Real fallbackResistance, const Input &input)
      {
        return { dp / fallbackResistance, 1.0 / fallbackResistance };
      }
  };

  struct PowerLaw
  {
    template <typename Input>
    static std::pair<Real, Real> flowLaw(
        Real dp, Real L, Real radius, Real fallbackResistance, const Input &input)
    {
      if (L <= 0.0 || radius <= 0.0)
        return {dp / fallbackResistance, 1.0 / fallbackResistance};

      const Real adp = std::abs(dp);
      const Real sgn = (dp >= 0.0) ? 1.0 : -1.0;

      if (adp < SolverConfig::pressureDropTolerance)
        return {dp / fallbackResistance, 1.0 / fallbackResistance};

      const Real m = input.m;
      const Real n = input.n;
      const Real invN = 1.0 / n;
      const Real tauW = (radius * adp) / (2.0 * L);

      const Real piR3 = std::numbers::pi_v<Real> * std::pow(radius, 3.0);
      const Real qAbs = (n * piR3 / (3.0 * n + 1.0)) * std::pow(tauW / m, invN);

      const Real dqAbs = invN * (qAbs / adp);

      if (!std::isfinite(qAbs) || dqAbs <= 0.0)
        return {dp / fallbackResistance, 1.0 / fallbackResistance};

      return {sgn * qAbs, dqAbs};
    }
  };

  struct Quemada
  {
    template <typename Input>
    static std::pair<Real, Real> flowLaw(
        Real dp, Real L, Real radius, Real fallbackResistance, const Input &input)
    {
      if (L <= 0.0 || radius <= 0.0)
        return {dp / fallbackResistance, 1.0 / fallbackResistance};

      const Real adp = std::abs(dp);
      const Real sgn = (dp >= 0.0) ? 1.0 : -1.0;

      if (adp < SolverConfig::pressureDropTolerance)
        return {dp / fallbackResistance, 1.0 / fallbackResistance};

      const Real tauW = (adp * radius) / (2.0 * L);
      const Real k0 = input.k_0;
      const Real kInf = input.k_Inf;
      const Real phi = input.phi_quemada;
      const Real gammaC = input.gamma_c;
      const Real muPl = input.mu_plasma;
      const Real muInfBlood = input.mu_Inf;

      const Real oneMinusHalfKInfPhi = 1.0 - 0.5 * kInf * phi;
      const Real tau_0 = muInfBlood * gammaC *
                         std::pow(0.5 * phi * (k0 - kInf), 2.0) /
                         std::pow(oneMinusHalfKInfPhi, 4.0);
      const Real nu_inf = muInfBlood / std::pow(oneMinusHalfKInfPhi, 2.0);
      const Real lambda = gammaC * std::pow((1.0 - 0.5 * k0 * phi) / oneMinusHalfKInfPhi, 2.0);

      const Real sqrtTau0 = std::sqrt(tau_0);
      const Real sqrtNuInfLambda = std::sqrt(nu_inf * lambda);

      const Real q = (sqrtTau0 - sqrtNuInfLambda) / std::max(sqrtTau0 + sqrtNuInfLambda, 1e-12);
      const Real alpha = (sqrtTau0 + sqrtNuInfLambda) / std::max(std::sqrt(tauW), 1e-12);

      Real P[8];
      {
        const Real q2 = q * q;
        const Real q3 = q2 * q;
        const Real q4 = q3 * q;
        const Real q5 = q4 * q;
        const Real q6 = q5 * q;
        const Real q7 = q6 * q;

        P[0] = -(1.0 / 7.0) * (q + 8.0);
        P[1] = -(1.0 / 42.0) * (13.0 * q2 - 8.0 * q - 7.0);
        P[2] = -(1.0 / 210.0) * (143.0 * q3 - 88.0 * q2 - 113.0 * q + 48.0);
        P[3] = -(1.0 / 840.0) * (1287.0 * q4 - 792.0 * q3 - 1342.0 * q2 + 632.0 * q + 175.0);
        P[4] = -(1.0 / 840.0) * (3003.0 * q5 - 1848.0 * q4 - 3894.0 * q3 + 1944.0 * q2 + 1011.0 * q - 256.0);
        P[5] = -(1.0 / 1680.0) * (15015.0 * q6 - 9240.0 * q5 - 23331.0 * q4 + 12096.0 * q3 + 9081.0 * q2 - 3179.0 * q - 525.0);
        P[6] = -(1.0 / 1680.0) * (45045.0 * q7 - 27720.0 * q6 - 82005.0 * q5 + 43680.0 * q4 + 42819.0 * q3 - 17304.0 * q2 - 5619.0 * q + 1024.0);
        P[7] = -(1.0 / 16.0) * std::pow(1.0 - q, 2.0) * (1.0 + q) * (429.0 * q5 + 165.0 * q4 - 330.0 * q3 - 90.0 * q2 + 45.0 * q + 5.0);
      }

      Real F = 0.0, dFdAlpha = 0.0;

      if (alpha > 1.0)
      {
        const Real common = (2.0 / 9.0) * std::pow(1.0 - q, 2.0) * (1.0 + q);
        F = 0.5 * (0.25 * std::pow(1.0 - q, 2.0) + common / alpha);
        dFdAlpha = 0.5 * (-common / (alpha * alpha));
      }
      else
      {
        const Real a = std::max(alpha, 1e-10);
        const Real a2 = a * a;
        const Real a7 = a2 * a2 * a2 * a;
        const Real a8 = a7 * a;

        Real S = 0.0, derS = 0.0;
        Real aPow = a;
        for (int i = 0; i < 7; ++i)
        {
          S += aPow * P[i];
          derS += (i + 1) * (aPow / a) * P[i];
          aPow *= a;
        }

        const Real T = std::max(0.0, 1.0 - 2.0 * a * q + a2);
        const Real sqrtT = std::sqrt(T);

        const Real safeDenom = std::max(1.0 - q, 1e-10);
        const Real argLog = (1.0 - a * q + sqrtT) / (a * safeDenom);
        const Real logT = (argLog <= 1e-12) ? 0.0 : std::log(argLog);

        const Real bracket = 1.0 - (8.0 / 7.0) * a * (1.0 + q) + (4.0 / 3.0) * a2
                             - a8 * P[6] + (1.0 + S) * sqrtT + a8 * P[7] * logT;
        F = 0.5 * bracket;

        const Real eps = 1e-12;
        const Real dSqrtT = (a - q) / std::max(sqrtT, eps);
        const Real dLogT = ((dSqrtT - q) / std::max(1.0 - a * q + sqrtT, eps)) - (1.0 / a);

        Real dBracket = -(8.0 / 7.0) * (1.0 + q) + (8.0 / 3.0) * a - 8.0 * a7 * P[6];
        dBracket += (derS * sqrtT + (1.0 + S) * dSqrtT);
        dBracket += P[7] * (8.0 * a7 * logT + a8 * dLogT);

        dFdAlpha = 0.5 * dBracket;
      }

      const Real constPart = std::pow(oneMinusHalfKInfPhi, 2.0) / (4.0 * muPl);
      const Real IQ = std::pow(tauW, 4.0) * constPart * F;

      const Real geomFactor = (8.0 * std::numbers::pi_v<Real> * std::pow(L, 3.0)) / std::pow(adp, 3.0);
      const Real qAbs = geomFactor * IQ;


      const Real dqAbs = (qAbs / adp) * (1.0 - 0.5 * alpha * dFdAlpha / std::max(F, 1e-12));

      if (!std::isfinite(qAbs) || !std::isfinite(dqAbs) || dqAbs <= 0.0)
        return {dp / fallbackResistance, 1.0 / fallbackResistance};

      return {sgn * qAbs, dqAbs};
    }
  };

  struct Cross
  {
    template <typename Input>
    static std::pair<Real, Real> flowLaw(
        Real dp, Real L, Real radius, Real fallbackResistance, const Input &input)
    {
      if (L <= 0.0 || radius <= 0.0)
        return {dp / fallbackResistance, 1.0 / fallbackResistance};

      const Real mu0 = input.mu_0;
      const Real muInf = input.mu_Inf;
      const Real lambda = input.lambda;

      // In the Cross law, input.n is the Cross transition exponent.
      const Real n = input.n;
      const Real delta = mu0 - muInf;

      const Real sgn = (dp >= 0.0) ? 1.0 : -1.0;
      const Real adp = std::abs(dp);

      const Real R0 =
          8.0 * mu0 * L / (std::numbers::pi_v<Real> * std::pow(radius, 4.0));

      if (adp < SolverConfig::pressureDropTolerance)
        return {dp / R0, 1.0 / R0};

      const Real tauW = radius * adp / (2.0 * L);

      auto mu = [&](Real g) -> Real
      {
        const Real x = std::pow(lambda * g, n);
        return muInf + delta / (1.0 + x);
      };

      auto dmu = [&](Real g) -> Real
      {
        if (g <= 0.0)
          return 0.0;

        const Real x = std::pow(lambda * g, n);
        return -delta * n * std::pow(lambda, n) * std::pow(g, n - 1.0) /
               std::pow(1.0 + x, 2.0);
      };

      auto tauMinusTauW = [&](Real g) -> std::pair<Real, Real>
      {
        const Real visc = mu(g);
        const Real dvisc = dmu(g);
        return {g * visc - tauW, visc + g * dvisc};
      };

      Math::RootFinding::NewtonRaphson<Real> rootFinder(
          SolverConfig::shearAbsoluteTolerance,
          SolverConfig::shearRelativeTolerance, SolverConfig::shearStepTolerance,
          SolverConfig::shearMaxIterations);

      Real gHi = std::max<Real>(tauW / muInf, SolverConfig::minShearRate);

      for (int k = 0; k < SolverConfig::maxBracketIterations && tauMinusTauW(gHi).first < 0.0; ++k)
        gHi *= 2.0;

      if (tauMinusTauW(gHi).first < 0.0)
      {
        std::cerr << "Warning: failed to bracket wall shear rate. "
                  << "Using Poiseuille fallback.\n";
        return {dp / R0, 1.0 / R0};
      }

      const auto gammaRoot = rootFinder.solve(
          tauMinusTauW, 0.5 * gHi, SolverConfig::shearStepTolerance, gHi);

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

        const Real visc = mu(g);
        const Real dvisc = dmu(g);
        const Real dtau = visc + g * dvisc;

        return std::pow(g, 3.0) * visc * visc * dtau;
      };

      Math::RungeKutta::RK4 integrator;

      const int steps = SolverConfig::integralSteps;
      const Real h = gammaW / static_cast<Real>(steps);

      Real I = 0.0;

      auto rhs = [&](Real g, Real y) -> Real
      {
        (void)y;
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

      const Real qAbs = std::numbers::pi_v<Real> * std::pow(radius, 3.0) * I /
                        std::pow(tauW, 3.0);

      const Real dqAbs =
          (std::numbers::pi_v<Real> * std::pow(radius, 3.0) * gammaW -
           3.0 * qAbs) /
          adp;

      if (!std::isfinite(qAbs) || !std::isfinite(dqAbs) || dqAbs <= 0.0)
      {
        std::cerr << "Warning: invalid WRMS flow derivative. "
                  << "Using Poiseuille fallback.\n";
        return {dp / R0, 1.0 / R0};
      }

      return {sgn * qAbs, dqAbs};
    }
  };

  struct CarreauYasuda
  {
    template <typename Input>
    static std::pair<Real, Real> flowLaw(Real dp, Real L, Real radius,
                                         Real fallbackResistance,
                                         const Input &input)
    {
      if (L <= 0.0 || radius <= 0.0)
        return {dp / fallbackResistance, 1.0 / fallbackResistance};

      const Real mu0 = input.mu_0;
      const Real muInf = input.mu_Inf;
      const Real lambda = input.lambda;
      const Real n = input.n;
      const Real yasuda = input.yasuda;
      const Real delta = mu0 - muInf;

      const Real sgn = (dp >= 0.0) ? 1.0 : -1.0;
      const Real adp = std::abs(dp);

      const Real R0 =
          8.0 * mu0 * L / (std::numbers::pi_v<Real> * std::pow(radius, 4.0));

      if (adp < SolverConfig::pressureDropTolerance)
        return {dp / R0, 1.0 / R0};

      const Real tauW = radius * adp / (2.0 * L);

      auto mu = [&](Real g) -> Real
      {
        return muInf + delta * std::pow(1.0 + std::pow(lambda * g, yasuda),
                                        (n - 1.0) / yasuda);
      };

      auto dmu = [&](Real g) -> Real
      {
        const Real base = 1.0 + std::pow(lambda * g, yasuda);

        return delta * (n - 1.0) * std::pow(base, (n - 1.0 - yasuda) / yasuda) *
               std::pow(lambda, yasuda) * std::pow(g, yasuda - 1.0);
      };

      auto tauMinusTauW = [&](Real g) -> std::pair<Real, Real>
      {
        const Real m = mu(g);
        const Real dm = dmu(g);
        return {g * m - tauW, m + g * dm};
      };

      Math::RootFinding::NewtonRaphson<Real> rootFinder(
          SolverConfig::shearAbsoluteTolerance,
          SolverConfig::shearRelativeTolerance, SolverConfig::shearStepTolerance,
          SolverConfig::shearMaxIterations);

      Real gHi = std::max<Real>(tauW / muInf, SolverConfig::minShearRate);

      for (int k = 0; k < SolverConfig::maxBracketIterations &&
                      tauMinusTauW(gHi).first < 0.0;
           ++k)
        gHi *= 2.0;

      if (tauMinusTauW(gHi).first < 0.0)
      {
        std::cerr << "Warning: failed to bracket wall shear rate. "
                  << "Using Poiseuille fallback.\n";
        return {dp / R0, 1.0 / R0};
      }

      const auto gammaRoot = rootFinder.solve(
          tauMinusTauW, 0.5 * gHi, SolverConfig::shearStepTolerance, gHi);

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

        const Real m = mu(g);
        const Real dm = dmu(g);
        const Real dtau = m + g * dm;

        return std::pow(g, 3.0) * m * m * dtau;
      };

      Math::RungeKutta::RK4 integrator;

      const int steps = SolverConfig::integralSteps;
      const Real h = gammaW / static_cast<Real>(steps);

      Real I = 0.0;

      auto rhs = [&](Real g, Real y) -> Real
      {
        (void)y;
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

      const Real qAbs = std::numbers::pi_v<Real> * std::pow(radius, 3.0) * I /
                        std::pow(tauW, 3.0);

      const Real dqAbs =
          (std::numbers::pi_v<Real> * std::pow(radius, 3.0) * gammaW -
           3.0 * qAbs) /
          adp;

      if (!std::isfinite(qAbs) || !std::isfinite(dqAbs) || dqAbs <= 0.0)
      {
        std::cerr << "Warning: invalid WRMS flow derivative. "
                  << "Using Poiseuille fallback.\n";
        return {dp / R0, 1.0 / R0};
      }

      return {sgn * qAbs, dqAbs};
    }
  };
} // namespace Rheology

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
      explicit WindkesselOutflowEvaluator(const Input &input) : m_input(input)
      {}

      template <class EvalData> void evaluate(EvalData &data) const
      {
        using Scalar = decltype(data.y);

        const Scalar aorticFlowRate = m_input.Kar * (data.pvMid - data.parMid);
        data.windkesselOutflow =
            (aorticFlowRate > Scalar(0)) ? aorticFlowRate : Scalar(0);

        const Scalar valveDeriv =
            (aorticFlowRate > Scalar(0)) ? m_input.Kar : Scalar(0);
        data.dWindkesselOutflow_dPv = valveDeriv;
        data.dWindkesselOutflow_dPar = -valveDeriv;

        const auto [qp, dqp] =
            RheologyModel::flowLaw(data.parMid - data.pdMid, m_input.proximalLength,
                                   m_input.proximalRadius, m_input.Rp, m_input);

        const auto [qd, dqd] =
            RheologyModel::flowLaw(data.pSvMid - data.pdMid, m_input.distalLength,
                                   m_input.distalRadius, m_input.Rd, m_input);

        data.windkesselflowP = qp - data.windkesselOutflow;
        data.windkesselflowD = -qp - qd;

        data.dWindkesselflowP_dPar = dqp - data.dWindkesselOutflow_dPar;
        data.dWindkesselflowP_dPd = -dqp;
        data.dWindkesselflowD_dPar = -dqp;
        data.dWindkesselflowD_dPd = dqp + dqd;
      }

    private:
      const Input &m_input;
  };
} // namespace Rodin::Heart::CCMLC2014::Physics

#endif
