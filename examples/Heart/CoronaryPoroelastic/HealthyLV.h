/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef EXAMPLES_HEART_CORONARYPOROELASTIC_HEALTHYLV_H
#define EXAMPLES_HEART_CORONARYPOROELASTIC_HEALTHYLV_H

/**
 * @file HealthyLV.h
 * @brief Healthy-adult calibration of the 0D poroelastic left ventricle.
 *
 * The resistive parameters are not tuned: they follow from physiological
 * targets (Table below) and two structural assumptions.
 *
 * 1. Windkessel branches (Murray). In a symmetric bifurcating tree obeying
 *    Murray's law, @f$ Q_k \propto r_k^3 @f$, the nominal wall shear rate
 *    @f$ 4Q_k/(\pi r_k^3) @f$ is the same in every generation. For any
 *    generalized-Newtonian fluid every generation then shares one apparent
 *    viscosity, and a tree of @f$ N @f$ generations with @f$ L_k = c\,r_k @f$
 *    is exactly equivalent, for all flows, to ONE tube of radius @f$ r_0 @f$
 *    and length @f$ L = N c\, r_0 @f$ (each generation carries the same
 *    pressure drop @f$ 2 c\,\tau_w @f$). The root radius follows from the
 *    Murray wall shear stress @f$ \tau^* @f$ at the mean branch flow,
 *    @f$ r_0^3 = 4\bar Q / (\pi \dot\gamma_{nom}(\tau^*)) @f$, and the length
 *    from the target resistance, @f$ L = r_0 R \bar Q / (2\tau^*) @f$.
 *    Proximal and distal branches share @f$ r_0 @f$; their lengths split as
 *    @f$ R_p : R_d @f$.
 *
 * 2. Perfusion conductances. The mean myocardial blood flow
 *    @f$ \bar Q_{cor} = \mathrm{MBF}\,\rho_t V_{w0} @f$ crosses the arterial
 *    and venous limbs in series; the venular pressure fraction @f$ f_v @f$
 *    fixes the split,
 *    @f$ \gamma_{ar} = \bar Q_{cor}/((1-f_v)\Delta P) @f$,
 *    @f$ \gamma_{ven} = \bar Q_{cor}/(f_v \Delta P) @f$,
 *    @f$ \Delta P = \bar p_{ar} - p_{sv} @f$. With @f$ \kappa = 2 @f$ the
 *    conductances scale as @f$ (\Phi/\phi_0)^2 @f$ (Poiseuille bed of fixed
 *    vessel count, @f$ R \propto V^{-2} @f$).
 *
 * Tuned by hand (12-cycle periodic regime): contractility @f$ \sigma_0 @f$,
 * passive exponents @f$ C_1 = C_3 @f$, aortic valve conductance and the two
 * compliances. Periodic regime (HR 70.6 bpm): EDV 118 mL, ESV 46 mL,
 * EF 61 %, CO 5.1 L/min, LVP 116 mmHg, EDP 7.6 mmHg, aortic 112/79 mmHg;
 * coronary inflow 116 mL/min (0.9 mL/min/g), systolic/diastolic mean inflow
 * 0.66, venous outflow systolic-dominant, @f$ p_f @f$ 7-62 mmHg,
 * @f$ \Delta\Phi \approx 0.008 @f$ (about 1 mL of intramyocardial blood
 * expelled per beat).
 *
 * References: Caruel et al., BMMB 2014 (doi:10.1007/s10237-013-0544-6);
 * Chapelle et al., Comput Mech 2010 (doi:10.1007/s00466-009-0452-x);
 * Barnafi et al., SIAM J Appl Math 2022 (doi:10.1137/21M1424482);
 * Murray, PNAS 1926; Sherman, J Gen Physiol 1981.
 */

#include <algorithm>
#include <cmath>
#include <numbers>
#include <type_traits>

#include "Rodin/Heart/PoroelasticSphere.h"

namespace Rodin::Examples::Heart::CoronaryPoroelastic
{
  using Model = Rodin::Heart::PoroelasticSphereT<>;

  /// Physiological targets of a healthy adult at rest (SI units).
  struct HealthyTargets
  {
      Real mmHg = 133.322;

      Real cardiacOutput = 5.0e-3 / 60.0;   ///< 5 L/min.
      Real meanArterialPressure = 96.0 * 133.322; ///< Design MAP (achieved ~ 94 mmHg).
      Real venousPressure = 1.0e3;         ///< p_sv (coronary sinus / RA).
      Real proximalResistanceFraction = 0.08; ///< R_p / (R_p + R_d).

      Real murrayWallShearStress = 1.5; ///< tau* [Pa] (arterial, Kamiya-Togawa).
      Real lengthToRadius = 20.0;       ///< c = L_k / r_k (diagnostic N only).

      Real myocardialBloodFlow = 0.8 / 60.0 * 1.0e-3; ///< 0.8 mL/min/g, in m^3/s/kg.
      Real tissueDensity = 1053.0;      ///< kg/m^3.
      Real venularPressureFraction = 0.2; ///< f_v: resting p_f = p_sv + f_v (MAP - p_sv) ~ 25 mmHg.
      Real perfusionConductanceExponent = 2.0; ///< kappa of gamma(Phi).
  };

  /// Equivalent tube of a Murray tree (radius, length, generations).
  struct MurrayBranch
  {
      Real radius = 0.0;
      Real length = 0.0;
      Real generations = 0.0;
      Real apparentViscosity = 0.0;
  };

  /// Nominal wall shear rate 4Q/(pi r^3) at the wall shear stress tauW, for
  /// the Windkessel rheology actually used by the 0D branches.
  inline Real nominalShearRate(const Model::Input& in, Real tauW)
  {
    namespace Rh = Rodin::Heart::CCMLC2014::Physics::Rheology;
    using WR = Rodin::Heart::CCMLC2014::Model::WindkesselRheology;
    // Unit tube: dp = 2 L tau_w / r with r = 1, L = 1.
    const Real dp = 2.0 * tauW;
    const Real fallback = 8.0 * 0.0035 / std::numbers::pi_v<Real>;
    Real q = dp / fallback;
    switch (in.windkesselRheology)
    {
      case WR::Cross:
        q = Rh::Cross::flowLaw(dp, 1.0, 1.0, fallback, in).first;
        break;
      case WR::CarreauYasuda:
        q = Rh::CarreauYasuda::flowLaw(dp, 1.0, 1.0, fallback, in).first;
        break;
      case WR::PowerLaw:
        q = Rh::PowerLaw::flowLaw(dp, 1.0, 1.0, fallback, in).first;
        break;
      case WR::Quemada:
        q = Rh::Quemada::flowLaw(dp, 1.0, 1.0, fallback, in).first;
        break;
      default:
        break;
    }
    return 4.0 * q / std::numbers::pi_v<Real>;
  }

  /// Murray equivalent tube carrying the mean flow Q with resistance R.
  inline MurrayBranch murrayBranch(
    const Model::Input& in, const HealthyTargets& tg, Real Q, Real R)
  {
    MurrayBranch b;
    const Real tau = tg.murrayWallShearStress;
    const Real g = nominalShearRate(in, tau);
    b.apparentViscosity = tau / g;
    b.radius = std::cbrt(4.0 * Q / (std::numbers::pi_v<Real> * g));
    b.length = b.radius * R * Q / (2.0 * tau);
    b.generations = b.length / (tg.lengthToRadius * b.radius);
    return b;
  }

  /// Diagnostics of the resistive calibration.
  struct HealthyCalibration
  {
      Real wallVolume = 0.0;
      Real coronaryFlow = 0.0;
      Real windkesselFlow = 0.0;
      Real Rp = 0.0;
      Real Rd = 0.0;
      MurrayBranch proximal;
      MurrayBranch distal;
  };

  inline Real smoothstep(Real s)
  {
    s = std::clamp<Real>(s, 0.0, 1.0);
    return s * s * (3.0 - 2.0 * s);
  }

  /// Active drive u(t) (Caruel et al. 2014 shape), period 0.85 s, smoothed edges.
  inline Real activation(Real t, Real period = 0.85)
  {
    const Real tau = t - period * std::floor(t / period);
    const Real tRampStart = 0.15, tRampEnd = 0.21, tPlateauEnd = 0.36;
    const Real tRelaxEnd = 0.45, tNegativeEnd = 0.6;
    const Real positiveValue = 35.0, negativeValue = -20.0;

    if (tau < tRampStart)
      return 0.0;
    if (tau < tRampEnd)
      return positiveValue * smoothstep((tau - tRampStart) / (tRampEnd - tRampStart));
    if (tau < tPlateauEnd)
      return positiveValue;
    if (tau < tRelaxEnd)
      return positiveValue +
        (negativeValue - positiveValue) *
        smoothstep((tau - tPlateauEnd) / (tRelaxEnd - tPlateauEnd));
    if (tau < tNegativeEnd)
      return negativeValue;
    return negativeValue *
      (1.0 - smoothstep((tau - tNegativeEnd) / (period - tNegativeEnd)));
  }

  inline Real relaxationTarget(Real ec)
  {
    const Real lowEc = 0.0, highEc = 2.0, lowValue = 1.6, highValue = 1.0;
    if (ec <= lowEc)
      return lowValue;
    if (ec >= highEc)
      return highValue;
    const Real s = (ec - lowEc) / (highEc - lowEc);
    return (1.0 - s) * lowValue + s * highValue;
  }

  inline Real relaxationTargetDerivative(Real ec)
  {
    const Real lowEc = 0.0, highEc = 2.0, lowValue = 1.6, highValue = 1.0;
    if (ec <= lowEc || ec >= highEc)
      return 0.0;
    return (highValue - lowValue) / (highEc - lowEc);
  }

  /// Left-atrial pressure (5 -> 7.5 -> 9.4 mmHg), period 0.85 s.
  inline Real atrialPressure(Real t, Real period = 0.85)
  {
    const Real tau = t - period * std::floor(t / period);
    const Real minValue = 500.0, maxValue = 1000.0, secondThreshold = 1250.0;
    const Real t1 = 0.02, t2 = 0.15, t3 = 0.17, t4 = 0.56, t5 = 0.62, t6 = 0.85;
    auto ramp = [](Real a, Real b, Real s) { return a + (b - a) * smoothstep(s); };

    if (tau < t1)
      return ramp(minValue, maxValue, tau / t1);
    if (tau < t2)
      return maxValue;
    if (tau < t3)
      return ramp(maxValue, minValue, (tau - t2) / (t3 - t2));
    if (tau < t4)
      return ramp(minValue, secondThreshold, (tau - t3) / (t4 - t3));
    if (tau < t5)
      return secondThreshold;
    if (tau < t6)
      return ramp(secondThreshold, minValue, (tau - t5) / (t6 - t5));
    return minValue;
  }

  /**
   * @brief Healthy 0D poroelastic LV input.
   *
   * @param[in] tg Physiological targets.
   * @param[out] cal Optional resistive-calibration diagnostics.
   */
  inline Model::Input makeHealthyInput(
    const HealthyTargets& tg = HealthyTargets(), HealthyCalibration* cal = nullptr)
  {
    Model::Input in;

    // Geometry: unloaded cavity 58 mL, wall 1.1 cm (LV mass ~ 128 g).
    in.R0 = 2.4e-2;
    in.d0 = 1.1e-2;
    in.phi0 = 0.1;

    // Active and viscous laws (CCMLC2014).
    in.Es = 3.0e6;
    in.mu = 70.0;
    in.eta = 70.0;
    in.alpha = 1.5;
    in.alphaR = 0.12;
    in.k0 = 1.0e5;
    in.sigma0 = 2.0e5;

    // Intramyocardial blood storage: C_im = V_w0 / K_Phi.
    in.KPhi = 2.0e5;

    // Arterial compliances (total 1.2e-8 m^3/Pa = 1.6 mL/mmHg).
    in.Cp = 5.0e-9;
    in.Cd = 7.0e-9;

    // Blood rheology of the Windkessel branches (Cross).
    in.mu_0 = 5.35;
    in.mu_Inf = 0.0033;
    in.lambda = 14.445;
    in.n = 0.8;
    in.m = 0.003;
    in.yasuda = 0.62;
    in.mu_plasma = 0.0032704;
    in.k_0 = 3.5678;
    in.gamma_c = 10.2754;
    in.k_Inf = 1.5352;
    in.windkesselRheology = Rodin::Heart::CCMLC2014::Model::WindkesselRheology::Cross;

    in.Kat = 6.0e-7;
    in.Kp = 5.0e-11;
    in.Kar = 5.0e-7; // peak systolic LV-Ao gradient ~ 5 mmHg
    in.cavityCapacity = 5.0e-12;

    in.absRegularization = 1e-14;
    in.wallQuadraturePoints = 8;

    in.initFibDef = 0.0;
    in.initActiveStiffness = 0.0;
    in.initActiveStress = 0.0;

    const Real pSv = tg.venousPressure;
    in.pSv = [pSv](Real) { return pSv; };
    in.pAt = [](Real t) { return atrialPressure(t); };
    in.u = [](Real t) { return activation(t); };
    in.m0 = relaxationTarget;
    in.dm0 = relaxationTargetDerivative;

    using PassiveEnergy = std::decay_t<decltype(in.passiveEnergy)>;
    typename PassiveEnergy::Parameters hp;
    hp.mu1 = 0.0;
    hp.mu2 = 0.0;
    hp.C0 = 1.9e3;
    hp.C1 = 0.5;
    hp.C2 = 1.9e3;
    hp.C3 = 0.5;
    in.passiveEnergy = PassiveEnergy(hp);

    // ---- Resistive calibration -----------------------------------------
    const Real Rin = in.R0;
    const Real Rout = in.R0 + in.d0;
    const Real Vw0 =
      4.0 / 3.0 * std::numbers::pi_v<Real> * (Rout * Rout * Rout - Rin * Rin * Rin);
    const Real dP = tg.meanArterialPressure - pSv;

    const Real Qcor = tg.myocardialBloodFlow * tg.tissueDensity * Vw0;
    const Real fv = tg.venularPressureFraction;
    in.gammaAr = Qcor / ((1.0 - fv) * dP);
    in.gammaVen = Qcor / (fv * dP);
    in.gammaArExponent = tg.perfusionConductanceExponent;
    in.gammaVenExponent = tg.perfusionConductanceExponent;

    const Real Qwk = tg.cardiacOutput - Qcor;
    const Real Rtot = dP / Qwk;
    in.Rp = tg.proximalResistanceFraction * Rtot;
    in.Rd = (1.0 - tg.proximalResistanceFraction) * Rtot;

    const MurrayBranch prox = murrayBranch(in, tg, Qwk, in.Rp);
    const MurrayBranch dist = murrayBranch(in, tg, Qwk, in.Rd);
    in.proximalRadius = prox.radius;
    in.proximalLength = prox.length;
    in.distalRadius = dist.radius;
    in.distalLength = dist.length;

    if (cal)
    {
      cal->wallVolume = Vw0;
      cal->coronaryFlow = Qcor;
      cal->windkesselFlow = Qwk;
      cal->Rp = in.Rp;
      cal->Rd = in.Rd;
      cal->proximal = prox;
      cal->distal = dist;
    }
    return in;
  }

  /// Initial state close to end-diastole of the periodic regime.
  inline Model::State makeHealthyInitialState(const Model::Input& in)
  {
    Model::State s0;
    s0.t = 0.0;
    s0.y = 0.0;
    s0.phi = in.phi0;
    s0.pv = in.pAt(0.0) - 100.0;
    s0.par = 11000.0;
    s0.pd = 10500.0;
    s0.ec = in.initFibDef;
    s0.gamma = std::sqrt(std::max<Real>(in.initActiveStiffness, 0.0));
    s0.beta = (s0.gamma > 0.0) ? (in.initActiveStress / s0.gamma) : 0.0;
    s0.kc = s0.gamma * s0.gamma;
    s0.tauc = s0.gamma * s0.beta;
    s0.w = in.m0(s0.ec);
    return s0;
  }
}

#endif
