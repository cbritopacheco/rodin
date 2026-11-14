/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file PoroelasticLV.h
 * @brief Shared parameters, forcings and study scenarios of the reduced
 * poroelastic left ventricle (PoroelasticSphere_0D and
 * CoronaryArtery_Explicit_PoroElastic).
 *
 * One blood rheology (Carreau-Yasuda, Table 1 of the manuscript) feeds the
 * resolved coronary flow, the outlet compartments, the pore viscosity of the
 * venous drainage and the systemic Windkessel; a scenario scales it, or the
 * cardiac parameters, consistently at every scale.
 */
#ifndef RODIN_EXAMPLES_HEART_POROELASTICLV_H
#define RODIN_EXAMPLES_HEART_POROELASTICLV_H

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "Rodin/Heart/PoroelasticSphere.h"

namespace Rodin::Examples::Heart
{
  using PoroelasticLV = Rodin::Heart::PoroelasticSphereT<>;

  /// Carreau-Yasuda blood viscosity, mu = muInf + (mu0 - muInf)(1 + (lambda g)^a)^((n-1)/a).
  struct BloodRheology
  {
      Real mu0 = 0.056;
      Real muInf = 3.45e-3;
      Real lambda = 1.902;
      Real n = 0.22;
      Real yasuda = 1.25;
      Real gammaRegularization = 1.0e-3;

      Real operator()(Real shearRate) const
      {
        return muInf +
          (mu0 - muInf) *
          std::pow(1.0 + std::pow(lambda * shearRate, yasuda), (n - 1.0) / yasuda);
      }
  };

  /**
   * @brief Study scenario. Factors multiply the healthy reference values.
   *
   * - healthy:          reference.
   * - aortic_stenosis:  valvular flow restriction, aortic conductance x0.2.
   * - hyperviscosity:   blood viscosity x1.6 (both plateaus) in the 3D flow,
   *                     outlet compartments, pore drainage and Windkessel.
   * - lvh:              concentric hypertrophy: wall thickness x1.35, passive
   *                     stiffness x2.5, vascular compliance x0.6, permeability x0.5.
   */
  struct Scenario
  {
      std::string name = "healthy";
      Real aorticConductance = 1.0;
      Real bloodViscosity = 1.0;
      Real wallThickness = 1.0;
      Real passiveStiffness = 1.0;
      Real storageCompliance = 1.0;
      Real permeability = 1.0;

      static Scenario make(const std::string& name)
      {
        Scenario s;
        s.name = name;
        if (name == "healthy")
          return s;
        if (name == "aortic_stenosis")
        {
          s.aorticConductance = 0.2;
          return s;
        }
        if (name == "hyperviscosity")
        {
          s.bloodViscosity = 1.6;
          return s;
        }
        if (name == "lvh")
        {
          s.wallThickness = 1.35;
          s.passiveStiffness = 2.5;
          s.storageCompliance = 0.6;
          s.permeability = 0.5;
          return s;
        }
        throw std::runtime_error("Unknown scenario '" + name +
          "' (healthy, aortic_stenosis, hyperviscosity, lvh).");
      }

      BloodRheology blood() const
      {
        BloodRheology b;
        b.mu0 *= bloodViscosity;
        b.muInf *= bloodViscosity;
        return b;
      }
  };

  inline Real smoothstep(Real s)
  {
    s = s < 0.0 ? 0.0 : (s > 1.0 ? 1.0 : s);
    return s * s * (3.0 - 2.0 * s);
  }

  inline Real periodicActivation(Real t)
  {
    const Real T = 0.85;
    const Real tau = t - T * std::floor(t / T);
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
    return negativeValue * (1.0 - smoothstep((tau - tNegativeEnd) / (T - tNegativeEnd)));
  }

  inline Real loadDependentRelaxationM0(Real ec)
  {
    const Real lowEc = 0.0, highEc = 2.0, lowValue = 1.6, highValue = 1.0;
    if (ec <= lowEc)
      return lowValue;
    if (ec >= highEc)
      return highValue;
    const Real s = (ec - lowEc) / (highEc - lowEc);
    return (1.0 - s) * lowValue + s * highValue;
  }

  inline Real loadDependentRelaxationDM0(Real ec)
  {
    const Real lowEc = 0.0, highEc = 2.0, lowValue = 1.6, highValue = 1.0;
    if (ec <= lowEc || ec >= highEc)
      return 0.0;
    return (highValue - lowValue) / (highEc - lowEc);
  }

  inline Real atrialPressure(Real t)
  {
    const Real T = 0.85;
    const Real tau = t - T * std::floor(t / T);
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
   * @brief Reduced poroelastic LV input for a scenario.
   *
   * @param scenario Study scenario.
   * @param gammaAr Lumped myocardial supply conductance (0 when the supply
   *   comes entirely from a resolved coronary tree).
   */
  inline PoroelasticLV::Input makePoroelasticLVInput(const Scenario& scenario, Real gammaAr)
  {
    PoroelasticLV::Input in;

    // Shell: mid-surface radius and thickness.
    in.R0 = 2.4e-2;
    in.d0 = 1.45e-2 * scenario.wallThickness;
    in.phi0 = 0.1;

    // Active (CCMLC2014) and skeletal viscosity.
    in.Es = 3.0e6;
    in.mu = 70.0;
    in.eta = 70.0;
    in.alpha = 1.5;
    in.alphaR = 0.12;
    in.k0 = 1.0e5;
    in.sigma0 = 1.25e5;

    // Passive Usyk-type energy (fitted to the former reduced-invariant law).
    {
      using PassiveEnergy = std::decay_t<decltype(in.passiveEnergy)>;
      typename PassiveEnergy::Parameters pp;
      pp.c0 = 100.0 * scenario.passiveStiffness;
      pp.brr = 1.24;
      pp.btt = 8.4;
      pp.brt = 0.0;
      in.passiveEnergy = PassiveEnergy(pp);
    }

    // Pore storage dPhi/dp_P = C0 exp(-p_P/p0), Delta Phi_max = C0 p0.
    in.storageCompliance = 5.0e-6 * scenario.storageCompliance;
    in.storagePressureScale = 1.0e4;
    in.storageSaturationGuard = 0.95;

    // Hydraulics: K, alpha_pm, G_ven,0 = 7e-10 m^3/(Pa s) x nu(300 1/s).
    in.gammaAr = gammaAr;
    in.permeability = 1.0e-11 * scenario.permeability;
    in.poreShearFactor = std::sqrt(2.0);
    in.venousConductance = 2.67e-12;
    {
      const BloodRheology b = scenario.blood();
      in.poreViscosity.eta0 = b.mu0;
      in.poreViscosity.etaInf = b.muInf;
      in.poreViscosity.lambda = b.lambda;
      in.poreViscosity.n = b.n;
      in.poreViscosity.a = b.yasuda;
      in.poreViscosity.shearRegularization = b.gammaRegularization;
    }

    // Systemic Windkessel (Cross rheology of the proximal and distal tubes).
    in.Rp = 5.0e7;
    in.Cp = 6e-9;
    in.Rd = 1.0e8;
    in.Cd = 1.0e-9;
    in.mu_0 = 5.35 * scenario.bloodViscosity;
    in.mu_Inf = 0.0033 * scenario.bloodViscosity;
    in.lambda = 14.445;
    in.n = 0.8;
    in.m = 0.003;
    in.yasuda = 0.62;
    in.mu_plasma = 0.0032704;
    in.k_0 = 3.5678;
    in.gamma_c = 10.2754;
    in.k_Inf = 1.5352;
    in.proximalRadius = 0.0125;
    in.proximalLength = 0.35;
    in.distalRadius = 0.002;
    in.distalLength = 0.55;
    in.windkesselRheology = Rodin::Heart::CCMLC2014::Model::WindkesselRheology::Cross;

    // Valves.
    in.Kat = 6.0e-7;
    in.Kp = 5.0e-11;
    in.Kar = 1.0e-7 * scenario.aorticConductance;
    in.cavityCapacity = 5.0e-12;

    in.absRegularization = 1e-14;
    in.initFibDef = 0.0;
    in.initActiveStiffness = 0.0;
    in.initActiveStress = 0.0;

    in.pSv = [](Real) { return 1.0e3; };
    in.pAt = atrialPressure;
    in.u = periodicActivation;
    in.m0 = loadDependentRelaxationM0;
    in.dm0 = loadDependentRelaxationDM0;
    return in;
  }

  /// Initial state at end diastole, rest pore storage.
  inline PoroelasticLV::State makePoroelasticLVInitialState(const PoroelasticLV::Input& in)
  {
    PoroelasticLV::State s0;
    s0.t = 0.0;
    s0.y = 0.0;
    s0.phi = in.phi0;
    s0.pv = in.pAt(0.0) - 100.0;
    s0.par = 11000.0;
    s0.pd = 10000.0;
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
