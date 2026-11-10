/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Coagulation.h
 * @brief Constants of the P-T-A-G reduced coagulation model and its Brinkman drag.
 *
 * Species (mol/m^3, per plasma volume) obey
 * @f$ \mathcal L c=\partial_t c|_\chi+(u-w)\cdot\nabla c-\nabla\cdot(D\nabla c) @f$
 * on the moving fluid domain, with mesh velocity @f$ w @f$:
 *
 *   L P = -R           prothrombin, the precursor
 *   L T =  R - R_I     thrombin, the enzyme
 *   L A =    - R_I     antithrombin, the stoichiometric inhibitor
 *   L G =      - R_F   fibrinogen, the substrate
 *
 *   R   = k_a T^2/(K_A + T) P/P_ref     lumped feedback amplification
 *   R_I = k_on A T                       1:1 irreversible inhibition
 *   R_F = k_F T G/(K_G + G)              Michaelis-Menten conversion
 *
 * Invariants of the continuous problem:
 *   P + T - A = P_0 - A_0,   P + T <= P_0,   A <= A_0,   0 <= G <= G_0.
 * Diagnostics: F = G_0 - G (pre-gel fibrin surrogate), I = A_0 - A.
 */
#ifndef RODIN_ATRIUM_FLUID_COAGULATION_H
#define RODIN_ATRIUM_FLUID_COAGULATION_H

#include <algorithm>
#include <cmath>
#include <Rodin/Types.h>

namespace Rodin::Examples::Heart::Atrium
{
  /**
   * @brief Constants of the P-T-A-G kinetics (SI: mol, m, s).
   *
   * 1 uM = 1e-3 mol/m^3. k_on = 7.1e3 M^-1 s^-1 = 7.1 m^3 mol^-1 s^-1, hence
   * k_on A_0 = 0.020 s^-1. k_F and K_G follow Higgins et al. (1983). K_A is
   * fixed by the frozen-pool activation threshold
   * T_th = k_on A_0 K_A/(k_a P_0/P_ref - k_on A_0) = 2.5 nM. k_a and j_0 are
   * the two calibration constants.
   */
  struct ConservingKinetics
  {
      Real diffusivityProtein = 4.6e-11; ///< D_P = D_T = D_A (m^2/s)
      Real diffusivityFibrinogen = 2.0e-11; ///< D_G (m^2/s)

      Real prothrombin0 = 1.4e-3; ///< P_0, inflow and initial (mol/m^3)
      Real antithrombin0 = 2.8e-3; ///< A_0 (mol/m^3)
      Real fibrinogen0 = 7.35e-3; ///< G_0: 2.5 g/L; 4 g/L = 11.76e-3
      Real prothrombinRef = 1.4e-3; ///< P_ref, fixed across cases
      Real inletThrombin = 0.0; ///< T at the veins (mol/m^3)

      Real amplificationRate = 0.10; ///< k_a (1/s), calibration constant
      Real amplificationHalf = 1.0e-5; ///< K_A (mol/m^3), 10 nM
      Real inhibitionRate = 7.1; ///< k_on (m^3 mol^-1 s^-1)
      Real conversionRate = 42.0; ///< k_F (1/s)
      Real conversionHalf = 3.6e-3; ///< K_G (mol/m^3)

      Real wallFlux = 2.0e-12; ///< j_0 (mol m^-2 s^-1) at E = 1, P = P_ref
      Real activationShearStress = 0.4; ///< Endothelial switch (Pa)
      Real activationShearWidth = 0.15; ///< Width of the switch (Pa)

    /**
     * @brief Thickness of the endocardial reaction layer (m).
     *
     * The wall activation E is extended into the lumen by
     * @f$ -\delta^2\Delta\tilde E+\tilde E=0 @f$, @f$ \tilde E=E @f$ on the
     * wall, and initiation becomes the volumetric channel
     * @f$ s=(j_0/\delta)\tilde E P/P_\mathrm{ref} @f$, scaled so that its
     * volume integral equals the wall dose. Zero restores the wall flux.
     */
      Real layerThickness = 5.0e-4;

      Real stabilizationScale = 1.0; ///< SUPG multiplier of the transport

      struct State
      {
          Real P, T, A, G;
      };

    /**
     * @brief One nodal reaction step, production explicit and destruction
     * implicit (modified Patankar-Euler).
     * @param s Nodal state after transport.
     * @param dt Time increment (s).
     * @param wallRate Layer initiation rate s/P (1/s) at this node.
     * @returns Reacted state; nonnegative, with P + T - A preserved for any dt.
     */
      State advance(const State& s, Real dt, Real wallRate = 0.0) const
      {
      // Transport does not guarantee positivity; undershoots are removed
      // before reacting. The resulting invariant drift is exported.
        const Real P = std::max<Real>(s.P, 0.0), T = std::max<Real>(s.T, 0.0),
                   A = std::max<Real>(s.A, 0.0), G = std::max<Real>(s.G, 0.0);
        const Real Pn = P /
          (1.0 +
            dt *
              (amplificationRate * T * T / ((amplificationHalf + T) * prothrombinRef) +
                wallRate));
        const Real dP = P - Pn;
        const Real a = dt * inhibitionRate;
        const Real b = 1.0 + a * (A + T + dP);
        const Real c = a * A * (T + dP);
        const Real dI =
          2.0 * c / (b + std::sqrt(std::max<Real>(b * b - 4.0 * a * c, 0.0)));
        const Real Tn = T + dP - dI;
        const Real An = A - dI;
        const Real Gn = G / (1.0 + dt * conversionRate * Tn / (conversionHalf + G));
        return {Pn, Tn, An, Gn};
      }
  };

  /// @brief Fibrin-dependent Brinkman drag, sigma_B = c F^2/(K_F^2 + F^2).
  struct BrinkmanDrag
  {
      bool enabled = true;
      Real coefficient = 1.0e7; ///< mu_inf/kappa_0 (kg m^-3 s^-1)
      Real gelThreshold = 2.0e-3; ///< K_F (mol/m^3)

      Real operator()(Real F) const
      {
        if (!enabled || F <= 0.0)
          return 0.0;
        return coefficient * F * F / (gelThreshold * gelThreshold + F * F);
      }
  };
}

#endif
