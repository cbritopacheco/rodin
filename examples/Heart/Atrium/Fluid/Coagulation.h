/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Coagulation.h
 * @brief Constants of the P-T-A-G-X reduced coagulation model and its Brinkman drag.
 *
 * Species (mol/m^3, per plasma volume) obey
 * @f$ \mathcal L c=\partial_t c|_\chi+(u-w)\cdot\nabla c-\nabla\cdot(D\nabla c) @f$
 * on the moving fluid domain, with mesh velocity @f$ w @f$:
 *
 *   L P = -R           prothrombin, the precursor
 *   L T =  R - R_I     thrombin, the enzyme
 *   L A =    - R_I     antithrombin, the stoichiometric inhibitor
 *   L G =      - R_F   fibrinogen, the substrate
 *   L X = -k_dX X      factor Xa, emitted by the activated endocardium
 *   d_t F = R_F        fibrin, immobile on the mesh
 *
 *   R   = (k_a T^2/(K_A + T) + kappa X) P/P_ref   amplification and initiation
 *   R_I = k_on A T                                 1:1 irreversible inhibition
 *   R_F = k_F T G/(K_G + G)                        Michaelis-Menten conversion
 *
 * Factor Xa is a catalyst and is not consumed by @f$ R @f$. Invariants of the
 * continuous problem:
 *   P + T - A = P_0 - A_0,   P + T <= P_0,   A <= A_0,   0 <= G <= G_0,
 * and the reaction conserves G + F. Without factor Xa, initiation is the
 * direct wall channel P -> T; without immobile fibrin, F = G_0 - G is a
 * diagnostic advected with the fibrinogen. I = A_0 - A is the
 * thrombin-antithrombin complex.
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
   * the two calibration constants. The factor Xa emission defaults to the
   * stasis-matched value j_X = j_0 k_dX/kappa, for which X tends to
   * j_X E/(delta k_dX) without flow and kappa X equals the direct wall
   * channel; the two initiations then differ only where the flow washes X out.
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
       * @brief Includes the oscillatory factor in the endothelial activation.
       *
       * True gives E = s(TAWSS) min(2 OSI, 1); false gives E = s(TAWSS), so
       * that unidirectional stasis also activates the wall.
       */
      bool activationUsesOSI = true;

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

      bool factorXa = true; ///< Initiation through transported factor Xa

      /// @brief j_X (mol m^-2 s^-1) at E = 1; negative selects j_0 k_dX/kappa.
      Real xaFlux = -1.0;
      Real xaRate = 3.0; ///< kappa (1/s), X P/P_ref -> T
      Real xaDecay = 1.0e-2; ///< k_dX (1/s), plasma inhibition of Xa

      /// @brief F is a nodal field with d_t F = R_F instead of G_0 - G.
      bool immobileFibrin = true;

      /// @brief Wall emission of factor Xa at E = 1 (mol m^-2 s^-1).
      Real resolvedXaFlux() const
      {
        return xaFlux < 0.0 ? wallFlux * xaDecay / xaRate : xaFlux;
      }

      struct State
      {
          Real P, T, A, G, X, F;
      };

    /**
     * @brief One nodal reaction step, production explicit and destruction
     * implicit (modified Patankar-Euler).
     * @param s Nodal state after transport.
     * @param dt Time increment (s).
     * @param wallRate Layer initiation rate s/P (1/s) of the direct channel.
     * @param xaSource Volumetric emission of factor Xa (mol m^-3 s^-1).
     * @returns Reacted state; nonnegative, with P + T - A and G + F preserved
     *   for any dt.
     *
     * With factor Xa, X* = (X + dt xaSource)/(1 + dt k_dX) is advanced first
     * and kappa X* joins the amplification of the prothrombin step.
     */
      State advance(
        const State& s, Real dt, Real wallRate = 0.0, Real xaSource = 0.0) const
      {
      // Transport does not guarantee positivity; undershoots are removed
      // before reacting. The resulting invariant drift is exported.
        const Real P = std::max<Real>(s.P, 0.0), T = std::max<Real>(s.T, 0.0),
                   A = std::max<Real>(s.A, 0.0), G = std::max<Real>(s.G, 0.0);
        const Real X = std::max<Real>(s.X, 0.0), F = std::max<Real>(s.F, 0.0);
        const Real Xn = factorXa ? (X + dt * xaSource) / (1.0 + dt * xaDecay) : X;
        const Real xaTerm = factorXa ? xaRate * Xn / prothrombinRef : 0.0;
        const Real Pn = P /
          (1.0 +
            dt *
              (amplificationRate * T * T / ((amplificationHalf + T) * prothrombinRef) +
                wallRate + xaTerm));
        const Real dP = P - Pn;
        const Real a = dt * inhibitionRate;
        const Real b = 1.0 + a * (A + T + dP);
        const Real c = a * A * (T + dP);
        const Real dI =
          2.0 * c / (b + std::sqrt(std::max<Real>(b * b - 4.0 * a * c, 0.0)));
        const Real Tn = T + dP - dI;
        const Real An = A - dI;
        const Real Gn = G / (1.0 + dt * conversionRate * Tn / (conversionHalf + G));
        return {Pn, Tn, An, Gn, Xn, F + (G - Gn)};
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
