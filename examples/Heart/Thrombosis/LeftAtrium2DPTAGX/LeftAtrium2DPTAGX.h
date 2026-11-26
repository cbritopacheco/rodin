// LeftAtrium2DPTAGX.h
//
// LeftAtrium2DPTAG plus the switches below (all on by default). With
// -ptag_xa false -ptag_immobile_fibrin false -ptag_fibrin_monomer false
// -ptag_k_tm 0 -ptag_skew false the driver is LeftAtrium2DPTAG, and with
// -ptag_eps_x 1 -ptag_kappa_x 3 -ptag_fibrin_monomer false -ptag_k_tm 0 it is
// the PTAGX of the first round, operation by operation.
//
//   factorXa        initiation through transported factor Xa instead of the
//                   direct wall P -> T channel:
//                     L X = -k_dX X,   -D_P d_n X = j_X E   (endocardium)
//                     R   = (k_a T^2/(K_A + T) + kappa X theta_X) P/P_ref,
//                     theta_X = eps + (1 - eps) T/(K_A + T).
//                   X is factor Xa itself, not the zymogen: the wall flux
//                   j_X = Gamma_TF k_cat is the tissue factor density in flux
//                   units, so TF and FVIIa never appear as fields (FX is at
//                   170 nM and loses < 1 nM over a run, so carrying it would
//                   add a field and a Michaelis constant for nothing).
//                   theta_X is the cofactor Va, which thrombin itself makes:
//                   kappa is the turnover of assembled prothrombinase and eps
//                   the activity without it, so initiation has to bootstrap
//                   its own cofactor and is a threshold, while the wave, where
//                   theta_X -> 1, is not. eps = 1 removes the requirement.
//                   With a reaction layer (delta > 0) the wall emission is the
//                   volumetric source (j_X/delta) chi E~ in the nodal step, with
//                   the same dose normalisation chi as the legacy channel.
//                   X is a catalyst (not consumed), so P + T - A is unchanged.
//                   In stasis X -> j_X chi E~/(delta k_dX), and the default
//                   j_X = j_0 k_dX/kappa makes the two initiations identical
//                   there; they differ only where the flow washes X out.
//   immobileFibrin  F is a nodal field, d_t F = R_F, instead of G_0 - G
//                   (which the flow advects). G + F is conserved by the nodal
//                   step; the Brinkman drag reads this F.
//   fibrinMonomer   the cleavage makes a soluble monomer f that the flow
//                   carries, and only the polymer is a gel:
//                     L f = R_F - R_pol,  d_t F = R_pol,
//                     R_pol = k_pol f (f + F).
//                   G + f + F is conserved. Thrombin can then cleave
//                   fibrinogen away from the lesion without leaving gel
//                   there: at nM the monomer is washed out before it gels, at
//                   uM it gels in seconds. Without flow it changes nothing.
//   thrombomodulin  the endocardium that is NOT activated removes thrombin,
//                   at rate k_TM (1 - E) in the same layer. The same E
//                   therefore raises the tissue factor and lowers the
//                   protection, which is what inflammatory activation does
//                   (Yamashita 2003: rapid pacing lowers thrombomodulin and
//                   TFPI in the atrial endocardium; Nakamura 2003: tissue
//                   factor rises in the atrial endothelium in AF). Thrombin
//                   caught there leaves as a serpin complex, so antithrombin
//                   pays for it and P + T - A stays exact; tfpiLoss does the
//                   TFPI half, k_dX -> k_dX (1 - lambda E~).
//
// -ptag_skew (default true) uses the skew-symmetric advection in the species
// transport; see Config::skewTransport.
//
// ---------------------------------------------------------------------------
// Two-dimensional left atrium with a rectangular appendage: the stabilised
// P1/P1 non-Newtonian flow of LeftAtrium2D, driven by two measured pressure
// waveforms, coupled to a precursor- and inhibitor-conserving reduced
// coagulation model (the "P-T-A-G" model) and to a fibrin-dependent Brinkman
// drag.
//
// Species (mol/m^3, per plasma volume), transported by L c = d_t c + u.grad c
// - div(D grad c):
//
//   L P = -R           prothrombin, the precursor
//   L T =  R - R_I     thrombin, the enzyme
//   L A =    - R_I     antithrombin, the stoichiometric inhibitor
//   L G =      - R_F   fibrinogen, the substrate (FPA-cleavage capacity)
//   L X = -k_dX X      factor Xa, the initiation enzyme (factorXa)
//
//   L f =      - R_pol + R_F   soluble fibrin monomer (fibrinMonomer)
//
//   R   = k_a T^2/(K_A + T) P/P_ref     lumped feedback amplification
//   R_I = k_on A T                       1:1 irreversible inhibition
//   R_F = k_F T G/(K_G + G)              Michaelis-Menten conversion
//
// Initiation (endocardium): E in [0,1] is the endothelial activation built
// from the cycle wall-shear indices, extended into an endocardial reaction
// layer of thickness delta (the trabeculated wall) as E~; the P -> T channel
// s = (j_0/delta) E~ P/P_ref carries the wall dose j_0 int E dA into the
// layer. With delta = 0 it reduces to the wall flux D_T d_n T = j_0 E P/P_ref,
// D_P d_n P = -j. Diagnostics: F = G_0 - G (pre-gel fibrin surrogate) and
// I = A_0 - A (thrombin-antithrombin complex).
//
// Why this and not a thrombin/fibrinogen/fibrin triple with a constant
// source: every molecule produced is taken from a tracked pool, so
//   P + T - A = P_0 - A_0,   P + T <= P_0,   A <= A_0,   0 <= G <= G_0
// are invariant regions of the continuous problem and nothing can grow without
// bound; R'(0) = 0, so a small contamination decays at rate k_on A instead of
// being amplified; and the thrombin burst ends by precursor exhaustion, as a
// thrombogram does.
//
// Flow feedback: the momentum equation carries -sigma_B(F) u with
// sigma_B = (mu_inf/kappa_0) F^2/(K_F^2 + F^2) >= 0, assembled implicitly with
// F lagged one step. The walls are rigid, so the gel scaffold is at rest and
// the term is exactly dissipative: -int sigma_B |u|^2 <= 0 at every step.
//
// Time integration: backward Euler for the flow and for the four transport
// equations (velocity lagged, SUPG + Codina crosswind, Robin term of P
// implicit, the T load using the P just solved so the wall exchange is
// conservative to round-off), followed by a nodal production-destruction step
// that is positive and budget-exact for any dt (modified Patankar-Euler). Lie
// splitting: with a first-order transport step, Strang would not raise the
// order.
//
// Equal order P1/P1, VMS (convection, grad-div) and PSPG. Runs under MPI.
#ifndef EXAMPLES_HEART_LEFTATRIUM2DPTAGX_H
#define EXAMPLES_HEART_LEFTATRIUM2DPTAGX_H

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/MPI.h>
#include <Rodin/PETSc.h>
#include <Rodin/Solver.h>
#include <Rodin/Types.h>
#include <Rodin/Variational.h>

#include "CoronaryArtery/ThrombosisModel.h"
#include "CoronaryArtery/VMSConvectionIntegrator.h"

namespace Rodin::Examples::Heart
{
  /// @brief Tabulated, periodic pressure waveform: two columns t [s], p [Pa],
  ///        '#' comments, linear interpolation, periodic extension.
  class PeriodicWaveform
  {
    public:
      using Real = Rodin::Real;

      void load(const std::string& path);
      Real operator()(Real t) const;

      Real getPeriod() const noexcept { return m_period; }
      Real getMinimum() const noexcept { return m_min; }
      Real getMaximum() const noexcept { return m_max; }
      const std::string& getPath() const noexcept { return m_path; }

    private:
      std::string m_path;
      std::vector<Real> m_t, m_p;
      Real m_period = 0.0, m_min = 0.0, m_max = 0.0;
  };

  /**
   * @brief Constants of the P-T-A-G kinetics (SI: mol, m, s).
   *
   * @details Literature values unless marked. 1 uM = 1e-3 mol/m^3.
   *          k_on = 7.1e3 M^-1 s^-1 = 7.1 m^3 mol^-1 s^-1, so
   *          k_on A_0 = 0.020 s^-1 (thrombin half-life ~35 s in plasma).
   *          k_F, K_G are Higgins et al. (1983) in fibrinogen-equivalent units
   *          (two A-alpha chains per molecule). K_A is fixed by the frozen-pool
   *          activation threshold T_th = k_on A_0 K_A/(k_a P_0/P_ref - k_on A_0)
   *          = 2.5 nM. k_a and j_0 are the two calibration constants.
   */
  struct ConservingKinetics
  {
      using Real = Rodin::Real;

      Real diffusivityProtein = 4.6e-11;      ///< D_P = D_T = D_A (m^2/s)
      Real diffusivityFibrinogen = 2.0e-11;   ///< D_G (m^2/s)

      Real prothrombin0 = 1.4e-3;             ///< P_0, inflow/initial (mol/m^3)
      Real antithrombin0 = 2.8e-3;            ///< A_0 (mol/m^3)
      Real fibrinogen0 = 7.35e-3;             ///< G_0: 2.5 g/L (SR); 4 g/L = 11.76e-3 (AF)
      Real prothrombinRef = 1.4e-3;           ///< P_ref, fixed across cases
      Real inletThrombin = 0.0;               ///< T at the veins (mol/m^3)

      Real amplificationRate = 0.10;          ///< k_a (1/s), calibrate
      Real amplificationHalf = 1.0e-5;        ///< K_A (mol/m^3), 10 nM
      Real inhibitionRate = 7.1;              ///< k_on (m^3 mol^-1 s^-1)
      Real conversionRate = 42.0;             ///< k_F (1/s)
      Real conversionHalf = 3.6e-3;           ///< K_G (mol/m^3)

      Real wallFlux = 2.0e-12;                ///< j_0 (mol m^-2 s^-1) at E = 1, P = P_ref
      Real activationShearStress = 0.4;       ///< endothelial switch (Pa)
      Real activationShearWidth = 0.15;       ///< width of the switch (Pa)
      /// @brief E = sigma(TAWSS) min(2 OSI, 1) if true (LeftAtrium2DPTAG);
      ///        E = sigma(TAWSS) if false, so that steady stasis (AF) is not
      ///        discounted for being unidirectional.
      bool activationUsesOSI = true;
      /// @brief Prescribed activation: E = this on the appendage labels and 0
      ///        on the other walls, instead of the shear indices. Negative
      ///        keeps the shear-based E. Used by the benchmarks, where the
      ///        activator is a chosen piece of wall.
      Real activationFixed = -1.0;

      /**
       * @brief Thickness of the endocardial reaction layer (m); 0 restores
       *        the wall-flux (Robin) form.
       *
       * @details The atrial endocardium is trabeculated (pectinate muscles,
       *          0.5-1 mm), so the activated surface is a layer, not a line.
       *          The activation E is extended into the lumen by
       *            -delta^2 Laplace(E~) + E~ = 0,  E~ = E on the wall,
       *          i.e. E~ ~ E exp(-d/delta), and the initiation becomes the
       *          volumetric P -> T channel s = (j_0/delta) E~ P/P_ref, scaled
       *          so that int s dV equals the wall dose j_0 int E dA exactly.
       *          With a Robin flux on no-slip nodes the products can leave
       *          the wall only by molecular diffusion (h^2/D ~ 2000 s per
       *          cell); in the layer they are born where the fluid moves.
       */
      Real layerThickness = 5.0e-4;

      Real stabilizationScale = 1.0;          ///< SUPG multiplier

      // ---- Factor Xa (LeftAtrium2DPTAGX) ----------------------------------
      bool factorXa = true;                   ///< initiation through transported X
      /// @brief j_X (mol m^-2 s^-1) at E = 1; < 0 selects the stasis-matched
      ///        value j_0 k_dX/kappa (resolved in initialize()).
      Real xaFlux = -1.0;
      Real xaRate = 3.0;                      ///< kappa (1/s): X P/P_ref -> T
      Real xaDecay = 1.0e-2;                  ///< k_dX (1/s): plasma inhibition of Xa
      /// @brief F as an immobile nodal field (d_t F = R_F) instead of G_0 - G.
      bool immobileFibrin = true;
      /// @brief sigma (m): procoagulant membrane per unit area of activated
      ///        endocardium. In the layer, theta = 1 + (sigma/delta) chi E~,
      ///        so that int (theta - 1) dV = sigma int E dA. 0 = off.
      Real surfaceAmplification = 0.0;

      /**
       * @brief Soluble fibrin monomer f between fibrinogen and the gel.
       *
       * @details Thrombin cleaves fibrinogen into a monomer that is still
       *          soluble and is carried by the flow; only the polymer is a
       *          gel. R_pol = k_pol f (f + F) (nucleation and growth on the
       *          fibres already there, one constant: Weisel & Nagaswami 1992
       *          give ~1e6 1/(M s)). Monomer made at nM is washed out of a
       *          swept region before it gels, at uM it gels in seconds, so the
       *          cleavage can happen away from the lesion without leaving gel
       *          there. With it off, every cleaved molecule gels on the spot.
       *          G + f + F is conserved by the nodal step either way.
       */
      bool fibrinMonomer = true;
      Real polymerizationRate = 1.0e3;        ///< k_pol (m^3 mol^-1 s^-1)

      /**
       * @brief k_TM (m/s): thrombin removal by the thrombomodulin of the
       *        endocardium that is NOT activated.
       *
       * @details Inflammatory activation of the endothelium raises tissue
       *          factor and lowers thrombomodulin, EPCR, TFPI and eNOS at the
       *          same time. The same E therefore does two things here: it
       *          turns the tissue-factor source on, j_X E, and it turns the
       *          protection off, k_TM (1 - E). Thrombin caught by
       *          thrombomodulin is inactivated by antithrombin about twenty
       *          times faster, so it leaves as a serpin complex: the nodal
       *          step pays for it with antithrombin and P + T - A stays exact,
       *          which also bounds the protection by the inhibitor pool.
       *          k_TM ~ Gamma_TM k_on ~ 1e-6 - 1e-5 m/s for 1e2 - 1e3
       *          thrombomodulin/um^2 and a near diffusion-limited k_on.
       *          0 = off (the wall has no protection anywhere).
       */
      Real thrombomodulinRate = 1.0e-5;

      /**
       * @brief lambda in [0, 1): share of the inhibition of Xa that the
       *        activated endocardium loses with its TFPI.
       *
       * @details Rapid atrial pacing downregulates thrombomodulin AND tissue
       *          factor pathway inhibitor in the atrial endocardium (Yamashita
       *          et al., Circulation 108 (2003) 2450), so the same E that
       *          raises the tissue factor also lets the Xa it makes live
       *          longer: k_dX -> k_dX (1 - lambda E~) in the layer. Of the
       *          ~1e-2 1/s at which plasma inhibits Xa, antithrombin accounts
       *          for ~6e-3 and TFPI for the rest, so lambda <= ~0.4. Default 0:
       *          the effect is a sharpening of j_X E, which is already there,
       *          and it is reported as a sensitivity instead.
       */
      Real tfpiLoss = 0.0;

      /**
       * @brief epsilon: activity of Xa that is NOT bound to the cofactor Va,
       *        relative to assembled prothrombinase.
       *
       * @details Free Xa is a poor prothrombin activator; the Xa.Va complex on
       *          a phospholipid membrane is orders of magnitude faster, and the
       *          Va comes from thrombin itself. So kappa is the turnover of the
       *          assembled complex (~10-30 1/s) and the Xa channel carries the
       *          cofactor that is available at the thrombin there,
       *            theta_X = epsilon + (1 - epsilon) T/(K_A + T),
       *          with the same half-saturation as the other thrombin feedback.
       *          epsilon = 1 is the model without a cofactor requirement (what
       *          kappa = 3 1/s lumped). Below 1 initiation has to bootstrap its
       *          own cofactor, so it becomes a threshold: that is what makes
       *          the lag sensitive to the tissue-factor dose and to the flow
       *          while the propagating wave, where theta_X -> 1, is not.
       */
      Real xaCofactorFloor = 1.0;          ///< epsilon (-)


      struct State { Real P, T, A, G, X, F, Fm; };

      /**
       * @brief One nodal reaction step, production explicit / destruction
       *        implicit (modified Patankar-Euler).
       *
       * @param wallRate The legacy initiation rate s/P (1/s) at this node.
       * @param xaSource Volumetric emission of X (mol m^-3 s^-1) at this node.
       * @details With factorXa: X* = (X + dt xaSource)/(1 + dt k_dX) first,
       *          and kappa X* joins the amplification in P*. F* = F + G - G*.
       * @details P* = P/(1 + dt [k_a T^2/((K_A+T) P_ref) + wallRate]),  dP = P - P*;
       *          dI = smaller root of dI = dt k_on (A - dI)(T + dP - dI);
       *          T* = T + dP - dI,  A* = A - dI,
       *          G* = G/(1 + dt k_F T/(K_G + G)).
       *          For any dt: all values stay >= 0, P* + T* - A* = P + T - A
       *          exactly, and the prothrombin consumed equals the thrombin
       *          produced.
       */
      State advance(const State& s, Real dt, Real wallRate = 0.0,
                    Real xaSource = 0.0, Real surfaceTheta = 0.0,
                    Real tmRate = 0.0, Real xaDecayLoss = 0.0) const
      {
        // The transport step (SUPG) does not guarantee positivity; an
        // undershoot is removed here, before it can react. The price is a
        // drift of P + T - A that the CSV reports, so it is never silent.
        const Real P = std::max<Real>(s.P, 0.0), T = std::max<Real>(s.T, 0.0),
                   A = std::max<Real>(s.A, 0.0), G = std::max<Real>(s.G, 0.0);
        const Real X = std::max<Real>(s.X, 0.0), F = std::max<Real>(s.F, 0.0),
                   Fm = std::max<Real>(s.Fm, 0.0);

        const Real Xn = factorXa
          ? (X + dt * xaSource) / (1.0 + dt * std::max<Real>(xaDecay - xaDecayLoss, 0.0)) : X;
        const Real thetaX = xaCofactorFloor
                          + (1.0 - xaCofactorFloor) * T / (amplificationHalf + T);
        const Real xaTerm = factorXa ? xaRate * Xn * thetaX / prothrombinRef : 0.0;
        const Real Pn = P / (1.0 + dt * (amplificationRate * (1.0 + surfaceTheta) * T * T /
                                     ((amplificationHalf + T) * prothrombinRef) + wallRate
                                     + xaTerm));
        const Real dP = P - Pn;

        const Real a = dt * inhibitionRate;
        const Real b = 1.0 + a * (A + T + dP);
        const Real c = a * A * (T + dP);
        const Real dI = 2.0 * c / (b + std::sqrt(std::max<Real>(b * b - 4.0 * a * c, 0.0)));

        Real Tn = T + dP - dI;
        Real An = A - dI;
        // Thrombomodulin of the intact endocardium: the thrombin it catches
        // ends in a serpin complex, so antithrombin pays for it and
        // P + T - A is still exact. Capped by the inhibitor left.
        if (tmRate > 0.0)
        {
          const Real dTM = std::min<Real>(An, Tn * dt * tmRate / (1.0 + dt * tmRate));
          Tn -= dTM;
          An -= dTM;
        }

        const Real Gn = G / (1.0 + dt * conversionRate * Tn / (conversionHalf + G));
        Real Fn, Fmn = 0.0;
        if (immobileFibrin && fibrinMonomer)
        {
          const Real Fmt = Fm + (G - Gn);
          Fmn = Fmt / (1.0 + dt * polymerizationRate * (Fmt + F));
          Fn = F + (Fmt - Fmn);
        }
        else
          Fn = F + (G - Gn);
        return { Pn, Tn, An, Gn, Xn, Fn, Fmn };
      }
  };

  /// @brief Fibrin-dependent Brinkman drag, sigma_B = c F^2/(K_F^2 + F^2).
  struct BrinkmanDrag
  {
      using Real = Rodin::Real;
      bool enabled = true;
      Real coefficient = 1.0e7;   ///< mu_inf/kappa_0 (kg m^-3 s^-1)
      Real gelThreshold = 2.0e-3; ///< K_F (mol/m^3)

      Real operator()(Real F) const
      {
        if (!enabled || F <= 0.0)
          return 0.0;
        return coefficient * F * F / (gelThreshold * gelThreshold + F * F);
      }
  };

  class LeftAtrium2DPTAGX
  {
    public:
      using Real = Rodin::Real;
      using Attribute = Rodin::Geometry::Attribute;
      using MeshType = Rodin::Geometry::Mesh<Rodin::Context::MPI>;
      using AttributeSet = Rodin::FlatSet<Attribute>;

      using VectorFESType =
        Rodin::Variational::H1<1, Rodin::Math::SpatialVector<Real>, MeshType>;
      using ScalarFESType = Rodin::Variational::H1<1, Real, MeshType>;

      using VectorGridFunctionType =
        Rodin::PETSc::Variational::GridFunction<VectorFESType>;
      using ScalarGridFunctionType =
        Rodin::PETSc::Variational::GridFunction<ScalarFESType>;

      using VectorTrialFunctionType = Rodin::PETSc::Variational::TrialFunction<
        VectorGridFunctionType, VectorFESType>;
      using ScalarTrialFunctionType = Rodin::PETSc::Variational::TrialFunction<
        ScalarGridFunctionType, ScalarFESType>;
      using VectorTestFunctionType =
        Rodin::PETSc::Variational::TestFunction<VectorFESType>;
      using ScalarTestFunctionType =
        Rodin::PETSc::Variational::TestFunction<ScalarFESType>;

      using LinearSystemType = Rodin::PETSc::Math::LinearSystem;

      using FlowProblemType = Rodin::Variational::Problem<LinearSystemType,
        VectorTrialFunctionType, ScalarTrialFunctionType,
        VectorTestFunctionType, ScalarTestFunctionType>;
      using VectorProblemType = Rodin::Variational::Problem<LinearSystemType,
        VectorTrialFunctionType, VectorTestFunctionType>;
      using ScalarProblemType = Rodin::Variational::Problem<LinearSystemType,
        ScalarTrialFunctionType, ScalarTestFunctionType>;
      using FluxFormType = Rodin::Variational::LinearForm<ScalarFESType, ::Vec>;

      /// @brief Pointwise scalar coefficient (never L2-projected: the
      ///        stabilisation parameters and the cycle indices are nodal
      ///        quantities, and a projection can make them negative).
      using ScalarCoefficientType = Rodin::Variational::RealFunction<
        std::function<Real(const Rodin::Geometry::Point&)>>;

      struct CarreauYasuda
      {
          Real mu0 = 0.0720;
          Real muInf = 0.0049;
          Real lambda = 6.1529;
          Real n = 0.2081;
          Real yasuda = 1.4173;
          Real gammaRegularization = 1.0e-3;
      };

      /// @brief Boundary labels of LA2D_rectLAA.mesh (see LeftAtrium2D.h).
      struct Labels
      {
          std::array<Attribute, 4> inlets{{6, 8, 10, 12}};
          Attribute outlet = 14;
          std::array<Attribute, 6> wall{{1, 5, 7, 9, 11, 13}};
          std::array<Attribute, 3> appendage{{2, 3, 4}};
      };

      struct Config
      {
          std::string meshPath = "../resources/examples/Heart/LA2D_rectLAA_2D.mesh";
          std::string inletPressurePath = "../resources/examples/Heart/presion_PV_SR.dat";
          std::string outletPressurePath = "../resources/examples/Heart/presion_MV_SR.dat";
          std::string xdmfBasename = "LeftAtrium2DPTAGX";
          std::string csvPath = "LeftAtrium2DPTAGX.csv";
          /// @brief Box xmin, xmax, ymin, ymax (m, after meshScale) for the
          ///        appendage integrals of the CSV; xmin >= xmax disables them.
          std::array<Real, 4> laaBox{{0.0, 0.0, 0.0, 0.0}};
          /// @brief File with the appendage polygon (one "x y" per line, m);
          ///        replaces laaBox in the CSV integrals when given.
          std::string laaPolygonPath;

          // ---- Benchmark set-ups (Qureshi et al., CinC 2021) ---------------
          /// @brief Transpiration on the appendage walls instead of no-slip:
          ///        u = -W sin(2 pi t/T + phase) n (n outward), imposed by
          ///        penalty. W > 0 enables it (linearised wall contraction on
          ///        a fixed mesh).
          /// @brief Prescribed inflow speed (m/s) at the veins instead of the
          ///        pressure boundary: u = -V(t) n, with V(t) the inlet
          ///        waveform rescaled so that its maximum is this value.
          ///        <= 0 keeps the pressure boundary.
          Real inletVelocity = 0.0;
          Real laaWallVelocity = 0.0;
          Real laaWallPhase = 0.0;
          Real wallPenalty = 1.0e5;   ///< Pa s/m
          /// @brief Thrombin seed x, y, radius (m), value (mol/m^3) at t = 0;
          ///        radius <= 0 disables it. Smooth edge over 20 % of r.
          std::array<Real, 4> seed{{0.0, 0.0, 0.0, 0.0}};

          Labels labels;
          Real meshScale = 1.0;

          Real rho = 1060.0;
          CarreauYasuda viscosity;
          Real pressurePenalty = 1.0e-12;
          Real pressureOffset = 0.0;

          /// @brief Inlet regularisations (Pa s/m); see LeftAtrium2D.h for
          ///        their sizing against the 30 Pa driving head.
          Real inletImpedance = 20.0;
          Real inletTangentialDamping = 50.0;
          /// @brief Directional do-nothing weight; 1 is the only value with an
          ///        unconditional energy estimate.
          Real backflowStabilization = 1.0;

          Real dt = 1.0e-2;
          Real period = 0.8;
          int flowCycles = 2;
          int speciesCycles = 10;
          int outputEvery = 20;
          /// @brief CSV row (and diagnostics) every csvEvery clock steps, and
          ///        at the end of every cycle.
          int csvEvery = 1;

          Real vmsScale = 1.0;
          Real gradDivScale = 1.0;
          Real pspgScale = 1.0;
          bool useVMS = true;
          /// @brief Quasi-static (OSS) instead of dynamic velocity subscales.
          bool quasiStaticSubscales = false;
          /// @brief Use sqrt(tau_C(u^{n-1})) in the grad-div projection, the
          ///        tau that was implicit on the previous step (telescopic).
          bool lagGradDivTau = false;
          /// @brief Weight of the transient, convective and drag parts of the
          ///        PSPG residual: 1 consistent, 0 pressure-gradient penalty.
          Real pspgResidualScale = 1.0;
          Real crosswindC = 0.7;

          /// @brief GMRES + block-Jacobi/ILU for the four transport solves
          ///        instead of the global direct solver.
          bool iterativeTransport = false;
          /// @brief Skew-symmetric advection about the inflow state in the
          ///        species transport: adds (1/2) int div(u_h) (c - c_in) v.
          /// @details P1/P1 VMS leaves div u_h = O(10) 1/s at the short wall
          ///          gaps between two inlets. With the plain u.grad c form a
          ///          deviation from the inflow state grows there at rate
          ///          ~|div u_h|/2 (SR, coarse mesh: P and A depleted to 15 %,
          ///          G to 12 %, so G_0 - G read as fibrin, within 10 cycles);
          ///          the plain skew form instead forces c_in by div u_h/2
          ///          (P to 3.6 P_0). About c_in, constants stay exact and
          ///          deviations are not amplified. false restores
          ///          LeftAtrium2DPTAG.
          bool skewTransport = true;

          bool solveKinetics = true;
          ConservingKinetics kinetics;
          BrinkmanDrag brinkman;

          /// @brief Cycle replay. After the warm-up the last flow cycle is
          ///        stored and replayed for the species; the flow is solved
          ///        again only when the fibrin has changed enough to alter it.
          /// @details Exact while sigma_B = 0 (the flow is periodic and the
          ///          species do not act on it); afterwards the error is the
          ///          one-cycle lag of a field that evolves over tens of
          ///          cycles. Cost per step: four scalar transports.
          bool replay = false;
          /// @brief Flow cycles solved at each refresh (the last one is stored).
          int replayFlowCycles = 2;
          /// @brief Refresh when max F > K_F/2 and the fibrin mass has grown by
          ///        this fraction since the stored cycle.
          Real replayRefresh = 0.2;
          /// @brief Advance the species every k flow steps with dt_q = k dt.
          /// @details The kinetics have time constants of 10-50 s and the
          ///          transport is implicit, so k = 2-4 is a time-step study,
          ///          not a stability question. Must divide the steps per
          ///          cycle. Applies in every phase, replay or not.
          int speciesSubsteps = 1;

          /// @brief Divergence guard on max|u| (m/s).
          Real maxVelocity = 1.0;
      };

      LeftAtrium2DPTAGX(const Rodin::Context::MPI& context, const Config& cfg);
      ~LeftAtrium2DPTAGX();

      LeftAtrium2DPTAGX(const LeftAtrium2DPTAGX&) = delete;
      LeftAtrium2DPTAGX& operator=(const LeftAtrium2DPTAGX&) = delete;

      LeftAtrium2DPTAGX& initialize();
      int run();

      Config& getConfig() noexcept { return m_cfg; }

    private:
      static MeshType makeMesh(const Rodin::Context::MPI& context, const Config& cfg);
      static Real cellSize(const Rodin::Geometry::Point& p);
      static void axpy(Real a, const ::Vec& x, ::Vec& y);

      bool isRoot() const;

      void setupSpaces();
      void setupFlow();
      void setupSpecies();
      void setupWallShear();
      void setupLayer();
      /// @brief Extend E into the endocardial layer and fix the dose scale.
      void extendActivation();

      /// @brief Backward-Euler SUPG transport of one species with the lagged
      ///        velocity, a Robin wall term robinScale*E*c and a wall load
      ///        robinScale*E*load, both zero for the species without wall
      ///        exchange.
      void assignTransport(ScalarProblemType& problem, ScalarTrialFunctionType& c,
        ScalarTestFunctionType& v, ScalarGridFunctionType& cur,
        ScalarGridFunctionType& prev, Real diffusivity, Real inlet,
        Real robinScale, Real loadScale, const ScalarGridFunctionType* load);
      /// @brief int_Omega f dV, or over the appendage box when inBox.
      Real integrate(const std::function<Real(const Rodin::Geometry::Point&)>& f,
        bool inBox);
      bool inLAA(const Rodin::Geometry::Point& p) const;

      template <class Expression>
      void project(const Expression& expr, ScalarGridFunctionType& out)
      {
        m_scalarProjection = Rodin::Variational::Integral(m_sTrial, m_sTest) -
          Rodin::Variational::Integral(expr, m_sTest);
        m_scalarProjection.solve(m_scalarProjectionKSP);
        out.setData(m_sTrial.getSolution().getData());
      }

      template <class Expression>
      void projectVector(const Expression& expr, VectorGridFunctionType& out)
      {
        m_vectorProjection = Rodin::Variational::Integral(m_wTrial, m_wTest) -
          Rodin::Variational::Integral(expr, m_wTest);
        m_vectorProjection.solve(m_vectorProjectionKSP);
        out.setData(m_wTrial.getSolution().getData());
      }

      bool solveFlow();
      /// @brief Copy u^{n+1} into slot k of the stored cycle.
      void storeSnapshot(int k);
      /// @brief Load slot k (and k-1) as the lagged velocities.
      void loadSnapshot(int k);
      /// @brief int_Omega F dV, the refresh monitor.
      Real fibrinMass();
      void solveSpecies();
      void reactionStep();
      void updateDiagnostics();
      void computeWallShear();
      void closeCycle(Real elapsed);
      void computeFluxes();
      Real boundaryMeasure(const AttributeSet& tags);

      Real viscosityAt(const Rodin::Geometry::Point& p,
        const VectorGridFunctionType& u) const;
      /// @brief Fibrin at a point: the immobile field, or G_0 - G clamped
      ///        to [0, G_0].
      Real fibrinAt(const Rodin::Geometry::Point& p) const;
      /// @brief sigma_B(F) at a point, from the lagged fibrinogen.
      Real dragAt(const Rodin::Geometry::Point& p) const;
      /// @brief Codina tau_1 with the Brinkman reaction, on the local viscosity.
      Real tau1At(const Rodin::Geometry::Point& p, const VectorGridFunctionType& u) const;
      Real vmsTauAt(const Rodin::Geometry::Point& p) const;
      Real sqrtTauCAt(const Rodin::Geometry::Point& p) const;
      Real sqrtTauCAt(const Rodin::Geometry::Point& p, const VectorGridFunctionType& u) const;
      Real tauCAt(const Rodin::Geometry::Point& p) const;
      Real tauPAt(const Rodin::Geometry::Point& p) const;
      Real crosswind(const Rodin::Geometry::Point& p, Real cur, Real prev,
        const Rodin::Math::SpatialVector<Real>& gradient, Real diffusivity) const;

      void writeCSVHeader();
      void writeCSVRow(int cycle);

      Config m_cfg;
      PeriodicWaveform m_inletWave, m_outletWave;

      MeshType m_mesh;
      Rodin::IO::XDMF m_xdmf;
      AttributeSet m_inletSet, m_wallSet;
      /// @brief Walls with no-slip, and appendage walls (transpiration).
      AttributeSet m_noSlipSet, m_appendageSet;
      Real m_laaWall = 0.0, m_inletSpeed = 0.0;
      std::vector<std::array<Real, 2>> m_laaPoly;
      Real laaWallAt(Real t) const;

      VectorFESType m_vh;
      ScalarFESType m_sh;

      // ---- Flow ------------------------------------------------------------
      VectorTrialFunctionType m_u;
      ScalarTrialFunctionType m_p;
      VectorTestFunctionType m_v;
      ScalarTestFunctionType m_q;
      VectorGridFunctionType m_uOld, m_uOldOld;

      ScalarTrialFunctionType m_sTrial;
      ScalarTestFunctionType m_sTest;
      VectorTrialFunctionType m_wTrial;
      VectorTestFunctionType m_wTest;

      ScalarCoefficientType m_tauFn, m_tauCFn, m_sqrtTauCFn, m_sqrtTauCOldFn,
        m_tauPFn, m_dragFn, m_pspgDragFn;
      ScalarGridFunctionType m_piTilde;
      VectorGridFunctionType m_convProjection, m_sub, m_subOld;

      // ---- Species: P, T, A, G ---------------------------------------------
      ScalarTrialFunctionType m_pT, m_tT, m_aT, m_gT;
      ScalarTestFunctionType m_vP, m_vT, m_vA, m_vG;
      ScalarGridFunctionType m_pCur, m_tCur, m_aCur, m_gCur;
      ScalarGridFunctionType m_pPrev, m_tPrev, m_aPrev, m_gPrev;
      ScalarGridFunctionType m_pNext, m_tNext, m_aNext, m_gNext;
      /// @brief Diagnostics F = G_0 - G, I = A_0 - A, sigma_B and the drift
      ///        of the invariant P + T - A - (P_0 - A_0).
      ScalarGridFunctionType m_fibrin, m_tat, m_drag, m_drift;

      // ---- Factor Xa and immobile fibrin -------------------------------------
      ScalarTrialFunctionType m_xT;
      ScalarTestFunctionType m_vX;
      ScalarGridFunctionType m_xCur, m_xPrev, m_xNext;
      /// @brief Soluble fibrin monomer (transported, used when fibrinMonomer).
      ScalarTrialFunctionType m_fmT;
      ScalarTestFunctionType m_vFm;
      ScalarGridFunctionType m_fmCur, m_fmPrev, m_fmNext;
      /// @brief Immobile fibrin (used when kinetics.immobileFibrin).
      ScalarGridFunctionType m_fCur, m_fNext;
      /// @brief Scratch field for the CSV integrals.
      ScalarGridFunctionType m_work;

      // ---- Endocardial layer ------------------------------------------------
      ScalarTrialFunctionType m_eT;
      ScalarTestFunctionType m_vE;
      /// @brief E~, the activation extended into the lumen; zero if the layer
      ///        thickness is zero.
      ScalarGridFunctionType m_layer;
      /// @brief Extension of 1: the layer itself. The thrombomodulin of the
      ///        intact wall lives on L~ - E~, the extension of 1 - E; note
      ///        that 1 - E~ is not a layer, it is ~1 in the whole lumen.
      ScalarGridFunctionType m_layerOne;
      /// @brief Dose normalisation, delta int_Gw E dA / int_Omega E~ dV.
      Real m_layerScale = 1.0;
      /// @brief The same for 1 - E: delta int_Gw (1-E) dA / int (L~ - E~) dV.
      Real m_protectionScale = 1.0;

      // ---- Wall shear indices ----------------------------------------------
      VectorGridFunctionType m_wss, m_symRec0, m_symRec1, m_netShear;
      ScalarGridFunctionType m_absShear, m_shearMagnitude, m_tawss, m_osi,
        m_activation;

      // ---- Fluxes ----------------------------------------------------------
      ScalarTestFunctionType m_qFlux;
      ScalarGridFunctionType m_one;
      FluxFormType m_flux;

      // ---- Problems and solvers --------------------------------------------
      FlowProblemType m_flow;
      Rodin::Solver::KSP m_flowKSP;
      ScalarProblemType m_transportP, m_transportT, m_transportA, m_transportG,
        m_transportX, m_transportFm;
      Rodin::Solver::KSP m_kspP, m_kspT, m_kspA, m_kspG, m_kspX, m_kspFm;
      ScalarProblemType m_extension;
      Rodin::Solver::KSP m_extensionKSP;
      ScalarProblemType m_scalarProjection;
      Rodin::Solver::KSP m_scalarProjectionKSP;
      VectorProblemType m_vectorProjection;
      Rodin::Solver::KSP m_vectorProjectionKSP;
      VectorTrialFunctionType m_wssTrial;
      VectorTestFunctionType m_wssTest;
      VectorProblemType m_wssProjection;
      Rodin::Solver::KSP m_wssKSP;

      /// @brief One stored cycle of u, one ghosted Vec per step.
      std::vector<::Vec> m_cycleU;

      /// @brief Species time step, speciesSubsteps * dt.
      Real m_dtChem = 0.0;
      Real m_t = 0.0, m_pIn = 0.0, m_pOut = 0.0;
      Real m_outletMeasure = 1.0, m_inletMeasure = 1.0;
      Real m_qIn = 0.0, m_qOut = 0.0, m_speed = 0.0, m_velocityScale = 0.0;
      bool m_flowFieldSplitsSet = false;
      bool m_initialized = false;
      bool m_indicesReady = false;
      int m_step = 0;

      std::ofstream m_csv;
  };
}

#endif
