/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file PTAGXBench.h
 * @brief P-T-A-G coagulation model with an immobile fibrin gel, a
 *        thrombomodulin wall sink and three platelet states, transported by
 *        a prescribed velocity (SUPG + Codina crosswind). Final form of the
 *        model discussed in the notes of November 2026.
 *
 * Species (mol/m^3 for proteins, 1/m^3 for platelets). Transported, with
 * L c = d_t c + u.grad c - div(D grad c):
 *
 *   L P = -R theta                  prothrombin
 *   L T =  R theta - R_I            thrombin
 *   L A =          - R_I            antithrombin
 *   L G =          - R_F            fibrinogen
 *   L N = -R_act                    resting platelets
 *   L M =  R_act - R_agg            activated mobile platelets
 *
 * Immobile (nodal ODEs):
 *
 *   d_t F = R_F                     fibrin gel; legacy option F = G_0 - G
 *   d_t S = R_agg (+ wall adhesion)  deposited platelets
 *
 * Rates:
 *
 *   R     = k_a T^2/(K_A + T) P/P_ref            amplification (prothrombinase)
 *   R_I   = k_on A T                             inhibition (antithrombin)
 *   R_F   = k_F T G/(K_G + G) [1 - F/F_cap]_+    polymerisation into the gel
 *   R_act = k_p min(1, T/(T_p + T) + S~/S_max) N  activation: thrombin or contact
 *   R_agg = k_agg (S~/S_max)(1 - S/S_max) M      aggregation onto the deposit
 *   theta = theta_0 + (M + S)/N_0                procoagulant membrane
 *   S~    = (I - l^2 Laplace)^-1 S               deposit seen one platelet away
 *
 * The stabilised core S_s = S F^2/(K_F^2 + F^2) is a diagnostic only (it no
 * longer enters theta). With the platelets switched off, theta = theta_0
 * (theta_0 = 1 is the original P-T-A-G model, i.e. a CAT with phospholipid).
 *
 * Wall conditions. Three wall phenotypes, each an attribute set:
 *
 *   reactive (tissue factor, activation E):
 *     D_P d_n P = -beta E P,  D_T d_n T = +beta E P,     beta = j_0/P_ref,
 *     D_pl d_n M = -k_w E (1 - S/S_max) M;
 *   collagen (exposed subendothelium):
 *     D_pl d_n M = -k_w (1 - S/S_max) M,
 *     D_pl d_n N = -k_wN (1 - S/S_max) N;
 *   thrombomodulin (intact endothelium):
 *     D_T d_n T = -k_TM T.
 *
 * What the platelets lose at the wall is added to S at the wall nodes, so
 * N + M + S is conserved; what P loses at the wall enters T, so P + T - A is
 * conserved except for the thrombomodulin loss, which is reported.
 *
 * Time integration: backward Euler for the six transports (velocity and
 * coefficients lagged, SUPG on the transient and convective residual, Codina
 * crosswind on the full residual of the previous step), then the wall
 * deposition, the screened-Poisson smoothing of S and a nodal modified
 * Patankar-Euler step that is positive and conserving for any dt. Lie
 * splitting.
 */
#ifndef EXAMPLES_HEART_THROMBOSIS_PTAGXBench_H
#define EXAMPLES_HEART_THROMBOSIS_PTAGXBench_H

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

namespace Rodin::Examples::Heart::Thrombosis
{
  /**
   * @brief Constants of the P-T-A-G + platelet model (SI: mol, m, s).
   *
   * @details Protein values as in AtriumPTAG.h (calibrated on the CAT).
   *          Platelet values from Wu et al. (2017, Sci. Rep. 7:42720) unless
   *          marked: the activation time 1/k_p ~ 1 s and the thrombin level
   *          T_p of their agonist threshold, the packing limit
   *          S_max = PLT_s,max/d_p, the adhesion velocities, the aggregation
   *          velocity k_aa divided by the platelet diameter, and the
   *          shear-enhanced diffusivity D_0 + D_1 gamma. The gel capacity
   *          F_cap is the densest fibrin measured by Wufsus et al. (2013,
   *          Biophys. J. 104:1812), ~100 mg/mL, in fibrinogen-equivalent
   *          units. The thrombomodulin velocity k_TM = k_on,TM Gamma_TM is the
   *          binding rate of thrombin to the thrombomodulin of one square
   *          metre of endothelium; with Da_TM = k_TM l_w/D_T ~ 1-10 the
   *          intact wall is close to a diffusion-limited sink, so the result
   *          depends weakly on its exact value.
   */
  struct PlateletKinetics
  {
      using Real = Rodin::Real;

      Real diffusivityProtein = 4.6e-11;    ///< D_P = D_T = D_A (m^2/s)
      Real diffusivityFibrinogen = 2.0e-11; ///< D_G (m^2/s)

      Real prothrombin0 = 1.4e-3;           ///< P_0 (mol/m^3)
      Real antithrombin0 = 2.8e-3;          ///< A_0 (mol/m^3)
      Real fibrinogen0 = 7.35e-3;           ///< G_0 (mol/m^3)
      Real prothrombinRef = 1.4e-3;         ///< P_ref (mol/m^3)

      Real amplificationRate = 0.10;        ///< k_a (1/s)
      Real amplificationHalf = 1.0e-5;      ///< K_A (mol/m^3)
      Real inhibitionRate = 7.1;            ///< k_on (m^3 mol^-1 s^-1)
      Real conversionRate = 42.0;           ///< k_F (1/s)
      Real conversionHalf = 3.6e-3;         ///< K_G (mol/m^3)

      Real wallFlux = 2.0e-12;              ///< j_0 at E = 1, P = P_ref (mol m^-2 s^-1)
      Real thrombomodulinRate = 1.0e-5;     ///< k_TM on the thrombomodulin set (m/s)

      /// @brief Fibrin is a gel (d_t F = R_F) or, legacy, F = G_0 - G.
      bool immobileFibrin = true;
      Real gelThreshold = 2.0e-3;           ///< K_F (mol/m^3), gel for Brinkman and S_s
      /// @brief Densest gel, fibrinogen-equivalent (mol/m^3); 0 disables the cap.
      Real fibrinCap = 0.3;                 ///< F_cap ~ 100 mg/mL ~ 40 G_0

      bool platelets = true;
      /// @brief Membrane availability that does not come from the modelled
      ///        platelets: 1 reproduces P-T-A-G (and a CAT with phospholipid),
      ///        0 is blood in which all procoagulant surface is platelet. It
      ///        describes the assay, not the blood.
      Real membraneBackground = 0.0;        ///< theta_0 (-)
      Real plateletCount = 2.5e14;          ///< N_0 (1/m^3), 250e9/L
      Real activationRate = 1.0;            ///< k_p (1/s)
      Real activationHalf = 3.0e-8;         ///< T_p (mol/m^3)
      Real plateletMax = 2.5e16;            ///< S_max (1/m^3)
      Real adhesionActivated = 4.0e-4;      ///< k_w (m/s), activated platelets
      Real adhesionResting = 4.0e-5;        ///< k_wN (m/s), resting platelets on collagen
      Real aggregationRate = 11.0;          ///< k_agg (1/s)
      Real aggregationLength = 2.78e-6;     ///< l, smoothing of S (m), one platelet diameter
      Real plateletDiffusivity0 = 1.58e-13; ///< D_0 (m^2/s)
      Real plateletDiffusivity1 = 6.0e-13;  ///< D_1 (m^2), times the shear rate

      Real stabilizationScale = 1.0;        ///< SUPG multiplier

      /// @brief Ablation switches (PTAGXBench_balance.pdf). Defaults: model as is.
      bool contactActivation = true;        ///< R_act stimulus includes S~/S_max
      bool linearFibrinogen = false;        ///< R_F = k_F T G/(K_G + G_0): no saturation in G
      Real membraneNorm = 0.0;              ///< N_ref normalising theta; 0 means N_0
      // ---- Surface amplification (PTAGXBench) -------------------------------
      /// @brief sigma (m): procoagulant membrane per unit area of the reactive
      ///        patch, spread over a layer of thickness delta_s next to it:
      ///        theta += (sigma/delta_s) chi w(x), w = exp(-d/delta_s),
      ///        chi such that int (sigma/delta_s) chi w dV = sigma |patch|.
      ///        Same constant as -ptag_sigma_s of LeftAtrium2DPTAGX.
      Real surfaceAmplification = 0.0;
      Real surfaceLayer = 10.0e-6;          ///< delta_s (m)

      /**
       * @brief Optional initiation enzyme X (factor Xa made by the TF wall).
       *
       * @details With it the TF wall emits X (Neumann load j_X E) instead of
       *          converting P into T directly; X converts prothrombin in the
       *          volume, R_X = kappa X P/P_ref, and is inactivated by
       *          antithrombin and TFPI with rate k_dX. kappa ~ 3 1/s is the
       *          turnover of a partly assembled prothrombinase (Xa with the Va
       *          available near the wall; full Xa.Va has k_cat ~ 30 1/s), k_dX
       *          ~ 1e-2 1/s the plasma lifetime of Xa (~100 s). Only the
       *          product kappa j_X matters for the thrombin made in quiescent
       *          plasma; j_X = k_cat,TF Gamma_TF (~5 1/s per TF-VIIa) is the
       *          TF density in flux units. At rest the layer of X has thickness
       *          sqrt(D/k_dX) ~ 70 um and builds up in ~1/k_dX, which gives the
       *          lag time and the fast initial phase of Thrombodynamics; under
       *          flow X is swept out of the Leveque layer before it acts, which
       *          is why TF patches under flow need far higher densities.
       */
      bool factorXa = true;
      Real xaFlux = 1.0e-13;                ///< j_X at E = 1 (mol m^-2 s^-1); Thrombodynamics-like activator
      Real xaRate = 3.0;                    ///< kappa (1/s), thrombin per Xa at P = P_0
      Real xaDecay = 1.0e-2;                ///< k_dX (1/s)

      /**
       * @brief Fraction phi of the feedback amplification routed through Xa
       *        (the intrinsic FXI -> FIXa -> FXa branch) instead of making
       *        thrombin at once.
       *
       * @details Thrombin's own feedback has two branches: it activates FV
       *          (prothrombinase, immediate) and FXI, whose FIXa:FVIIIa tenase
       *          makes more Xa after a delay. Here a fraction phi of the rate
       *          k_a theta T^2/(K_A + T) is emitted as Xa, (k_dX/kappa) times
       *          it, so that in a steady state kappa X reproduces exactly the
       *          same rate: phi redistributes the feedback in time (relaxation
       *          1/k_dX ~ 100 s) without changing its gain, and phi = 0 is the
       *          model as it stands. It stays quadratic in T at the origin, so
       *          R'(0) = 0 holds and a trace of thrombin still decays.
       */
      Real intrinsicFraction = 0.0;         ///< phi (-)

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


      /**
       * @brief Optional soluble fibrin monomer f between fibrinogen and gel.
       *
       * @details Thrombin makes monomer (R_F without the cap); the monomer
       *          polymerises bimolecularly, R_pol = k_pol f (f + F)[1-F/F_cap]_+
       *          (nucleation and growth on existing fibres, one constant,
       *          Weisel & Nagaswami 1992: ~1e6 1/(M s)). Monomer at nM is
       *          washed away before it gels; at uM it gels in seconds. This
       *          is the fibrin threshold of Colace et al. (2012).
       */
      bool fibrinMonomer = true;
      Real polymerizationRate = 1.0e3;      ///< k_pol (m^3 mol^-1 s^-1)

      struct State { Real P, T, A, G, F, N, M, S, X, Fm; };

      /// @brief Stabilised (gel-anchored) part of the deposit. Diagnostic.
      Real stabilized(Real S, Real F) const
      {
        const Real K = gelThreshold;
        return S * F * F / (K * K + F * F);
      }

      /// @brief Procoagulant membrane factor theta = theta_0 + (M + S)/N_0.
      Real membrane(Real M, Real S) const
      {
        if (!platelets)
          return membraneBackground;
        const Real norm = membraneNorm > 0.0 ? membraneNorm : plateletCount;
        return membraneBackground + (M + S) / norm;
      }

      /// @brief Activation stimulus in [0, 1]: thrombin, or contact with the
      ///        deposit within one platelet diameter (paracrine agonists).
      Real stimulus(Real T, Real smoothS) const
      {
        const Real byThrombin = T / (activationHalf + T);
        const Real byContact =
          contactActivation ? std::max<Real>(smoothS, 0.0) / plateletMax : 0.0;
        return std::min<Real>(1.0, byThrombin + byContact);
      }

      /// @brief Room left in the gel, [1 - F/F_cap]_+ (1 without a cap).
      Real gelRoom(Real F) const
      {
        if (fibrinCap <= 0.0)
          return 1.0;
        return std::max<Real>(0.0, 1.0 - F / fibrinCap);
      }

      /**
       * @brief One nodal reaction step (modified Patankar-Euler).
       *
       * @param s State after transport and wall deposition.
       * @param dt Time step.
       * @param smoothS Smoothed deposit S~ at the node.
       * @returns The new state: non-negative for any dt, with
       *          P + T - A, G + f + F (immobile fibrin) and N + M + S unchanged.
       */
      State advance(const State& s, Real dt, Real smoothS, Real surfaceTheta = 0.0) const
      {
        const Real P = std::max<Real>(s.P, 0.0), T = std::max<Real>(s.T, 0.0),
                   A = std::max<Real>(s.A, 0.0), G = std::max<Real>(s.G, 0.0),
                   F = std::max<Real>(s.F, 0.0), N = std::max<Real>(s.N, 0.0),
                   M = std::max<Real>(s.M, 0.0), S = std::max<Real>(s.S, 0.0),
                   X = std::max<Real>(s.X, 0.0), Fm = std::max<Real>(s.Fm, 0.0);

        const Real theta = membrane(M, S) + surfaceTheta;
        const Real ampT = amplificationRate * theta * T * T / (amplificationHalf + T);
        // Feedback split: a fraction phi goes through Xa (delayed), the rest
        // makes thrombin at once. Initiation enzyme: produced, then decays.
        const Real phi = factorXa ? std::clamp<Real>(intrinsicFraction, 0.0, 1.0) : 0.0;
        const Real toX = phi * ampT * xaDecay / xaRate;
        const Real Xn = factorXa ? (X + dt * toX) / (1.0 + dt * xaDecay) : 0.0;

        const Real thetaX = factorXa
          ? xaCofactorFloor + (1.0 - xaCofactorFloor) * T / (amplificationHalf + T) : 0.0;
        const Real amplification = (1.0 - phi) * ampT
                                 + (factorXa ? xaRate * Xn * thetaX : 0.0);
        const Real Pn = P / (1.0 + dt * amplification / prothrombinRef);
        const Real dP = P - Pn;

        const Real a = dt * inhibitionRate;
        const Real b = 1.0 + a * (A + T + dP);
        const Real c = a * A * (T + dP);
        const Real dI = 2.0 * c / (b + std::sqrt(std::max<Real>(b * b - 4.0 * a * c, 0.0)));

        const Real Tn = T + dP - dI;
        const Real An = A - dI;
        const bool monomer = immobileFibrin && fibrinMonomer;
        // Without the monomer the cap acts on the cleavage itself.
        const Real room = (immobileFibrin && !monomer) ? gelRoom(F) : 1.0;
        const Real Kg = conversionHalf + (linearFibrinogen ? fibrinogen0 : G);
        const Real Gn = G / (1.0 + dt * conversionRate * room * Tn / Kg);
        Real Fn, Fmn = 0.0;
        if (monomer)
        {
          const Real Fmt = Fm + (G - Gn);
          Fmn = Fmt / (1.0 + dt * polymerizationRate * (Fmt + F) * gelRoom(F));
          Fn = F + (Fmt - Fmn);
        }
        else
        {
          Fn = immobileFibrin ? F + (G - Gn)
                              : std::clamp<Real>(fibrinogen0 - Gn, 0.0, fibrinogen0);
        }

        State out{ Pn, Tn, An, Gn, Fn, N, M, S, Xn, Fmn };
        if (platelets)
        {
          const Real Nn = N / (1.0 + dt * activationRate * stimulus(Tn, smoothS));
          const Real Mt = M + (N - Nn);
          const Real r = dt * aggregationRate * std::max<Real>(smoothS, 0.0) / plateletMax;
          Real Sn = (S + r * Mt) / (1.0 + r * Mt / plateletMax);
          if (Sn - S > Mt)
            Sn = S + Mt;
          out.N = Nn;
          out.S = Sn;
          out.M = Mt - (Sn - S);
        }
        return out;
      }
  };

  /**
   * @brief Driver of the model above on a 2D or 3D mesh with a prescribed
   *        velocity.
   *
   * @details The velocity is a planar Poiseuille profile along x,
   *          u = gamma(t) y (1 - y/H) e_x, with wall shear rate gamma(t) at
   *          y = 0 (gamma = 0 gives a quiescent fluid). A different velocity
   *          only needs another implementation of updateVelocity().
   */
  class PTAGXBench
  {
    public:
      using Real = Rodin::Real;
      using Attribute = Rodin::Geometry::Attribute;
      using MeshType = Rodin::Geometry::Mesh<Rodin::Context::MPI>;

      using VectorFESType =
        Rodin::Variational::H1<1, Rodin::Math::SpatialVector<Real>, MeshType>;
      using ScalarFESType = Rodin::Variational::H1<1, Real, MeshType>;

      using VectorGridFunctionType =
        Rodin::PETSc::Variational::GridFunction<VectorFESType>;
      using ScalarGridFunctionType =
        Rodin::PETSc::Variational::GridFunction<ScalarFESType>;
      using ScalarTrialFunctionType =
        Rodin::PETSc::Variational::TrialFunction<ScalarGridFunctionType, ScalarFESType>;
      using ScalarTestFunctionType =
        Rodin::PETSc::Variational::TestFunction<ScalarFESType>;

      using LinearSystemType = Rodin::PETSc::Math::LinearSystem;
      using ScalarProblemType = Rodin::Variational::Problem<LinearSystemType,
        ScalarTrialFunctionType, ScalarTestFunctionType>;
      using FormType = Rodin::Variational::LinearForm<ScalarFESType, ::Vec>;

      struct Config
      {
          std::string meshPath = "channel.mesh";
          std::string xdmfBasename = "PTAGXBench";
          std::string csvPath = "PTAGXBench.csv";
          Real meshScale = 1.0;

          /// @brief Dirichlet (plasma composition) boundary.
          FlatSet<Attribute> inlet{ 1 };
          /// @brief Wall exposing tissue factor (lesion, activated endocardium).
          FlatSet<Attribute> reactive{ 5 };
        /// @brief x-range of the reactive patch on the wall y = 0, for the
        ///        surface layer (PTAGXBench; meshes of ../Test/make_mesh.py).
        std::array<Real, 2> patch{{-1.0e30, 1.0e30}};
        /// @brief Axisymmetric quiescent problem: x is the radius r (axis at
        ///        x = 0), every transport integral carries the weight r. A
        ///        patch [0, a] on y = 0 is then a disc of radius a.
        bool axisymmetric = false;
          /// @brief Wall exposing collagen (denuded subendothelium); empty on
          ///        an intact or merely activated endocardium.
          FlatSet<Attribute> collagen{};
          /// @brief Intact endothelium carrying thrombomodulin.
          FlatSet<Attribute> thrombomodulin{};
          /// @brief Activation E of the reactive wall.
          Real activation = 1.0;

          /// @brief Poiseuille channel: full height H and wall shear rate,
          ///        switched from gamma_1 to gamma_2 at t_switch (< 0: never).
          Real channelHeight = 1.0e-3;
          Real wallShearRate = 0.0;
          Real wallShearRate2 = 0.0;
          Real switchTime = -1.0;

          Real dt = 0.5;
          Real finalTime = 600.0;
          /// @brief CSV row every this many seconds; XDMF frame every outputInterval.
          Real csvInterval = 2.0;
          Real outputInterval = 60.0;
          int logEvery = 100;

          /// @brief Stop as soon as max T reaches this value (mol/m^3); 0: never.
          Real stopThrombin = 0.0;

          /// @brief Codina crosswind constant (0 disables it).
          Real crosswindC = 0.7;

          PlateletKinetics kinetics;
      };

      PTAGXBench(const Rodin::Context::MPI& context, const Config& cfg);

      PTAGXBench(const PTAGXBench&) = delete;
      PTAGXBench& operator=(const PTAGXBench&) = delete;

      PTAGXBench& initialize();
      int run();

    private:
      static MeshType makeMesh(const Rodin::Context::MPI& context, const Config& cfg);
      static Real cellSize(const Rodin::Geometry::Point& p);

      bool isRoot() const;

      Real wallShearRate(Real t) const;
      void updateVelocity(Real t);
      /// @brief Cellwise h_K, |u|^2 and shear rate at the centroid, refreshed
      ///        with the velocity; volume coefficients read them by cell index.
      void updateCellCache();
      Index cellOf(const Rodin::Geometry::Point& p) const;
      Real plateletDiffusivityAt(const Rodin::Geometry::Point& p) const;
      Real crosswind(const Rodin::Geometry::Point& p, Real cur, Real prev,
        const Rodin::Math::SpatialVector<Real>& gradient, Real diffusivity) const;

      /// @brief Is the boundary point on one of the attributes of the set?
      bool onWall(const Rodin::Geometry::Point& p, const FlatSet<Attribute>& set) const;
      /// @brief Adhesion velocity of activated platelets at a wall point:
      ///        k_w E on the reactive set, k_w on the collagen set, times the
      ///        free fraction (1 - S/S_max).
      Real adhesionActivatedAt(const Rodin::Geometry::Point& p) const;
      /// @brief Adhesion velocity of resting platelets: k_wN on collagen only.
      Real adhesionRestingAt(const Rodin::Geometry::Point& p) const;

      /**
       * @brief Backward-Euler SUPG + crosswind transport with the lagged
       *        velocity and the wall terms of one species.
       *
       * @param diffusivity Pointwise diffusivity D(x).
       * @param inlet Dirichlet value on the inlet set.
       * @param wall Attribute set carrying the Robin term and the load.
       * @param robin Implicit Robin coefficient on that set.
       * @param load Explicit wall load on that set.
       * @param sink Implicit Robin coefficient on the thrombomodulin set.
       */
      template <class DFn, class RFn, class LFn>
      void assignTransport(int species, ScalarProblemType& problem,
        ScalarTrialFunctionType& c, ScalarTestFunctionType& v,
        ScalarGridFunctionType& cur, ScalarGridFunctionType& prev,
        const DFn& diffusivity, Real inlet, const FlatSet<Attribute>& wall,
        const RFn& robin, const LFn& load, Real sink);

      /// @brief Cellwise Codina coefficient of one species at the centroid,
      ///        from the residual of the step just completed.
      void updateCrosswind(int species, const ScalarGridFunctionType& cur,
        const ScalarGridFunctionType& prev);

      void setupSpecies();
      void solveTransport(int species, ScalarProblemType& problem, Rodin::Solver::KSP& ksp,
        ScalarTrialFunctionType& trial, ScalarGridFunctionType& cur,
        ScalarGridFunctionType& prev);
      void depositAtWall();
      void smoothDeposit();
      void reactionStep();

      Real integrate(const ScalarGridFunctionType& f);
      void writeCSVHeader();
      void writeCSVRow();

      Config m_cfg;
      MeshType m_mesh;
      Rodin::IO::XDMF m_xdmf;

      VectorFESType m_vh;
      ScalarFESType m_sh;

      VectorGridFunctionType m_u;
      std::vector<Real> m_cellH, m_cellSpeed2, m_cellShear;
      /// @brief Per species (P, T, A, G, N, M, X, f): pointwise diffusivity
      ///        and cellwise crosswind coefficient.
      static constexpr int SpeciesCount = 8;
      std::array<std::function<Real(const Rodin::Geometry::Point&)>, SpeciesCount>
        m_speciesDiffusivity;
      std::array<std::vector<Real>, SpeciesCount> m_cellCrosswind;
      /// @brief Nodal states of one reaction step, in the order in which the
      ///        projection visits the DOFs (once each, increasing global index).
      std::vector<PlateletKinetics::State> m_stateCache;
      size_t m_cacheIndex = 0;
      /// @brief Walls on which platelets adhere: reactive and collagen sets.
      FlatSet<Attribute> m_adhesive;

      // Transported: P, T, A, G, N, M, and optionally X (Xa) and f (monomer).
      ScalarTrialFunctionType m_pT, m_tT, m_aT, m_gT, m_nT, m_mT, m_xT, m_fmT;
      ScalarTestFunctionType m_vP, m_vT, m_vA, m_vG, m_vN, m_vM, m_vX, m_vFm;
      ScalarGridFunctionType m_p, m_t, m_a, m_g, m_n, m_m, m_x, m_fm;
      ScalarGridFunctionType m_pPrev, m_tPrev, m_aPrev, m_gPrev, m_nPrev, m_mPrev,
        m_xPrev, m_fmPrev;
      // Immobile: F, S; smoothed deposit; next-state buffers.
      ScalarGridFunctionType m_f, m_s, m_sSmooth;
      ScalarGridFunctionType m_pNext, m_tNext, m_aNext, m_gNext, m_fNext, m_nNext,
        m_mNext, m_sNext, m_xNext, m_fmNext;
      // Wall deposition: nodal lumped mass and boundary loss vector.
      ScalarGridFunctionType m_lumped, m_deposit;
      // Diagnostics.
      ScalarGridFunctionType m_one, m_work, m_stab, m_theta, m_tat;
      /// @brief Surface layer w(x) and its normalisation (PTAGXBench).
      ScalarGridFunctionType m_surf;
      Real m_surfScale = 0.0;
      Real surfaceThetaAt(const Rodin::Geometry::Point& p) const;
      void setupSurfaceLayer();

      ScalarTestFunctionType m_q;
      FormType m_form;

      ScalarProblemType m_transportP, m_transportT, m_transportA, m_transportG,
        m_transportN, m_transportM, m_transportX, m_transportFm;
      Rodin::Solver::KSP m_kspP, m_kspT, m_kspA, m_kspG, m_kspN, m_kspM, m_kspX, m_kspFm;

      ScalarTrialFunctionType m_smoothTrial;
      ScalarTestFunctionType m_smoothTest;
      ScalarProblemType m_smooth;
      Rodin::Solver::KSP m_smoothKSP;

      Real m_tTransport = 0.0, m_tAssembly = 0.0, m_tDeposit = 0.0, m_tReaction = 0.0,
        m_tDiagnostics = 0.0; ///< wall-clock seconds per phase
      Real m_time = 0.0;
      Real m_gamma = -1.0;
      Real m_tmLoss = 0.0;       ///< int_0^t int_tm k_TM T, per unit depth (mol/m)
      Real m_clipped = 0.0;      ///< deposit returned to M at saturated wall nodes
      /// @brief Mass added by clipping transport undershoots to zero, per budget.
      Real m_clipPTA = 0.0, m_clipG = 0.0, m_clipNMS = 0.0;
      Real m_minTransport = 0.0; ///< most negative T after transport (mol/m^3)
      int m_step = 0;
      bool m_initialized = false;

      std::ofstream m_csv;
  };
}

#endif
