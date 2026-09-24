// LeftAtrium2DPTAG.h
//
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
#ifndef EXAMPLES_HEART_LEFTATRIUM2DPTAG_H
#define EXAMPLES_HEART_LEFTATRIUM2DPTAG_H

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

      struct State { Real P, T, A, G; };

      /**
       * @brief One nodal reaction step, production explicit / destruction
       *        implicit (modified Patankar-Euler).
       *
       * @param wallRate The layer initiation rate s/P (1/s) at this node.
       * @details P* = P/(1 + dt [k_a T^2/((K_A+T) P_ref) + wallRate]),  dP = P - P*;
       *          dI = smaller root of dI = dt k_on (A - dI)(T + dP - dI);
       *          T* = T + dP - dI,  A* = A - dI,
       *          G* = G/(1 + dt k_F T/(K_G + G)).
       *          For any dt: all values stay >= 0, P* + T* - A* = P + T - A
       *          exactly, and the prothrombin consumed equals the thrombin
       *          produced.
       */
      State advance(const State& s, Real dt, Real wallRate = 0.0) const
      {
        // The transport step (SUPG) does not guarantee positivity; an
        // undershoot is removed here, before it can react. The price is a
        // drift of P + T - A that the CSV reports, so it is never silent.
        const Real P = std::max<Real>(s.P, 0.0), T = std::max<Real>(s.T, 0.0),
                   A = std::max<Real>(s.A, 0.0), G = std::max<Real>(s.G, 0.0);

        const Real Pn = P / (1.0 + dt * (amplificationRate * T * T /
                                     ((amplificationHalf + T) * prothrombinRef) + wallRate));
        const Real dP = P - Pn;

        const Real a = dt * inhibitionRate;
        const Real b = 1.0 + a * (A + T + dP);
        const Real c = a * A * (T + dP);
        const Real dI = 2.0 * c / (b + std::sqrt(std::max<Real>(b * b - 4.0 * a * c, 0.0)));

        const Real Tn = T + dP - dI;
        const Real An = A - dI;
        const Real Gn = G / (1.0 + dt * conversionRate * Tn / (conversionHalf + G));
        return { Pn, Tn, An, Gn };
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

  class LeftAtrium2DPTAG
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
          std::string xdmfBasename = "LeftAtrium2DPTAG";
          std::string csvPath = "LeftAtrium2DPTAG.csv";

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

      LeftAtrium2DPTAG(const Rodin::Context::MPI& context, const Config& cfg);
      ~LeftAtrium2DPTAG();

      LeftAtrium2DPTAG(const LeftAtrium2DPTAG&) = delete;
      LeftAtrium2DPTAG& operator=(const LeftAtrium2DPTAG&) = delete;

      LeftAtrium2DPTAG& initialize();
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
        Real robinScale, const ScalarGridFunctionType* load);

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
      /// @brief Fibrin surrogate F = G_0 - G at a point, clamped to [0, G_0].
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

      // ---- Endocardial layer ------------------------------------------------
      ScalarTrialFunctionType m_eT;
      ScalarTestFunctionType m_vE;
      /// @brief E~, the activation extended into the lumen; zero if the layer
      ///        thickness is zero.
      ScalarGridFunctionType m_layer;
      /// @brief Dose normalisation, delta int_Gw E dA / int_Omega E~ dV.
      Real m_layerScale = 1.0;

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
      ScalarProblemType m_transportP, m_transportT, m_transportA, m_transportG;
      Rodin::Solver::KSP m_kspP, m_kspT, m_kspA, m_kspG;
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
