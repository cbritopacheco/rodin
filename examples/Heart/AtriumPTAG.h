/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file AtriumPTAG.h
 * @brief 3D atrium flow, 0D left ventricle and the P-T-A-G coagulation model.
 *
 * The flow is the one of Atrium.h: stabilised non-Newtonian P1/P1 flow with
 * four pressure inlets and one outlet closed by the 0D left ventricle through
 * a mitral diode. The thrombin/fibrinogen/fibrin kinetics of that example are
 * replaced by the precursor- and inhibitor-conserving reduced coagulation
 * model of LeftAtrium2DPTAG (the "P-T-A-G" model) and by a fibrin-dependent
 * Brinkman drag in the momentum equation.
 *
 * Species (mol/m^3, per plasma volume), transported by L c = d_t c + u.grad c
 * - div(D grad c):
 *
 *   L P = -R           prothrombin, the precursor
 *   L T =  R - R_I     thrombin, the enzyme
 *   L A =    - R_I     antithrombin, the stoichiometric inhibitor
 *   L G =      - R_F   fibrinogen, the substrate (FPA-cleavage capacity)
 *
 *   R   = k_a T^2/(K_A + T) P/P_ref     lumped feedback amplification
 *   R_I = k_on A T                       1:1 irreversible inhibition
 *   R_F = k_F T G/(K_G + G)              Michaelis-Menten conversion
 *
 * Initiation (endocardium): E in [0,1] is the endothelial activation built
 * from the cycle wall-shear indices (TAWSS, OSI), extended into an endocardial
 * reaction layer of thickness delta as E~; the P -> T channel
 * s = (j_0/delta) E~ P/P_ref carries the wall dose j_0 int E dA into the
 * layer. With delta = 0 it reduces to the wall flux D_T d_n T = j_0 E P/P_ref,
 * D_P d_n P = -j. Diagnostics: F = G_0 - G (pre-gel fibrin surrogate) and
 * I = A_0 - A (thrombin-antithrombin complex).
 *
 * Invariants of the continuous problem:
 *   P + T - A = P_0 - A_0,   P + T <= P_0,   A <= A_0,   0 <= G <= G_0.
 *
 * Flow feedback: the momentum equation carries -sigma_B(F) u with
 * sigma_B = (mu_inf/kappa_0) F^2/(K_F^2 + F^2) >= 0, assembled implicitly with
 * F lagged one step; the walls are rigid, so the term is exactly dissipative.
 *
 * Time integration: backward Euler for the flow and for the four transport
 * equations (velocity lagged, SUPG + Codina crosswind, optional Temam
 * correction, Robin term of P implicit, the T load using the P just solved),
 * followed by a nodal production-destruction step that is positive and
 * budget-exact for any dt (modified Patankar-Euler). Lie splitting.
 *
 * Equal order P1/P1, VMS (convection, grad-div) and PSPG. Runs under MPI; the
 * mesh is partitioned on the root rank and scaled by Config::meshScale.
 */
#ifndef EXAMPLES_HEART_ATRIUMPTAG_H
#define EXAMPLES_HEART_ATRIUMPTAG_H

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include "Rodin/Heart/CCMLC2014.h"
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
  /**
   * @brief Constants of the P-T-A-G kinetics (SI: mol, m, s).
   *
   * @details Same definition as in LeftAtrium2DPTAG.h (branch
   *          test/forcedconvection); kept here so the 3D case is
   *          self-contained until the two examples share a header.
   *          Literature values unless marked. 1 uM = 1e-3 mol/m^3.
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
       *          0.5-1 mm), so the activated surface is a layer, not a
       *          surface. The activation E is extended into the lumen by
       *            -delta^2 Laplace(E~) + E~ = 0,  E~ = E on the wall,
       *          i.e. E~ ~ E exp(-d/delta), and the initiation becomes the
       *          volumetric P -> T channel s = (j_0/delta) E~ P/P_ref, scaled
       *          so that int s dV equals the wall dose j_0 int E dA exactly.
       *          With a Robin flux on no-slip nodes the products can leave
       *          the wall only by molecular diffusion; in the layer they are
       *          born where the fluid moves.
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

  class AtriumPTAG
  {
    public:
      using Real = Rodin::Real;
      using Model = Rodin::Heart::CCMLC2014T<>;
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

      using VectorTrialFunctionType =
        Rodin::PETSc::Variational::TrialFunction<VectorGridFunctionType, VectorFESType>;
      using ScalarTrialFunctionType =
        Rodin::PETSc::Variational::TrialFunction<ScalarGridFunctionType, ScalarFESType>;
      using VectorTestFunctionType =
        Rodin::PETSc::Variational::TestFunction<VectorFESType>;
      using ScalarTestFunctionType =
        Rodin::PETSc::Variational::TestFunction<ScalarFESType>;

      using LinearSystemType = Rodin::PETSc::Math::LinearSystem;

      using FlowProblemType =
        Rodin::Variational::Problem<LinearSystemType, VectorTrialFunctionType,
          ScalarTrialFunctionType, VectorTestFunctionType, ScalarTestFunctionType>;
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

      /// @brief Carreau-Yasuda blood viscosity.
      struct CarreauYasuda
      {
          Real mu0 = 0.0720;
          Real muInf = 0.0049;
          Real lambda = 6.1529;
          Real n = 0.2081;
          Real yasuda = 1.4173;
          Real gammaRegularization = 1.0e-3;
      };

      /// @brief Piecewise-linear LV activation waveform (as in Atrium.h).
      struct Activation
      {
          Real period = 0.85;
          Real tRampStart = 0.15;
          Real tRampEnd = 0.21;
          Real tPlateauEnd = 0.36;
          Real tRelaxEnd = 0.45;
          Real tNegativeEnd = 0.6;
          Real positiveValue = 35.0;
          Real negativeValue = -20.0;
      };

      /// @brief Piecewise-linear atrial pressure waveform driving the 0D LV.
      struct AtrialPressure
      {
          Real period = 0.85;
          Real minValue = 500.0;
          Real maxValue = 1000.0;
          Real secondThreshold = 1250.0;
          Real t1 = 0.02;
          Real t2 = 0.15;
          Real t3 = 0.17;
          Real t4 = 0.56;
          Real t5 = 0.62;
          Real t6 = 0.85;
      };

      /// @brief 0D left ventricle parameters and initial conditions (as in Atrium.h).
      struct LVModel
      {
          Real rho = 1.0e3;
          Real R0 = 2.45e-2;
          Real d0 = 1.4e-2;
          Real Es = 3.0e6;
          Real mu = 70.0;
          Real eta = 70.0;
          Real alpha = 1.5;
          Real alphaR = 0.12;
          Real k0 = 1.0e5;
          Real sigma0 = 1.2e5;
          Real Rp = 5.0e7;
          Real Cp = 5.0e-9;
          Real Rd = 1.0e8;
          Real Cd = 5.0e-10;
          Real Kat = 1.0e-6;
          Real Kp = 5.0e-10;
          Real Kar = 2.0e-7;
          Real cavityCapacity = 5.0e-12;
          Real localTolerance = 1.0e-12;
          int localMaxIterations = 50;
          Real localDamping = 1.0;
          Real absRegularization = 1.0e-14;
          Real passiveC0 = 1.9e3;
          Real passiveC1 = 1.1e-1;
          Real passiveC2 = 1.9e3;
          Real passiveC3 = 1.1e-1;
          Real relaxationM0Low = 1.6;
          Real relaxationM0High = 1.0;
          Real relaxationM0LowEc = 0.0;
          Real relaxationM0HighEc = 2.0;
          Real systemicVenousPressure = 1.0e3;
          int maxIterations = 200;
          Real absoluteTolerance = 1.0e-8;
          Real relativeTolerance = 1.0e-8;
          Real stepTolerance = 1.0e-10;
          Real dampingFactor = 1.0;
          Real initialY = 0.0;
          Real initialV = 0.0;
          Real initialPvOffset = -100.0;
          Real initialPar = 11000.0;
          Real initialPd = 10000.0;
          Real mu_0 = 0.301;
          Real mu_Inf = 0.0055;
          Real lambda = 16.152;
          Real n = 0.21;
          Real yasuda = 0.77;
          Real proximalRadius = 0.0125;
          Real proximalLength = 0.4;
          Real distalRadius = 0.00175;
          Real distalLength = 0.2;
      };

      struct Config
      {
          /// @brief Same mesh, labels and scale as Atrium.h.
          std::string meshPath = "../resources/examples/Heart/Atrium/Coarse.mesh";
          std::string xdmfBasename = "AtriumPTAG";
          std::string csvPath = "AtriumPTAG.csv";

          /// @brief No-slip wall (endocardium).
          Attribute wall = 4;
          /// @brief Pulmonary vein inlets.
          std::array<Attribute, 4> inlets{{1, 2, 3, 5}};
          /// @brief Mitral outlet.
          Attribute outlet = 6;

          /// @brief Mesh coordinate scale applied after partitioning (mm -> m).
          Real meshScale = 1.0e-3;

          Real rho = 1060.0;
          CarreauYasuda viscosity;
          /// @brief Pressure penalty of the equal-order pair.
          Real pressurePenalty = 1.0e-12;

          /// @brief Pulmonary venous pressure, p = mean + amplitude sin(2 pi t/T).
          Real inletPressureMean = 1200.0;
          Real inletPressureAmplitude = 200.0;

          /// @brief Normal impedance at the pressure-driven inlets (Pa s/m);
          ///        see Atrium.h for its role in the energy bound.
          Real inletImpedance = 1.0e3;
          /// @brief Tangential damping at the inlets (Pa s/m).
          Real inletTangentialDamping = 1.0e3;
          Real inletBackflowStabilization = 1.0;
          Real outletBackflowStabilization = 1.0;

          /// @brief Mitral resistances (Pa s/m^3). The valve is a diode: open
          ///        while the atrium is above the ventricle, closed otherwise.
          Real mitralResistanceOpen = 5.0e6;
          Real mitralResistanceClosed = 1.0e12;

          Real dt = 1.0e-3;
          Real period = 0.85;
          /// @brief Cycles solved before the species are switched on (the
          ///        indices need one closed cycle, so this is at least 1).
          int flowCycles = 2;
          int speciesCycles = 10;

          /// @brief XDMF output period, in steps. Every cycle boundary is
          ///        written regardless, so 0 keeps only those.
          int outputEvery = 10;

          Real vmsScale = 1.0;
          Real gradDivScale = 1.0;
          Real pspgScale = 1.0;
          bool useVMS = true;
          /// @brief Quasi-static (OSS) instead of dynamic velocity subscales.
          bool quasiStaticSubscales = false;
          /// @brief Weight of the transient, convective and drag parts of the
          ///        PSPG residual: 1 consistent, 0 pressure-gradient penalty.
          /// @details Atrium.h runs the pressure-gradient penalty only, so 0
          ///          reproduces that flow exactly (while sigma_B = 0).
          Real pspgResidualScale = 0.0;
          /// @brief Codina crosswind constant of the transport. 0 disables it.
          Real crosswindC = 0.7;

          /// @brief Scale of the Temam correction (1/2)(div u) c v of the four
          ///        transport equations; see Atrium.h for the derivation.
          /// @details Off by default: this is the form LeftAtrium2DPTAG was
          ///          calibrated with. With the equal-order pair the residual
          ///          divergence is a cellwise reaction coefficient of size
          ///          |u|/h, and where it is negative and dt |div u|/2 > 1
          ///          the implicit step amplifies instead of damping; on a
          ///          coarse 3D mesh that drove P, A and G to zero within a
          ///          few steps. Linear in c, so when enabled it leaves the
          ///          invariant P + T - A untouched.
          Real speciesTemamScale = 0.0;

          /// @brief GMRES + block-Jacobi/ILU for the four transport solves
          ///        instead of the global direct solver.
          bool iterativeTransport = false;

          bool solveKinetics = true;
          ConservingKinetics kinetics;
          BrinkmanDrag brinkman;

          /// @brief Cycle replay. After the warm-up the last flow cycle is
          ///        stored and replayed for the species; the flow (and the 0D
          ///        ventricle) is solved again only when the fibrin has
          ///        changed enough to alter it.
          /// @details Exact while sigma_B = 0; afterwards the error is the
          ///          one-cycle lag of a field that evolves over tens of
          ///          cycles. Memory: one velocity Vec per step of the cycle.
          bool replay = false;
          /// @brief Flow cycles solved at each refresh (the last one is stored).
          int replayFlowCycles = 2;
          /// @brief Refresh when max F > K_F/2 and the fibrin mass has grown by
          ///        this fraction since the stored cycle.
          Real replayRefresh = 0.2;
          /// @brief Advance the species every k flow steps with dt_q = k dt.
          ///        Must divide the steps per cycle.
          int speciesSubsteps = 1;

          /// @brief Divergence guard on max|u| (m/s).
          Real maxVelocity = 5.0;

          Activation activation;
          AtrialPressure atrialPressure;
          LVModel lv;
      };

      /// @brief Wall-clock cost of each phase of one step.
      struct Timing
      {
          Real zeroD = 0.0;
          Real vms = 0.0;
          Real assembly = 0.0;
          Real solve = 0.0;
          Real shear = 0.0;
          Real fluxes = 0.0;
          Real species = 0.0;
          Real output = 0.0;
          Real total = 0.0;
      };

      AtriumPTAG(const Rodin::Context::MPI& context, const Config& cfg);
      ~AtriumPTAG();

      AtriumPTAG(const AtriumPTAG&) = delete;
      AtriumPTAG& operator=(const AtriumPTAG&) = delete;

      AtriumPTAG& initialize();
      int run();

      Config& getConfig() noexcept
      {
        return m_cfg;
      }
      const Config& getConfig() const noexcept
      {
        return m_cfg;
      }

    private:
      static MeshType makeMesh(const Rodin::Context::MPI& context, const Config& cfg);
      static Model::Input makeInput(const Config& cfg);
      static Real activationWave(const Activation& cfg, Real t);
      static Real atrialWave(const AtrialPressure& cfg, Real t);
      static Real cellSize(const Rodin::Geometry::Point& p);
      /// @brief y <- y + a x, with the ghost values refreshed.
      static void axpy(Real a, const ::Vec& x, ::Vec& y);

      bool isRoot() const;

      void setupSpaces();
      void setupFlow();
      void setupSpecies();
      void setupWallShear();
      void setupLayer();
      /// @brief Whether the point is a wall node (see m_wallArea).
      bool onWallNode(const Rodin::Geometry::Point& p) const;
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

      /// @brief L2 projection of a scalar expression, reusing one mass matrix.
      template <class Expression>
      void project(const Expression& expr, ScalarGridFunctionType& out)
      {
        m_scalarProjection = Rodin::Variational::Integral(m_sTrial, m_sTest) -
          Rodin::Variational::Integral(expr, m_sTest);
        m_scalarProjection.solve(m_scalarProjectionKSP);
        out.setData(m_sTrial.getSolution().getData());
      }

      /// @brief L2 projection of a vector expression, reusing one mass matrix.
      template <class Expression>
      void projectVector(const Expression& expr, VectorGridFunctionType& out)
      {
        m_vectorProjection = Rodin::Variational::Integral(m_wTrial, m_wTest) -
          Rodin::Variational::Integral(expr, m_wTest);
        m_vectorProjection.solve(m_vectorProjectionKSP);
        out.setData(m_wTrial.getSolution().getData());
      }

      bool advance0D();
      bool solveFlow();
      /// @brief Copy u^{n+1} into slot k of the stored cycle.
      void storeSnapshot(int k);
      /// @brief Load slot k as the lagged velocity.
      void loadSnapshot(int k);
      /// @brief int_Omega F dV, the refresh monitor.
      Real fibrinMass();
      void solveSpecies();
      void reactionStep();
      void updateDiagnostics();
      void computeWallShear();
      void closeCycle(Real elapsed);
      void computeFluxes();
      Real boundaryArea(Attribute tag);

      /// @brief Local Carreau-Yasuda viscosity at the lagged velocity.
      Real viscosityAt(const Rodin::Geometry::Point& p) const;
      /// @brief Fibrin surrogate F = G_0 - G at a point, clamped to [0, G_0].
      Real fibrinAt(const Rodin::Geometry::Point& p) const;
      /// @brief sigma_B(F) at a point, from the lagged fibrinogen.
      Real dragAt(const Rodin::Geometry::Point& p) const;
      /// @brief Codina tau_1 with the Brinkman reaction, on the local viscosity.
      Real tau1At(const Rodin::Geometry::Point& p) const;
      Real vmsTauAt(const Rodin::Geometry::Point& p) const;
      Real sqrtTauCAt(const Rodin::Geometry::Point& p) const;
      Real tauCAt(const Rodin::Geometry::Point& p) const;
      Real tauPAt(const Rodin::Geometry::Point& p) const;
      Real crosswind(const Rodin::Geometry::Point& p, Real cur, Real prev,
        const Rodin::Math::SpatialVector<Real>& gradient, Real diffusivity) const;

      void writeCSVHeader();
      void writeCSVRow(int cycle);

      Config m_cfg;
      Model::Input m_input;
      Model m_model;

      MeshType m_mesh;
      Rodin::IO::XDMF m_xdmf;
      AttributeSet m_inletSet;

      VectorFESType m_vh;
      ScalarFESType m_sh;

      // ---- Flow ------------------------------------------------------------
      VectorTrialFunctionType m_u;
      ScalarTrialFunctionType m_p;
      VectorTestFunctionType m_v;
      ScalarTestFunctionType m_q;
      VectorGridFunctionType m_uOld;

      // ---- Reusable mass projections ---------------------------------------
      ScalarTrialFunctionType m_sTrial;
      ScalarTestFunctionType m_sTest;
      VectorTrialFunctionType m_wTrial;
      VectorTestFunctionType m_wTest;

      // ---- VMS, PSPG and drag coefficients (pointwise) ----------------------
      ScalarCoefficientType m_tauFn;
      ScalarCoefficientType m_tauCFn;
      ScalarCoefficientType m_sqrtTauCFn;
      ScalarCoefficientType m_tauPFn;
      ScalarCoefficientType m_dragFn;
      ScalarCoefficientType m_pspgDragFn;
      ScalarGridFunctionType m_piTilde;
      VectorGridFunctionType m_convProjection;
      VectorGridFunctionType m_sub;
      VectorGridFunctionType m_subOld;

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
      /// @brief int_Gw phi_i dA: positive exactly at the wall nodes. Nodal
      ///        indicator used to restrict the cycle indices to the wall.
      /// @details The wall quantities are projected over the cells, never
      ///          over Region::Boundary: that iterator visits only the owned
      ///          boundary faces, so a rank that holds a wall node as a ghost
      ///          without owning a wall face at it never writes that slot,
      ///          and the reverse scatter of GridFunction::flush() (ghost
      ///          wins) then overwrites the owner's value with the stale
      ///          one. The result was a dashed line of zeros along every
      ///          partition seam in TAWSS, OSI and E.
      ScalarGridFunctionType m_wallArea;
      Real m_wallAreaThreshold = 0.0;
      /// @brief Ghost-consistent copy of the recovered traction.
      VectorGridFunctionType m_wssRec;
      VectorGridFunctionType m_wss;
      VectorGridFunctionType m_gradRec0;
      VectorGridFunctionType m_gradRec1;
      VectorGridFunctionType m_gradRec2;
      VectorGridFunctionType m_netShear;
      ScalarGridFunctionType m_absShear;
      ScalarGridFunctionType m_shearMagnitude;
      ScalarGridFunctionType m_tawss;
      ScalarGridFunctionType m_osi;
      ScalarGridFunctionType m_activation;

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

      /// @brief One stored cycle of u, one ghosted Vec per step (replay).
      std::vector<::Vec> m_cycleU;

      // ---- Time-dependent coefficients read by the frozen forms ------------
      /// @brief Species time step, speciesSubsteps * dt.
      Real m_dtChem = 0.0;
      Real m_t = 0.0;
      Real m_pIn = 0.0;
      Real m_pOut = 0.0;
      Real m_mitralZ = 0.0;
      Real m_outletArea = 1.0;
      Real m_outletPressure = 0.0;
      Real m_qIn = 0.0;
      Real m_qOut = 0.0;
      Real m_speed = 0.0;
      bool m_valveOpen = true;
      bool m_flowFieldSplitsSet = false;
      bool m_initialized = false;
      bool m_indicesReady = false;

      /// @brief Index of the step being solved; the first ones are traced
      ///        phase by phase, so a stall is visible where it happens.
      int m_step = 0;
      Timing m_timing;

      std::ofstream m_csv;
  };
}

#endif
