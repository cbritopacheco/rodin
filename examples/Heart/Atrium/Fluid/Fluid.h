/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Fluid.h
 * @brief Distributed ALE Oseen P1/P1 atrial flow with split-VMS or SUPG/PSPG stabilization.
 *
 * Equal-order P1/P1 Oseen flow uses either term-by-term orthogonal subscales
 * (convection along u^n - w, grad-div, and pressure-gradient stabilization) or
 * residual-based SUPG/PSPG with grad-div. Standalone flow has pressure-driven
 * openings and a fixed no-slip wall. FSI supplies Robin interface data and
 * harmonic ALE motion. The lumped ventricle with its mitral diode, the
 * orthogonal subscales, and the P-T-A-G coagulation transport are available in
 * both fluid and FSI modes.
 * The fluid submesh retains the common MPI partition and parent identities;
 * its coordinates are converted to SI units with Parameters::meshScale.
 */
#ifndef RODIN_ATRIUM_FLUID_FLUID_H
#define RODIN_ATRIUM_FLUID_FLUID_H

#include <array>
#include <chrono>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "Rodin/Heart/CCMLC2014.h"
#include <Rodin/Geometry.h>
#include <Rodin/MPI.h>
#include <Rodin/PETSc.h>
#include <Rodin/Solver.h>
#include <Rodin/Types.h>
#include <Rodin/Variational.h>

#include "../../CoronaryArtery/ThrombosisModel.h"

#include "../RunParameters.h"
#include "Coagulation.h"

namespace Rodin::Examples::Heart::Atrium
{
  class FSI;
  class Fluid
  {
    public:
      class Output;

      using Real = Rodin::Real;
      using Model = Rodin::Heart::CCMLC2014T<>;
      using Attribute = Rodin::Geometry::Attribute;
      using ParentMeshType = Rodin::Geometry::Mesh<Rodin::Context::MPI>;
      using MeshType = Rodin::Geometry::SubMesh<Rodin::Context::MPI>;

      using VectorFESType =
        Rodin::Variational::H1<1, Rodin::Math::SpatialVector<Real>, ParentMeshType>;
      using ScalarFESType = Rodin::Variational::H1<1, Real, ParentMeshType>;

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

      /**
       * @brief A scalar coefficient evaluated pointwise.
       * @details The stabilisation parameters are element-local quantities and
       * are evaluated where they are used, never L2-projected onto a
       * continuous field: the projection of a positive function that
       * varies by a decade between neighbouring cells overshoots, and
       * a tau that comes out negative at a node turns a
       * positive-semidefinite stabilisation into an indefinite one.
       * It is also five mass solves per step cheaper.
       */
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

      /// @brief Piecewise-linear LV activation waveform.
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

      /// @brief 0D left ventricle parameters and initial conditions.
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

      /// @brief Scalar diagnostics from the most recent completed flow state, in SI units.
      struct Diagnostics
      {
        int cycle;
        Real inletPressure, outletPressure;
        Optional<Real> arterialPressure, distalPressure;
        bool valveOpen;
        Real inletFlux, outletFlux, maximumVelocity;
        Real maximumTAWSS, maximumOSI, maximumActivation;
        Real maximumThrombin, minimumFibrinogen, maximumFibrin;
        int linearIterations, preconditionerAge, preconditionerRefreshes;
        bool preconditionerRefreshed, preconditionerRetried;
      };

      struct Parameters
      {
        /// @brief No-slip wall.
        Attribute wall = 4;

        /// @brief Pulmonary vein inlets.
        std::array<Attribute, 4> inlets{{1, 2, 3, 5}};

        /// @brief Mitral outlet.
        Attribute outlet = 6;

        /// @brief Mesh coordinate scale applied after partitioning.
        Real meshScale = 1.0e-3;

        Real rho = 1060.0;
        CarreauYasuda viscosity;

        /// @brief Pressure penalty of the equal-order pair.
        Real pressurePenalty = 1.0e-12;

        /// @brief Pulmonary venous pressure, p = mean + amplitude sin(2 pi t/T).
        Real inletPressureMean = 1200.0;
        Real inletPressureAmplitude = 200.0;

        /**
         * @brief Normal impedance at the pressure-driven inlets (Pa s/m).
         * @details A pressure inlet with zero source impedance has no energy
         * bound: the convective boundary flux grows like |u|^3 while
         * only viscous dissipation opposes it.
         */
        Real inletImpedance = 1.0e3;

        /// @brief Tangential damping at the inlets (Pa s/m).
        Real inletTangentialDamping = 1.0e3;
        Real inletBackflowStabilization = 1.0;
        Real outletBackflowStabilization = 1.0;

        /**
         * @brief Mitral resistances (Pa s/m^3). The valve is a diode: open
         * while the atrium is above the ventricle, closed otherwise.
         */
        Real mitralResistanceOpen = 5.0e6;
        Real mitralResistanceClosed = 1.0e12;

        /// @brief Advance the lumped ventricle and use its outlet pressure.
        bool useVentricularModel = false;

        /// @brief Prescribed outlet pressure when the ventricle is disabled, in Pa.
        Real prescribedOutletPressure = 0.0;

        /// @brief Outlet resistance when the ventricle is disabled, in Pa s/m^3.
        Real prescribedOutletResistance = 0.0;

        /// @brief Numeric preconditioner refresh interval in timesteps; one refreshes every step.
        int preconditionerLag = 1;

        Real period = 0.85;

        /**
         * @brief Linear ramp time of the opening pressures, in seconds.
         * @details Zero applies the full venous and ventricular pressures from
         * the first step. A positive value avoids an impulsive load on a wall
         * whose reference configuration is unloaded.
         */
        Real pressureRamp = 0.0;

        /// @brief Flow cycles completed before the species are switched on.
        int flowCycles = 3;

        /**
         * @brief Per-term stabilisation scales, so each can be dialled down
         * on its own. 1 is the standard value, 0 removes the term.
         */
        Real vmsScale = 1.0;
        Real gradDivScale = 1.0;
        Real pspgScale = 1.0;

        /**
         * @brief Selects term-by-term orthogonal-subscale stabilization.
         * @details Convection, divergence, and pressure gradient are each
         * stabilized separately. The convective and divergence projections use
         * the previous time level. False selects residual-based SUPG/PSPG.
         */
        bool useVMS = false;

        /// @brief Quasi-static instead of dynamic convective subscales.
        bool quasiStaticSubscales = false;

        /// @brief Constants of the P-T-A-G coagulation model.
        ConservingKinetics kinetics;

        /// @brief Fibrin-dependent drag in the momentum equation.
        BrinkmanDrag brinkman;

        /// @brief Codina crosswind constant. 0 disables it.
        Real crosswindC = 0.7;

        /// @brief Transport and react the P-T-A-G species after the flow cycles.
        bool solveKinetics = false;

        /**
         * @brief Numerical divergence guard on maximum speed (m/s).
         * @details The default is three orders above a 1 m/s reference speed;
         * this heuristic is not a physiological acceptance criterion.
         */
        Real maxVelocity = 1e3;

        Activation activation;
        AtrialPressure atrialPressure;
        LVModel lv;
      };

      /// @brief Backward-compatible name for the atrium parameters.

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

      Fluid(const ParentMeshType& parent, const Parameters& cfg, const RunParameters& run,
        Attribute material = 1);
      ~Fluid();

      Fluid(const Fluid&) = delete;
      Fluid& operator=(const Fluid&) = delete;

      /**
       * @brief Initializes retained mesh-dependent fields, forms, solvers, and output.
       * @returns This simulation.
       * @note Collective on the atrium communicator. Configuration remains fixed
       *   after construction because forms retain their coefficients.
       */
      Fluid& initialize();

      /**
       * @brief Advances one fluid, 0D, and enabled transport time step.
       * @note Collective on the atrium communicator. A failed step cannot be
       *   retried on this instance. Call finish() after the last successful step.
       */
      void step();

      /**
       * @brief Finalizes output, writing the last completed state if necessary.
       * @note Collective on the atrium communicator. No steps may follow.
       */
      void finish();

      /// @brief Current physical time in seconds.
      Real getTime() const noexcept
      {
        return m_t;
      }

      /// @brief Number of successfully completed steps.
      int getStepCount() const noexcept
      {
        return m_step;
      }

      /**
       * @brief Returns the owned distributed reference mesh.
       * @returns Mesh in SI length units.
       */
      const MeshType& getMesh() const noexcept
      {
        return m_mesh;
      }

      /**
       * @brief Returns the fluid velocity field.
       * @returns Observed velocity in m/s; meaningful after a successful step.
       */
      const VectorGridFunctionType& getVelocity() const noexcept
      {
        return m_u.getSolution();
      }

      /**
       * @brief Returns the fluid pressure field.
       * @returns Observed pressure in Pa; meaningful after a successful step.
       */
      const ScalarGridFunctionType& getPressure() const noexcept
      {
        return m_p.getSolution();
      }

      /// @brief Returns cached flow diagnostics, or empty before the first completed step.
      const Optional<Diagnostics>& getDiagnostics() const noexcept
      {
        return m_diagnostics;
      }

      /// @brief Returns the parameters supplied to the simulation.
      const Parameters& getParameters() const noexcept
      {
        return m_cfg;
      }

    private:
      friend class FSI;

      /// @brief Freezes history and prepares the next physical time.
      void beginStep();

      /// @brief Commits the converged trial, cycle statistics, and output.
      void commitStep();

      /// @brief Adds the accepted wall shear to the cycle statistics.
      void accumulateWallShear();

      static MeshType makeMesh(
        const ParentMeshType& parent, const Parameters& cfg, Attribute material);
      static Model::Input makeInput(const Parameters& cfg);
      static Real activationWave(const Activation& cfg, Real t);
      static Real atrialWave(const AtrialPressure& cfg, Real t);

      bool isRoot() const;

      void setupSpaces();
      void setupFlow();
      void setupSpecies();
      void setupWallShear();

      /// @brief Refreshes the lagged orthogonal-subscale projections on the previous configuration.
      void prepareSubscales();

      /**
       * @brief Assigns the backward-Euler ALE transport of one species.
       * @param problem Retained transport problem.
       * @param c Trial concentration.
       * @param v Test function.
       * @param cur Current nodal concentration, moving with the mesh.
       * @param prev Previous nodal concentration for the crosswind residual.
       * @param diffusivity Molecular diffusivity (m^2/s).
       * @param inlet Concentration prescribed at the venous openings.
       * @param robinScale Wall exchange coefficient; zero without wall exchange.
       * @param load Field supplying the wall load, or null.
       */
      void assignTransport(ScalarProblemType& problem, ScalarTrialFunctionType& c,
        ScalarTestFunctionType& v, ScalarGridFunctionType& cur,
        ScalarGridFunctionType& prev, Real diffusivity, Real inlet, Real robinScale,
        const ScalarGridFunctionType* load);

      /// @brief Extends the wall activation into the endocardial layer and fixes its dose.
      void extendActivation();

      /// @brief Nodal production-destruction step of the four species.
      void reactionStep();

      /// @brief Updates fibrin, thrombin-antithrombin, drag, and invariant-drift fields.
      void updateCoagulationDiagnostics();

      /// @brief Returns whether the evaluation point is a wall node.
      bool onWallNode(const Rodin::Geometry::Point& p) const;

      /// @brief Fibrin surrogate F = G_0 - G at a point, clamped to [0, G_0].
      Real fibrinAt(const Rodin::Geometry::Point& p) const;

      /// @brief Brinkman coefficient from the lagged fibrinogen; zero before the species start.
      Real dragAt(const Rodin::Geometry::Point& p) const;

      /// @brief Cellwise viscosity and size shared by all stabilization coefficients.
      struct CellCoefficients
      {
          size_t epoch = 0;
          Real viscosity = 0, size = 0;
      };

      /**
       * @brief Returns the lagged viscosity and size of the cell containing a point.
       * @details Both are constant on an affine P1 tetrahedron. Entries are
       * reused until the lagged velocity or the geometry changes.
       */
      const CellCoefficients& cellCoefficients(const Rodin::Geometry::Point& p) const;

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

      /// @brief y <- y + a x, with the ghost values refreshed.
      static void axpy(Real a, const ::Vec& x, ::Vec& y);

      bool advance0D();
      bool solveFlow();
      void solveSpecies();

      /// @brief Reassembles and solves the three retained gradient-recovery projections.
      void recoverGradient();

      void computeWallShear();
      void closeCycle(Real elapsed);
      void computeFluxes();
      Real boundaryArea(Attribute tag);

      /// @brief h_K, the local element size.
      static Real cellSize(const Rodin::Geometry::Point& p);

      /**
       * @brief Local Carreau-Yasuda viscosity at the lagged velocity.
       * @details The stabilisation parameters must be built from the
       * viscosity the momentum equation actually sees. Sizing them
       * with the zero-shear plateau mu0 instead makes tau_1 too small
       * and tau_C = rho h^2 / (4 tau_1) too large by the same factor
       * -- on blood that is more than a decade -- and tau_C weighs an
       * explicitly lagged grad-div source.
       */
      Real viscosityAt(const Rodin::Geometry::Point& p) const;

      /// @brief Codina tau_1, c1 = 4, c2 = 2, k = 1 (P1), on the local viscosity.
      Real tau1At(const Rodin::Geometry::Point& p) const;

      /// @brief Convective VMS parameter.
      Real vmsTauAt(const Rodin::Geometry::Point& p) const;

      /// @brief Square root of the grad-div parameter.
      Real sqrtTauCAt(const Rodin::Geometry::Point& p) const;

      /// @brief Grad-div parameter, the exact square of the above.
      Real tauCAt(const Rodin::Geometry::Point& p) const;

      /// @brief PSPG parameter.
      Real tauPAt(const Rodin::Geometry::Point& p) const;

      /// @brief Codina crosswind coefficient of one species along u - w.
      Real crosswind(const Rodin::Geometry::Point& p, Real cur, Real prev,
        const Rodin::Math::SpatialVector<Real>& gradient, Real diffusivity) const;

      Parameters m_cfg;
      const RunParameters m_run;
      Optional<Diagnostics> m_diagnostics;
      Model::Input m_input;
      Model m_model;

      MeshType m_mesh;
      std::unique_ptr<Output> m_output;
      Rodin::FlatSet<Attribute> m_inletSet;

      VectorFESType m_vh;
      ScalarFESType m_sh;

      // ---- Flow ------------------------------------------------------------
      VectorTrialFunctionType m_u;
      ScalarTrialFunctionType m_p;
      VectorTestFunctionType m_v;
      ScalarTestFunctionType m_q;
      VectorGridFunctionType m_uOld;
      VectorGridFunctionType m_wallVelocity, m_wallDisplacement, m_wallTraction;
      VectorGridFunctionType m_meshVelocity, m_meshDisplacement;
      bool m_coupled = false;
      Real m_robinCoefficient = 0; ///< Reference-area impedance in Pa s/m.
      std::vector<Real> m_inverseInterfaceStretches;
      std::vector<Real> m_massRatios;

      /// @brief Nodal speed; its maximum is exact for affine P1 tetrahedral velocity.
      ScalarGridFunctionType m_speedMagnitude;

      // ---- Reusable mass projections ---------------------------------------
      ScalarTrialFunctionType m_sTrial;
      ScalarTestFunctionType m_sTest;
      VectorTrialFunctionType m_wTrial;
      VectorTestFunctionType m_wTest;

      // ---- VMS -------------------------------------------------------------
      // The taus are pointwise coefficients; only the two orthogonal-subscale
      // projections are genuine L2 projections.
      ScalarCoefficientType m_tauFn;
      ScalarCoefficientType m_tauCFn;
      ScalarCoefficientType m_sqrtTauCFn;
      ScalarCoefficientType m_tauPFn;
      ScalarCoefficientType m_viscosityFn;
      ScalarCoefficientType m_dragFn;
      mutable std::vector<CellCoefficients> m_cellCoefficients;
      size_t m_coefficientEpoch = 1;
      ScalarGridFunctionType m_piTilde;
      VectorGridFunctionType m_convProjection;
      VectorGridFunctionType m_sub;
      VectorGridFunctionType m_subOld;

      // ---- Species: P, T, A, G ---------------------------------------------
      ScalarTrialFunctionType m_pT, m_tT, m_aT, m_gT;
      ScalarTestFunctionType m_vP, m_vT, m_vA, m_vG;
      ScalarGridFunctionType m_pCur, m_tCur, m_aCur, m_gCur;
      ScalarGridFunctionType m_pPrev, m_tPrev, m_aPrev, m_gPrev;

      /// @brief F = G_0 - G, I = A_0 - A, sigma_B, and P + T - A - (P_0 - A_0).
      ScalarGridFunctionType m_fibrin, m_tat, m_drag, m_drift;

      // ---- Endocardial layer ------------------------------------------------
      ScalarTrialFunctionType m_eT;
      ScalarTestFunctionType m_vE;

      /// @brief Activation extended into the lumen; zero for zero layer thickness.
      ScalarGridFunctionType m_layer;

      /// @brief Dose normalization, delta int_Gw E dA / int_Omega E~ dV.
      Real m_layerScale = 1.0;

      /// @brief Assembled wall measure; positive exactly at wall nodes on every rank.
      ScalarGridFunctionType m_wallArea;
      Real m_wallAreaThreshold = 0.0;

      /// @brief Ghost-consistent copy of the recovered wall traction.
      VectorGridFunctionType m_wssRec;

      // ---- Wall shear indices ----------------------------------------------
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

      // ---- Time-dependent coefficients read by the frozen forms ------------
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

      /// @brief Cycle wall-shear indices are available.
      bool m_indicesReady = false;

      /// @brief Species have started; the Brinkman drag then enters the flow.
      bool m_speciesActive = false;

      /// @brief Step at which the retained flow preconditioner was last refreshed.
      Optional<int> m_lastPreconditionerRefreshStep;
      int m_preconditionerRefreshes = 0;
      int m_linearIterations = 0;
      bool m_preconditionerRefreshed = false;
      bool m_preconditionerRetried = false;
      bool m_initialized = false;
      bool m_finished = false;
      bool m_failed = false;
      int m_stepsPerCycle = 0;
      int m_lastOutputStep = 0;
      Real m_initialDt = 0;
      Real m_initialPeriod = 0;
      bool m_initialUseVentricularModel = true;
      Real m_cycleElapsed = 0;

      /// @brief Number of completed steps, or index of the step being solved.
      int m_step = 0;
      Timing m_timing;
      std::chrono::steady_clock::time_point m_stepStart;
  };
}

#endif
