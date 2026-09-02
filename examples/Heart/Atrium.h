/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Atrium.h
 * @brief Stabilised non-Newtonian atrium flow coupled to the 0D left ventricle.
 *
 * 3D atrium: stabilised non-Newtonian flow with four pressure inlets and one
 * outlet closed by the 0D left ventricle through a mitral diode, plus the
 * thrombin/fibrinogen/fibrin kinetics transported by that flow.
 *
 *   inlets  : p = p_pv(t)             (natural traction + impedance + damping)
 *   outlet  : p = p_v(t) from the 0D LV, with R_mv A (u.n)(v.n) assembled
 *             implicitly; R_mv switches between open and closed values.
 *   wall    : no slip, endothelial thrombin flux weighted by the cycle
 *             wall-shear indices.
 *
 * Equal order P1/P1, stabilised with VMS (convection, grad-div) and PSPG.
 * Runs under MPI; the mesh is partitioned on the root rank and scaled by
 * Config::meshScale.
 */
#ifndef EXAMPLES_HEART_ATRIUM_H
#define EXAMPLES_HEART_ATRIUM_H

#include <array>
#include <fstream>
#include <functional>
#include <string>

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
  class Atrium
  {
    public:
      using Real = Rodin::Real;
      using Model = Rodin::Heart::CCMLC2014T<>;
      using Attribute = Rodin::Geometry::Attribute;
      using MeshType = Rodin::Geometry::Mesh<Rodin::Context::MPI>;

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

      using SpeciesProblemType = Rodin::Variational::Problem<LinearSystemType,
        ScalarTrialFunctionType, ScalarTrialFunctionType, ScalarTrialFunctionType,
        ScalarTestFunctionType, ScalarTestFunctionType, ScalarTestFunctionType>;

      using VectorProblemType = Rodin::Variational::Problem<LinearSystemType,
        VectorTrialFunctionType, VectorTestFunctionType>;

      using ScalarProblemType = Rodin::Variational::Problem<LinearSystemType,
        ScalarTrialFunctionType, ScalarTestFunctionType>;

      using FluxFormType = Rodin::Variational::LinearForm<ScalarFESType, ::Vec>;

      /// @brief A scalar coefficient evaluated pointwise.
      /// @details The stabilisation parameters are element-local quantities and
      ///          are evaluated where they are used, never L2-projected onto a
      ///          continuous field: the projection of a positive function that
      ///          varies by a decade between neighbouring cells overshoots, and
      ///          a tau that comes out negative at a node turns a
      ///          positive-semidefinite stabilisation into an indefinite one.
      ///          It is also five mass solves per step cheaper.
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

      struct Config
      {
          std::string meshPath = "../resources/examples/Heart/Atrium/Coarse.mesh";
          std::string xdmfBasename = "Atrium";
          std::string csvPath = "Atrium.csv";

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

          /// @brief Normal impedance at the pressure-driven inlets (Pa s/m).
          /// @details A pressure inlet with zero source impedance has no energy
          ///          bound: the convective boundary flux grows like |u|^3 while
          ///          only viscous dissipation opposes it.
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
          /// @brief Cycles solved before the species are switched on.
          int flowCycles = 2;
          int speciesCycles = 10;

          /// @brief XDMF output period, in steps. Every cycle boundary is
          ///        written regardless, so 0 keeps only those.
          int outputEvery = 10;

          /// @brief Per-term stabilisation scales, so each can be dialled down
          ///        on its own. 1 is the standard value, 0 removes the term.
          Real vmsScale = 1.0;
          Real gradDivScale = 1.0;
          Real pspgScale = 1.0;

          /// @brief Enable the VMS convection and grad-div stabilisation.
          /// @details Disabling it leaves tau = tau_C = 0, so the custom
          ///          integrators contribute nothing; PSPG is kept, since the
          ///          equal-order pair is not stable without it. Diagnostic
          ///          switch: these are the only terms in this example that no
          ///          other MPI example in the tree exercises.
          bool useVMS = true;

          ThrombosisParameters thrombosis;
          /// @brief Codina crosswind constant. 0 disables it.
          Real crosswindC = 0.7;
          bool solveKinetics = true;

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

      Atrium(const Rodin::Context::MPI& context, const Config& cfg);
      ~Atrium();

      Atrium(const Atrium&) = delete;
      Atrium& operator=(const Atrium&) = delete;

      Atrium& initialize();
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

      bool isRoot() const;

      void setupSpaces();
      void setupFlow();
      void setupSpecies();
      void setupWallShear();

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
      void computeWallShear();
      void closeCycle(Real elapsed);
      void computeFluxes();
      Real boundaryArea(Attribute tag);

      /// @brief h_K, the local element size.
      static Real cellSize(const Rodin::Geometry::Point& p);

      /// @brief Local Carreau-Yasuda viscosity at the lagged velocity.
      /// @details The stabilisation parameters must be built from the
      ///          viscosity the momentum equation actually sees. Sizing them
      ///          with the zero-shear plateau mu0 instead makes tau_1 too small
      ///          and tau_C = rho h^2 / (4 tau_1) too large by the same factor
      ///          -- on blood that is more than a decade -- and tau_C weighs an
      ///          explicitly lagged grad-div source.
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

      /// @brief Codina crosswind coefficient of one species.
      Real crosswind(const Rodin::Geometry::Point& p, Real cur, Real prev,
        const Rodin::Math::SpatialVector<Real>& gradient, Real diffusivity,
        Real reaction) const;
      void writeCSVHeader();
      void writeCSVRow(int cycle);

      Config m_cfg;
      Model::Input m_input;
      Model m_model;

      MeshType m_mesh;
      Rodin::IO::XDMF m_xdmf;
      Rodin::FlatSet<Attribute> m_inletSet;

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

      // ---- VMS -------------------------------------------------------------
      // The taus are pointwise coefficients; only the two orthogonal-subscale
      // projections are genuine L2 projections.
      ScalarCoefficientType m_tauFn;
      ScalarCoefficientType m_tauCFn;
      ScalarCoefficientType m_sqrtTauCFn;
      ScalarCoefficientType m_tauPFn;
      ScalarGridFunctionType m_piTilde;
      VectorGridFunctionType m_convProjection;
      VectorGridFunctionType m_sub;
      VectorGridFunctionType m_subOld;

      // ---- Species ---------------------------------------------------------
      ScalarTrialFunctionType m_th;
      ScalarTrialFunctionType m_fg;
      ScalarTrialFunctionType m_fn;
      ScalarTestFunctionType m_vth;
      ScalarTestFunctionType m_vfg;
      ScalarTestFunctionType m_vfn;
      ScalarGridFunctionType m_thCur;
      ScalarGridFunctionType m_fgCur;
      ScalarGridFunctionType m_fnCur;
      ScalarGridFunctionType m_thPrev;
      ScalarGridFunctionType m_fgPrev;
      ScalarGridFunctionType m_fnPrev;

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
      SpeciesProblemType m_species;
      Rodin::Solver::KSP m_speciesKSP;
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
      bool m_initialized = false;

      /// @brief Index of the step being solved; the first ones are traced
      ///        phase by phase, so a stall is visible where it happens.
      int m_step = 0;
      Timing m_timing;

      std::ofstream m_csv;
  };
}

#endif
