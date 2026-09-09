// LeftAtrium2D.h
//
// Two-dimensional left atrium with a rectangular appendage: stabilised
// non-Newtonian flow driven by two *measured* pressure waveforms, plus the
// thrombin/fibrinogen/fibrin kinetics transported by that flow.
//
//   inlets (PV) : p = p_pv(t)   read from resources/.../presion_PV_SR.dat
//   outlet (MV) : p = p_v(t)    read from resources/.../presion_MV_SR.dat
//   wall + LAA  : no slip, endothelial thrombin flux weighted by the cycle
//                 wall-shear indices.
//
// This is the 2D counterpart of Atrium.{h,cpp}.  The differences that matter
// are not cosmetic:
//
//   1. The mesh is genuinely two-dimensional.  LA2D_rectLAA.mesh is written as
//      MEDIT "Dimension 3" with z = 0, and Rodin reads the space dimension
//      straight from that keyword, so loading it unchanged gives a 2D surface
//      living in 3D: a three-component velocity and a boundary normal that is
//      not the in-plane one.  LA2D/make_la2d_mesh.py flattens it once; this
//      example refuses anything whose space dimension is not 2.
//
//   2. There is no 0D ventricle.  Atrium closes the mitral outlet with the
//      CCMLC2014 model and a diode.  Here both pressures are tabulated, and
//      p_pv - p_mv is identically zero over the intervals where the valve is
//      shut, so the waveform pair already carries the valve.  Adding a diode
//      or a mitral resistance on top of it would impose the same closure
//      twice.
//
//   3. The coordinates are already in metres (the atrium is 56 mm across), so
//      meshScale is 1, not 1e-3.
//
//   4. The driving pressure difference peaks at 30 Pa.  Every boundary
//      impedance therefore has to be small compared with dp/U ~ 150 Pa s/m or
//      it becomes the throttle that sets the flow rate; see Config.
//
// Everything the cycle indices are built from is evaluated at the *nodes*.
// TAWSS, OSI and the activation weight are nonlinear functions of two
// accumulators, and taking them through an L2 projection -- i.e. through
// quadrature -- is what produces negative |tau_w|, negative TAWSS and an
// activation that turns the endothelial thrombin flux into a sink.  The taus
// are pointwise coefficients for the same reason.  The only genuine L2
// projections left are the gradient recovery and the two VMS orthogonal
// subscales, where a projection is what the method actually asks for.
//
// They are also written ONLY on wall nodes.  Off the wall tau_w decays to
// zero, the activation logistic saturates at 1/(1+exp(-tau_a/w)) = 0.935 and
// OSI degenerates into a 0/0 that drifts to 1/2, so a whole-domain index field
// peaks in the middle of the cavity where there is no endothelium -- 0.935 x 1
// = 0.933, which is what the first run reported as its maximum activation.
//
// Equal order P1/P1, stabilised with VMS (convection, grad-div) and PSPG.
// Runs under MPI; the mesh is partitioned on the root rank.
#ifndef EXAMPLES_HEART_LEFTATRIUM2D_H
#define EXAMPLES_HEART_LEFTATRIUM2D_H

#include <array>
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
  /**
   * @brief A tabulated, periodic pressure waveform.
   *
   * @details Two whitespace-separated columns, t [s] and p [Pa], with '#'
   *          comments. The samples must be strictly increasing in t. The
   *          signal is extended periodically with period t_last - t_first, and
   *          interpolated linearly in between: the files carry 401 samples at
   *          2 ms, so a linear reconstruction is well inside the sampling
   *          error of the data itself, and it cannot overshoot the way a cubic
   *          can at the valve-closure corners.
   */
  class PressureWaveform
  {
    public:
      using Real = Rodin::Real;

      PressureWaveform() = default;

      /// @brief Reads the file, or throws if it cannot be parsed.
      void load(const std::string& path);

      /// @brief p(t), extended periodically.
      Real operator()(Real t) const;

      Real getPeriod() const noexcept { return m_period; }
      Real getMinimum() const noexcept { return m_min; }
      Real getMaximum() const noexcept { return m_max; }
      Real getMean() const noexcept { return m_mean; }
      size_t getSampleCount() const noexcept { return m_t.size(); }
      const std::string& getPath() const noexcept { return m_path; }

    private:
      std::string m_path;
      std::vector<Real> m_t;
      std::vector<Real> m_p;
      Real m_period = 0.0;
      Real m_min = 0.0;
      Real m_max = 0.0;
      Real m_mean = 0.0;
  };

  class LeftAtrium2D
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

      /**
       * @brief Boundary labels of LA2D_rectLAA.mesh.
       *
       * @details Read off the file, not assumed. The boundary is a closed
       *          polyline split into fourteen tagged pieces; going clockwise
       *          from the mitral cut:
       *
       *            14  straight segment, x = +23.5 mm, 25.0 mm long -> MV
       *             1  arc from the MV to the base of the appendage
       *           2,3,4  the three straight sides of the rectangular LAA
       *             5  arc between the LAA and the first vein
       *          6,8,10,12  four arcs, 10.3/9.0/9.0/10.3 mm -> the four PVs
       *           7,9,11  the wall segments separating them
       *            13  the whole inferior arc back to the MV
       *
       *          6/12 and 8/10 are mirror pairs about y = 0 and 7/9/11 sit
       *          between them, which is what fixes the alternation: the other
       *          reading, veins on 5/7/9/11, puts a vein ostium flush against
       *          the appendage wall and leaves the leftmost 4.8 mm arc as a
       *          vein.
       */
      struct Labels
      {
          /// @brief Pulmonary vein ostia (pressure inlets).
          std::array<Attribute, 4> inlets{{6, 8, 10, 12}};
          /// @brief Mitral orifice (pressure outlet).
          Attribute outlet = 14;
          /// @brief Atrial body wall (no slip).
          std::array<Attribute, 6> wall{{1, 5, 7, 9, 11, 13}};
          /// @brief Left atrial appendage wall (no slip). Kept apart from the
          ///        body only so the appendage can be selected in a report;
          ///        both carry the same condition and the same thrombin flux.
          std::array<Attribute, 3> appendage{{2, 3, 4}};
      };

      struct Config
      {
          /// @brief The FLATTENED mesh, produced by LA2D/make_la2d_mesh.py.
          std::string meshPath =
            "../resources/examples/Heart/output.mesh";
          /// @brief Pulmonary venous pressure, imposed on the four PV ostia.
          std::string inletPressurePath =
            "../resources/examples/Heart/presion_PV_SR.dat";
          /// @brief Ventricular pressure, imposed on the mitral orifice.
          std::string outletPressurePath =
            "../resources/examples/Heart/presion_MV_SR.dat";

          std::string xdmfBasename = "LeftAtrium2D";
          std::string csvPath = "LeftAtrium2D.csv";

          Labels labels;

          /// @brief The mesh is already in metres.
          Real meshScale = 1.0;

          Real rho = 1060.0;
          CarreauYasuda viscosity;
          /// @brief Pressure penalty of the equal-order pair.
          Real pressurePenalty = 1.0e-12;

          /// @brief Constant offset added to both waveforms (Pa).
          /// @details Only the difference drives the flow, so this shifts the
          ///          pressure level without touching the solution.
          Real pressureOffset = 0.0;

          /// @brief Normal impedance at the pressure-driven inlets (Pa s/m).
          /// @details Scale first, then choose. The waveforms differ by at most
          ///          30 Pa and the expected inflow is a couple of tenths of a
          ///          metre per second, so anything approaching 30/0.24 = 125
          ///          Pa s/m stops being a regularisation and becomes the
          ///          resistance that sets the flow rate -- which is why
          ///          Atrium's 1e3 cannot be carried over. 20 spends about a
          ///          sixth of the head at 0.24 m/s and stands in for the
          ///          viscous resistance of the veins themselves. It is a
          ///          regularisation, not the energy bound: that is the
          ///          incoming-kinetic-energy term below.
          Real inletImpedance = 20.0;
          /// @brief Tangential damping at the inlets (Pa s/m).
          /// @details Light: it only discourages a vein ostium from acting as a
          ///          free-slip window, and at 0.1 m/s it costs 5 Pa.
          Real inletTangentialDamping = 50.0;
          /// @brief Weight of the incoming-kinetic-energy term on the inlets.
          /// @details 1 is the directional do-nothing condition and the only
          ///          value with an unconditional energy estimate. It also
          ///          turns the prescribed pressure into a TOTAL pressure:
          ///          the boundary then carries p + rho|u_n|^2/2. At the
          ///          velocity this problem targets that shift is worth 30 Pa,
          ///          the whole driving head, so it is a modelling choice and
          ///          not a numerical detail. Lower values weaken the bound;
          ///          0 recovers a pure static-pressure inlet, which has no
          ///          energy estimate at all.
          Real inletBackflowStabilization = 1.0;
          Real outletBackflowStabilization = 1.0;

          Real dt = 1.0e-3;
          /// @brief Cycle length (s). Checked against the waveform files.
          Real period = 0.8;
          /// @brief Cycles solved before the species are switched on.
          int flowCycles = 2;
          int speciesCycles = 10;

          /// @brief XDMF output period, in steps. Every cycle boundary is
          ///        written regardless, so 0 keeps only those.
          int outputEvery = 20;

          /// @brief Per-term stabilisation scales, so each can be dialled down
          ///        on its own. 1 is the standard value, 0 removes the term.
          Real vmsScale = 1.0;
          Real gradDivScale = 1.0;
          Real pspgScale = 1.0;

          /// @brief Enable the VMS convection and grad-div stabilisation.
          /// @details At false, tau and tau_C are returned as zero, so the four
          ///          VMS integrators contribute nothing and their projections
          ///          are skipped. PSPG is kept either way: the equal-order
          ///          pair is not stable without it. Diagnostic switch --
          ///          these are the only terms in this example that no other
          ///          MPI example in the tree exercises.
          bool useVMS = true;

          /// @brief Weight of the transient and convective halves of the PSPG
          ///        momentum residual. 1 is the consistent method, 0 the
          ///        pressure-gradient-only penalty.
          /// @details PSPG is tau_p grad(q).R_M with
          ///          R_M = rho (u - u^n)/dt + rho (grad u) u^n + grad p.
          ///          Keeping only the last term, as Atrium does, is a
          ///          Brezzi-Pitkaranta penalty: it stabilises the pair, but it
          ///          does not vanish on the exact solution, so it perturbs the
          ///          pressure instead of being consistent. The viscous part of
          ///          R_M is dropped because it involves second derivatives,
          ///          which vanish elementwise on P1.
          ///
          ///          A scale rather than a flag because the terms are written
          ///          into a single form-language expression that is assigned
          ///          once: at 0 they assemble to zero, which is what the
          ///          comparison against Atrium's form needs.
          Real pspgResidualScale = 1.0;

          ThrombosisParameters thrombosis;
          /// @brief Codina crosswind constant. 0 disables it.
          Real crosswindC = 0.7;
          bool solveKinetics = true;

          /// @brief Thrombin carried in by the pulmonary venous blood
          ///        (mol/m^3).
          /// @details Atrium sends in thrombin-free blood, so every molecule in
          ///          the domain comes from the wall flux. Here the veins carry
          ///          a small circulating level -- 1 nmol/m^3 -- so that
          ///          fibrinogen is consumed everywhere at a slow background
          ///          rate and the appendage signal is read against it rather
          ///          than against an exact zero.
          Real inletThrombin = 1.0e-6;
          /// @brief Fibrin carried in by the pulmonary venous blood (mol/m^3).
          Real inletFibrin = 0.0;

          /// @brief Divergence guard on max|u| (m/s).
          /// @details Judge it against sqrt(2 max|dp|/rho) = 0.24 m/s, which is
          ///          all the two waveforms can pay for; setupSpaces() prints
          ///          the ratio. 2.0 was 8x that -- the first run reached 1.22
          ///          m/s, a dynamic head of 794 Pa against a 30 Pa drive, and
          ///          carried on for hundreds of steps before tripping.
          Real maxVelocity = 1.0;
      };

      /// @brief Wall-clock cost of each phase of one step.
      struct Timing
      {
          Real vms = 0.0;
          Real assembly = 0.0;
          Real solve = 0.0;
          Real shear = 0.0;
          Real fluxes = 0.0;
          Real species = 0.0;
          Real output = 0.0;
          Real total = 0.0;
      };

      LeftAtrium2D(const Rodin::Context::MPI& context, const Config& cfg);
      ~LeftAtrium2D();

      LeftAtrium2D(const LeftAtrium2D&) = delete;
      LeftAtrium2D& operator=(const LeftAtrium2D&) = delete;

      LeftAtrium2D& initialize();
      int run();

      Config& getConfig() noexcept { return m_cfg; }
      const Config& getConfig() const noexcept { return m_cfg; }

    private:
      static MeshType makeMesh(const Rodin::Context::MPI& context, const Config& cfg);
      static AttributeSet makeInletSet(const Config& cfg);
      static AttributeSet makeWallSet(const Config& cfg);

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

      bool solveFlow();
      void solveSpecies();
      void computeWallShear();
      void closeCycle(Real elapsed);
      void computeFluxes();
      Real boundaryMeasure(const AttributeSet& tags);

      /// @brief h_K, the local element size.
      static Real cellSize(const Rodin::Geometry::Point& p);

      /// @brief Local Carreau-Yasuda viscosity at the lagged velocity.
      /// @details The stabilisation parameters must be built from the
      ///          viscosity the momentum equation actually sees. Sizing them
      ///          with the zero-shear plateau mu0 instead makes tau_1 too small
      ///          and tau_C = rho h^2 / (4 tau_1) too large by the same factor.
      Real viscosityAt(const Rodin::Geometry::Point& p) const;

      /// @brief Codina tau_1, c1 = 4, c2 = 2, k = 1 (P1), on the local viscosity.
      Real tau1At(const Rodin::Geometry::Point& p) const;
      /// @brief Convective VMS parameter.
      Real vmsTauAt(const Rodin::Geometry::Point& p) const;
      /// @brief Square root of the grad-div parameter.
      Real sqrtTauCAt(const Rodin::Geometry::Point& p) const;
      /// @brief Grad-div parameter, the exact square of the above.
      Real tauCAt(const Rodin::Geometry::Point& p) const;
      /// @brief PSPG parameter, tau_1/rho.
      Real tauPAt(const Rodin::Geometry::Point& p) const;

      /// @brief Codina crosswind coefficient of one species.
      Real crosswind(const Rodin::Geometry::Point& p, Real cur, Real prev,
        const Rodin::Math::SpatialVector<Real>& gradient, Real diffusivity,
        Real reaction) const;

      void writeCSVHeader();
      void writeCSVRow(int cycle);

      Config m_cfg;

      PressureWaveform m_inletWave;
      PressureWaveform m_outletWave;

      MeshType m_mesh;
      Rodin::IO::XDMF m_xdmf;
      AttributeSet m_inletSet;
      AttributeSet m_wallSet;

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
      // symRec0 and symRec1 are the two rows of 2 eps(u) recovered onto the
      // nodes, not the rows of grad u: the viscous traction is mu (grad u +
      // grad u^T) n, and although the transposed half has no tangential
      // component on a straight no-slip wall, it does on the LAA corners and
      // on the discrete wall in general.
      VectorGridFunctionType m_wss;
      VectorGridFunctionType m_symRec0;
      VectorGridFunctionType m_symRec1;
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
      Real m_outletMeasure = 1.0;
      Real m_inletMeasure = 1.0;
      Real m_outletPressure = 0.0;
      Real m_qIn = 0.0;
      Real m_qOut = 0.0;
      /// @brief Flux through each PV ostium separately.
      /// @details Reported because the aggregate cannot tell an evenly fed
      ///          atrium from one where two ostia carry everything, and the
      ///          first velocity field showed exactly two jets. Also the
      ///          cheapest check on the label assignment: a patch this example
      ///          calls a vein but the mesh calls wall carries no flux.
      std::array<Real, 4> m_qInPatch{{0.0, 0.0, 0.0, 0.0}};
      Real m_speed = 0.0;
      /// @brief sqrt(2 max|dp| / rho): the largest velocity the two waveforms
      ///        can account for. Printed at startup and used by the guard.
      Real m_velocityScale = 0.0;
      bool m_flowFieldSplitsSet = false;
      bool m_initialized = false;

      /// @brief Set by the first closeCycle(). The kinetics are not started
      ///        before it: the endothelial thrombin flux is weighted by the
      ///        activation field, and that field is a function of OSI and
      ///        TAWSS, which do not exist until a whole cycle has been
      ///        accumulated.
      bool m_indicesReady = false;

      /// @brief Index of the step being solved; the first ones are traced
      ///        phase by phase, so a stall is visible where it happens.
      int m_step = 0;
      Timing m_timing;

      std::ofstream m_csv;
  };
}

#endif
