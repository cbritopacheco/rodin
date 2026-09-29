/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file LeftAtrium2D_viscous_logarithmic_implicit.h
 * @brief LeftAtrium2D_viscous_logarithmic with the constitutive equation
 *        written in the log-conformation variable itself and solved by Newton
 *        (Knechtges, JNNFM 2015, fully-implicit log-conformation), instead of
 *        the exp(psi) form of Castillo et al. linearised by Picard.
 *
 * Same fluid, sigma = (eta_p/lambda0)(exp(psi) - I), sPTT factor f, same mesh,
 * data, P1 spaces, BDF1 everywhere and lagged split-OSS stabilisation. The
 * constitutive equation (21) of the exp form, divided by lambda/(2 lambda0)
 * and pulled through Dexp^{-1}, becomes
 *
 *   d psi/dt + (u.grad) psi = Dexp_psi^{-1}[L exp(psi) + exp(psi) L^T
 *       - 2 (1 - lambda0/lambda) eps(u)] - (f(psi)/lambda) (I - exp(-psi)),
 *
 * with Dexp_psi^{-1}[exp(psi) - I] = I - exp(-psi) exactly. Exponential growth
 * of the conformation is linear growth of psi here (the stretch term is
 * bounded by |grad u|, whatever psi is), and the relaxation I - exp(-psi) is
 * bounded and monotone: implicit Euler on this equation has no sign change at
 * 2 |grad u| dt > 1, unlike an implicit stretch in the exp form, and a Picard
 * fixed point on the exp form only contracts for dt |grad u| < 1, which the
 * ostium/wall junctions violate at peak inflow.
 *
 * Dexp^{-1} comes from the same eigen-decomposition as Dexp (1/F_ij instead of
 * F_ij), so there is no lambda_1 = lambda_2 special case as in the Fattal-
 * Kupferman decomposition. Newton: the pointwise nonlinear map
 * N(psi; L) = -Dexp^{-1}[...] + (f/lambda)(I - exp(-psi)) is differentiated in
 * psi by central differences (3 Voigt directions, 6 extra 2x2 eigenproblems per
 * quadrature point); its dependence on grad u is linear and kept exact and
 * implicit, and the advection is linearised as u^k.grad psi + u.grad psi^k -
 * u^k.grad psi^k. Momentum sees exp(psi) ~ exp(psi^k) + Dexp[psi - psi^k] as
 * before. The iteration stops on the relative psi increment; a step above
 * newtonMaxStep in max norm is scaled down to it.
 *
 * Stabilisation: S1 and S2 as before; the constitutive equation is now the psi
 * equation, which is 2x the sigma equation of the exp form in the linear
 * regime (psi/lambda0 - 2 eps(u) against psi/(2 lambda0) - eps(u)), so the
 * terms tested with chi carry that factor. The S3 convective pair acts on psi,
 * alpha_psi <P^perp[u.grad psi], u.grad chi> with alpha_psi = 2 k^2 s alpha3,
 * the weight the sigma pair had for sigma ~ s psi; the indefinite stretch pair
 * is dropped (the term-by-term form of Castillo & Codina 2016 does the same).
 * psi = 0 is imposed on the pulmonary vein ostia: the psi equation is
 * hyperbolic in psi and needs an inflow condition.
 */
#ifndef EXAMPLES_HEART_LEFTATRIUM2D_VISCOUS_LOGARITHMIC_IMPLICIT_H
#define EXAMPLES_HEART_LEFTATRIUM2D_VISCOUS_LOGARITHMIC_IMPLICIT_H

#include <array>
#include <cstdint>
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
#include "OrthogonalProjection.h"

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

      Real getPeriod() const noexcept
      {
        return m_period;
      }
      Real getMinimum() const noexcept
      {
        return m_min;
      }
      Real getMaximum() const noexcept
      {
        return m_max;
      }
      Real getMean() const noexcept
      {
        return m_mean;
      }
      size_t getSampleCount() const noexcept
      {
        return m_t.size();
      }
      const std::string& getPath() const noexcept
      {
        return m_path;
      }

    private:
      std::string m_path;
      std::vector<Real> m_t;
      std::vector<Real> m_p;
      Real m_period = 0.0;
      Real m_min = 0.0;
      Real m_max = 0.0;
      Real m_mean = 0.0;
  };

  class LeftAtrium2DViscousLogImplicit
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

      using VectorTrialFunctionType =
        Rodin::PETSc::Variational::TrialFunction<VectorGridFunctionType, VectorFESType>;
      using ScalarTrialFunctionType =
        Rodin::PETSc::Variational::TrialFunction<ScalarGridFunctionType, ScalarFESType>;
      using VectorTestFunctionType =
        Rodin::PETSc::Variational::TestFunction<VectorFESType>;
      using ScalarTestFunctionType =
        Rodin::PETSc::Variational::TestFunction<ScalarFESType>;

      using LinearSystemType = Rodin::PETSc::Math::LinearSystem;

      /// @brief (u, p, sigma; v, q, chi). sigma lives in a VectorFESType of
      ///        dimension d(d+1)/2 (Voigt), so it shares the vector types.
      using FlowProblemType =
        Rodin::Variational::Problem<LinearSystemType, VectorTrialFunctionType,
          ScalarTrialFunctionType, VectorTrialFunctionType, VectorTestFunctionType,
          ScalarTestFunctionType, VectorTestFunctionType>;

      using SpeciesProblemType = Rodin::Variational::Problem<LinearSystemType,
        ScalarTrialFunctionType, ScalarTrialFunctionType, ScalarTrialFunctionType,
        ScalarTestFunctionType, ScalarTestFunctionType, ScalarTestFunctionType>;

      using VectorProblemType = Rodin::Variational::Problem<LinearSystemType,
        VectorTrialFunctionType, VectorTestFunctionType>;

      using FluxFormType = Rodin::Variational::LinearForm<ScalarFESType, ::Vec>;

      /// @brief Cellwise (P0) stabilisation parameters.
      using CellFESType = Rodin::Variational::P0<Real, MeshType>;
      using CellGridFunctionType = Rodin::PETSc::Variational::GridFunction<CellFESType>;

      /// @brief Oldroyd-B blood: T = 2 eta_s eps(u) + sigma, eta_s = beta eta_0,
      ///        eta_p = (1 - beta) eta_0 (Castillo et al., Eqs. 3-4).
      /// @details Elhanafy, Guaily & Elsaid, Adv. Mech. Eng. 11 (2019), Table 1,
      ///          after Bodnar et al.: eta_0 = 3.59e-3 Pa s, beta = 0.889.
      struct OldroydB
      {
          Real etaS = 3.19e-3;  ///< solvent viscosity eta_s (Pa s)
          Real etaP = 4.0e-4;   ///< polymeric viscosity eta_p (Pa s)
          Real lambda = 0.06;   ///< relaxation time lambda (s)
          /// @brief lambda0 = max(k lambda, lambda0_min), Sec. 3.2; k = 1 is
          ///        Fattal-Kupferman, k small is reported to converge better.
          Real lambda0Factor = 1.0;
          Real lambda0Min = 0.0;
          /// @brief Linear sPTT extensibility parameter epsilon; 0 recovers
          ///        Oldroyd-B. Typical values 0.02-0.25.
          Real pttEpsilon = 0.1;
      };

      /**
       * @brief Boundary labels of RectLAA.mesh.
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
          /// @brief The FLATTENED mesh, produced by LeftAtrium2D/make_la2d_mesh.py.
          std::string meshPath =
            "../resources/examples/Heart/LeftAtrium2D/RectLAA2D.mesh";
          /// @brief Pulmonary venous pressure, imposed on the four PV ostia.
          std::string inletPressurePath =
            "../resources/examples/Heart/LeftAtrium2D/PulmonaryVeinPressureSR.dat";
          /// @brief Ventricular pressure, imposed on the mitral orifice.
          std::string outletPressurePath =
            "../resources/examples/Heart/LeftAtrium2D/MitralValvePressureSR.dat";

          std::string xdmfBasename = "LeftAtrium2D_viscous_logarithmic_implicit";
          std::string csvPath = "LeftAtrium2D_viscous_logarithmic_implicit.csv";

          Labels labels;

          /// @brief The mesh is already in metres.
          Real meshScale = 1.0;

          Real rho = 1060.0;
          OldroydB oldroydB;
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
          Real vmsScale = 1.0;       ///< S1, dynamic convection
          Real gradDivScale = 1.0;   ///< S2, alpha2
          Real pressureScale = 1.0;  ///< S1, alpha1 <P^perp[grad p], grad q>
          Real stressDivScale = 1.0; ///< S1, alpha1 <P^perp[div sigma], div chi>
          Real stressScale = 1.0;    ///< S3, alpha3

          /// @brief Enable the VMS convection and grad-div stabilisation.
          /// @details At false, tau_K and alpha2 are set to zero, so the S1
          ///          convection and S2 terms vanish. The pressure and stress
          ///          terms are kept: the
          ///          equal-order u-p and u-sigma pairs are not stable without
          ///          them.
          bool useVMS = true;

          /// @brief Maximum Newton iterations per step.
          int conformationIterations = 6;
          /// @brief Stop when max|psi^{k+1} - psi^k| / max(1, max|psi^k|) is
          ///        below this; a stall is reported only above 100x it.
          Real conformationTolerance = 1.0e-6;
          /// @brief A Newton step larger than this in max|dpsi| is scaled
          ///        down to it (a psi step of 2 is a factor e^2 in tau).
          Real newtonMaxStep = 2.0;
          /// @brief Impose psi = 0 (relaxed blood) on the pulmonary vein ostia.
          bool inletConformation = true;

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

      LeftAtrium2DViscousLogImplicit(const Rodin::Context::MPI& context, const Config& cfg);
      ~LeftAtrium2DViscousLogImplicit();

      LeftAtrium2DViscousLogImplicit(const LeftAtrium2DViscousLogImplicit&) = delete;
      LeftAtrium2DViscousLogImplicit& operator=(const LeftAtrium2DViscousLogImplicit&) = delete;

      LeftAtrium2DViscousLogImplicit& initialize();
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
      static AttributeSet makeInletSet(const Config& cfg);
      static AttributeSet makeWallSet(const Config& cfg);

      bool isRoot() const;

      void setupSpaces();
      void setupFlow();
      void setupSpecies();
      void setupWallShear();

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

      /// @brief Codina tau_1, c1 = 4, c2 = 2, k = 1 (P1), on the solvent
      ///        viscosity eta_s: the only viscosity the momentum equation sees.
      Real tau1At(const Rodin::Geometry::Point& p) const;
      /// @brief alpha3 = [c3 f/(2 eta_p) + c4 (lambda|u|/(2 eta_p h)
      ///        + lambda|grad u|/eta_p)]^-1, c3 = 4, c4 = 1/4, Eq. (42), with
      ///        the PTT factor f at psi^n on the relaxation term.
      Real alpha3At(const Rodin::Geometry::Point& p) const;

      /// @brief The alphas (P0) and every Pi_h, all at step n.
      void updateStabilization();

      /// @brief lambda0 = max(k lambda, lambda0_min).
      Real lambda0() const;

      /// @brief PTT factor f = 1 + epsilon (lambda/lambda0)(tr exp(psi) - d).
      Real pttFactor(Real traceExp) const;

      /// @brief sigma at the nodes, from m_psiIt.
      void updateConformation();

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
      /// @brief Tensor space: P1, Voigt order (xx, xy, yy).
      VectorFESType m_tauh;
      CellFESType m_ch;

      // ---- Flow ------------------------------------------------------------
      VectorTrialFunctionType m_u;
      ScalarTrialFunctionType m_p;
      VectorTestFunctionType m_v;
      ScalarTestFunctionType m_q;
      VectorGridFunctionType m_uOld;
      /// @brief Newton velocity u^k of the constitutive equation.
      VectorGridFunctionType m_uIt;
      /// @brief psi = log(tau) (Voigt), the stress test function chi, psi^n.
      VectorTrialFunctionType m_psi;
      VectorTestFunctionType m_chi;
      VectorGridFunctionType m_psiOld;
      /// @brief Linearisation point psi^k of the Newton iterations.
      VectorGridFunctionType m_psiIt;
      /// @brief Bumped whenever (u^k, psi^k) or (u^n, psi^n) change; keys the
      ///        cache of the pointwise log-conformation data.
      std::uint64_t m_psiRevision = 0;
      /// @brief sigma = (eta_p/lambda0)(exp(psi) - I) at the nodes (Voigt), for
      ///        output and the wall traction (updateConformation).
      VectorGridFunctionType m_sigma;

      // ---- Reusable mass projection (wall-shear recovery) -------------------
      VectorTrialFunctionType m_wTrial;
      VectorTestFunctionType m_wTest;

      // ---- VMS, split OSS (updateStabilization) ----------------------------
      CellGridFunctionType m_tauK;    ///< dynamic convective tau_K
      CellGridFunctionType m_alpha1;  ///< tau_1/rho
      CellGridFunctionType m_alpha2;  ///< rho h^2/(4 tau_1)
      CellGridFunctionType m_alpha3;
      CellGridFunctionType m_alphaPsi; ///< 2 k^2 s alpha3, the psi advection weight
      OrthogonalProjection<VectorFESType> m_piConv;     ///< Pi[(grad u^n) u^n]
      OrthogonalProjection<VectorFESType> m_sub;        ///< dynamic subscale u'
      VectorGridFunctionType m_subOld;
      OrthogonalProjection<ScalarFESType> m_piDiv;      ///< Pi[div u^n]
      OrthogonalProjection<VectorFESType> m_piGradP;    ///< Pi[grad p^n]
      OrthogonalProjection<VectorFESType> m_piDivSigma; ///< Pi[div sigma^n]
      OrthogonalProjection<VectorFESType> m_piEps;      ///< Pi[eps(u^n)], Voigt
      OrthogonalProjection<VectorFESType> m_piAdvPsi;   ///< Pi[(u^n.grad) psi^n], Voigt

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
      /// @brief max |sigma_ij| (Pa), reported with the flow.
      Real m_stress = 0.0;
      /// @brief Newton iterations of the last step and their last increment.
      int m_conformationIts = 0;
      Real m_psiIncrement = 0.0;
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
