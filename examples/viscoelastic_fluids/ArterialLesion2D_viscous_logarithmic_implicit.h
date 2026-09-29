/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ArterialLesion2D_viscous_logarithmic_implicit.h
 * @brief Pulsatile sPTT/Oldroyd-B flow through an idealised stenosis or
 *        fusiform aneurysm, in the fully-implicit log-conformation form of
 *        Heart/LeftAtrium2D_viscous_logarithmic_implicit.
 *
 * The fluid problem is the one of the left-atrium example, term by term:
 * T = 2 eta_s eps(u) + sigma, sigma = (eta_p/lambda0)(exp(psi) - I), the psi
 * equation solved by Newton, P1 spaces, BDF1, and the lagged split-OSS
 * stabilisation (S1 dynamic convection, pressure and div sigma; S2 grad-div;
 * S3 u-psi compatibility and psi advection). There are no coagulation species.
 * Only the domain and the boundary conditions change.
 *
 * Domain. A planar channel of width D whose walls follow
 *
 *   h(z)/R = 1 + a f(z),  f(z) = (1 + cos(pi z/ell))/2 for |z| <= ell,
 *
 * with a = sqrt(1 - S) - 1 (stenosis) or a = Gam - 1 (aneurysm); the meshes
 * are written in units of D by make_lesion_mesh.py and scaled by D here. The
 * formulation is planar Cartesian: it is NOT the axisymmetric pipe, which needs
 * the hoop components of eps(u) and psi and r-weighted integrals. In the plane
 * the same profile gives a width reduction 1 - sqrt(1 - S), not an area one.
 *
 * Inlet. u = U(y, t) e_x, the fully developed pulsatile channel flow of the
 * same fluid at the imposed flow rate Q(t) = Ubar D q(t). q is expanded in
 * Fourier modes; mode 0 is Poiseuille and mode n >= 1 is the Womersley channel
 * profile with the complex viscosity eta*(n omega) = eta_s + eta_p/(1 + i n
 * omega lambda). For unidirectional flow the upper-convected terms leave the
 * shear stress on the linear Maxwell equation, so the profile is exact for
 * Oldroyd-B and the linear-viscoelastic limit for sPTT. psi = 0 is imposed on
 * the inlet (relaxed fluid), as on the pulmonary veins of the parent example.
 *
 * Outlet. Traction p_out n with p_out constant, plus the directional
 * do-nothing (backflow) term. With rigid walls and a single outlet the
 * velocity does not depend on p_out: it shifts the pressure by a constant.
 *
 * Dimensionless groups, with eta_0 = eta_s + eta_p:
 *   Re = rho Ubar D/eta_0,  Wo = (D/2) sqrt(omega rho/eta_0),
 *   Wi = lambda Ubar/D,     De = lambda/T,
 * and Wi/De = Ubar T/D = (pi/2) Re/Wo^2. Given D, rho, eta_s and eta_p, the
 * run is set by (Re, Wo) and one of (Wi, De, lambda).
 */
#ifndef EXAMPLES_VISCOELASTICFLUIDS_ARTERIALLESION2D_VISCOUS_LOGARITHMIC_IMPLICIT_H
#define EXAMPLES_VISCOELASTICFLUIDS_ARTERIALLESION2D_VISCOUS_LOGARITHMIC_IMPLICIT_H

#include <array>
#include <complex>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/MPI.h>
#include <Rodin/PETSc.h>
#include <Rodin/Solver.h>
#include <Rodin/Types.h>
#include <Rodin/Variational.h>

#include "OrthogonalProjection.h"

namespace Rodin::Examples::ViscoelasticFluids
{
  /**
   * @brief Fully developed pulsatile flow of a Maxwell/Oldroyd-B fluid in a
   *        channel of half-width h, at a prescribed flow-rate waveform.
   *
   * @details The normalised flow rate is q(t) = sum_n c_n e^{i n omega t},
   *          c_0 = 1, c_{-n} = conj(c_n), truncated at N modes. Each mode has
   *          the profile of unit mean
   *
   *            phi_0(y) = (3/2)(1 - y^2/h^2),
   *            phi_n(y) = [1 - cosh(k_n y)/cosh(k_n h)]
   *                       / [1 - tanh(k_n h)/(k_n h)],
   *
   *          k_n^2 = i n omega rho / eta*(n omega), so that
   *          U(y, t)/Ubar = sum_n c_n phi_n(y) e^{i n omega t} and the flux
   *          of U is exactly Ubar 2h q(t).
   */
  class PulsatileInflow
  {
    public:
      using Real = Rodin::Real;
      using Complex = std::complex<Real>;

      PulsatileInflow() = default;

      /// @brief q(t) = 1 + A sin(omega t).
      PulsatileInflow& setSinusoidal(Real amplitude);

      /// @brief q(t) from a two-column file (t, Q), normalised by its mean;
      ///        only the shape is kept, the time column is rescaled to the
      ///        period of the run.
      PulsatileInflow& load(const std::string& path, size_t harmonics);

      /// @brief Fixes the fluid and the channel; computes every k_n.
      PulsatileInflow& setFluid(Real rho, Real etaS, Real etaP, Real lambda,
        Real omega, Real halfWidth);

      /// @brief q(t).
      Real flowRate(Real t) const;

      /// @brief U(y, t)/Ubar.
      Real velocity(Real y, Real t) const;

      /// @brief max_t q(t), sampled.
      Real getPeak() const;

      size_t getHarmonicCount() const noexcept
      {
        return m_c.size() - 1;
      }

    private:
      /// @brief phi_n(y), n >= 1.
      Complex mode(size_t n, Real y) const;

      std::vector<Complex> m_c{ Complex(1.0) };
      std::vector<Complex> m_k;
      Real m_omega = 0.0;
      Real m_h = 1.0;
  };

  class ArterialLesion2DViscousLogImplicit
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

      /// @brief (u, p, psi; v, q, chi). psi lives in a VectorFESType of
      ///        dimension d(d+1)/2 (Voigt), so it shares the vector types.
      using FlowProblemType =
        Rodin::Variational::Problem<LinearSystemType, VectorTrialFunctionType,
          ScalarTrialFunctionType, VectorTrialFunctionType, VectorTestFunctionType,
          ScalarTestFunctionType, VectorTestFunctionType>;

      using VectorProblemType = Rodin::Variational::Problem<LinearSystemType,
        VectorTrialFunctionType, VectorTestFunctionType>;

      using FluxFormType = Rodin::Variational::LinearForm<ScalarFESType, ::Vec>;

      /// @brief Cellwise (P0) stabilisation parameters.
      using CellFESType = Rodin::Variational::P0<Real, MeshType>;
      using CellGridFunctionType = Rodin::PETSc::Variational::GridFunction<CellFESType>;

      /// @brief sPTT/Oldroyd-B blood, as in the parent example.
      /// @details Elhanafy, Guaily & Elsaid, Adv. Mech. Eng. 11 (2019), Table 1,
      ///          after Bodnar et al.: eta_0 = 3.59e-3 Pa s, beta = 0.889.
      ///          lambda is overridden when Config::weissenberg or
      ///          Config::deborah is set.
      struct OldroydB
      {
          Real etaS = 3.19e-3;  ///< solvent viscosity eta_s (Pa s)
          Real etaP = 4.0e-4;   ///< polymeric viscosity eta_p (Pa s)
          Real lambda = 0.06;   ///< relaxation time lambda (s)
          /// @brief lambda0 = max(k lambda, lambda0_min).
          Real lambda0Factor = 1.0;
          Real lambda0Min = 0.0;
          /// @brief Linear sPTT extensibility parameter; 0 recovers Oldroyd-B.
          Real pttEpsilon = 0.1;
      };

      /**
       * @brief Boundary labels written by make_lesion_mesh.py.
       *
       * @details 1 inlet (x = -Lu D), 2 outlet (x = Ld D), 3 parent-vessel
       *          wall, 5 lesion wall (|x| <= ell). Both walls are no-slip; the
       *          lesion is tagged apart so its shear can be integrated alone.
       */
      struct Labels
      {
          Attribute inlet = 1;
          Attribute outlet = 2;
          std::array<Attribute, 2> wall{{3, 5}};
          Attribute lesion = 5;
      };

      struct Config
      {
          /// @brief Planar MEDIT "Dimension 2" mesh in units of D.
          std::string meshPath =
            "../resources/examples/viscoelastic_fluids/S75_planar_coarse.mesh";

          std::string xdmfBasename = "ArterialLesion2D_viscous_logarithmic_implicit";
          std::string csvPath = "ArterialLesion2D_viscous_logarithmic_implicit.csv";

          Labels labels;

          /// @brief Parent-vessel width D (m); the mesh is scaled by it.
          Real diameter = 4.0e-3;

          Real rho = 1060.0;
          OldroydB oldroydB;
          /// @brief Pressure penalty of the equal-order pair.
          Real pressurePenalty = 1.0e-12;

          /// @brief Re = rho Ubar D/eta_0; sets the mean inlet velocity.
          Real reynolds = 300.0;
          /// @brief Wo = (D/2) sqrt(omega rho/eta_0); sets the period.
          Real womersley = 4.0;
          /// @brief q(t) = 1 + A sin(omega t); A = 0 is steady inflow.
          Real amplitude = 0.5;
          /// @brief Optional flow-rate shape (t, Q); replaces the sinusoid.
          std::string flowWaveformPath;
          /// @brief Fourier modes kept from the file.
          int harmonics = 10;

          /// @brief Wi = lambda Ubar/D. If > 0, it sets lambda.
          Real weissenberg = 0.0;
          /// @brief De = lambda/T. If > 0, it sets lambda (and wins over Wi).
          Real deborah = 0.0;

          /// @brief The inflow is ramped over this many cycles,
          ///        r(t) = (1 - cos(pi t/T_r))/2, to start from rest without an
          ///        impulse. 0 disables it.
          Real rampCycles = 1.0;

          /// @brief Outlet traction level (Pa). Does not affect u.
          Real outletPressure = 0.0;
          Real outletBackflowStabilization = 1.0;

          int stepsPerCycle = 400;
          int cycles = 5;

          /// @brief XDMF output period, in steps. Every cycle boundary is
          ///        written regardless, so 0 keeps only those.
          int outputEvery = 20;

          /// @brief Per-term stabilisation scales, as in the parent example.
          Real vmsScale = 1.0;       ///< S1, dynamic convection
          Real gradDivScale = 1.0;   ///< S2, alpha2
          Real pressureScale = 1.0;  ///< S1, alpha1 <P^perp[grad p], grad q>
          Real stressDivScale = 1.0; ///< S1, alpha1 <P^perp[div sigma], div chi>
          Real stressScale = 1.0;    ///< S3, alpha3
          bool useVMS = true;

          /// @brief Newton on the coupled system, as in the parent example.
          int conformationIterations = 6;
          Real conformationTolerance = 1.0e-6;
          Real newtonMaxStep = 2.0;
          /// @brief Impose psi = 0 (relaxed fluid) on the inlet.
          bool inletConformation = true;

          /// @brief Divergence guard: max|u| above this multiple of Ubar
          ///        stops the run.
          Real maxVelocityFactor = 20.0;
      };

      /// @brief Wall-clock cost of each phase of one step.
      struct Timing
      {
          Real vms = 0.0;
          Real assembly = 0.0;
          Real solve = 0.0;
          Real shear = 0.0;
          Real fluxes = 0.0;
          Real output = 0.0;
          Real total = 0.0;
      };

      ArterialLesion2DViscousLogImplicit(const Rodin::Context::MPI& context, const Config& cfg);
      ~ArterialLesion2DViscousLogImplicit();

      ArterialLesion2DViscousLogImplicit(const ArterialLesion2DViscousLogImplicit&) = delete;
      ArterialLesion2DViscousLogImplicit& operator=(const ArterialLesion2DViscousLogImplicit&) = delete;

      ArterialLesion2DViscousLogImplicit& initialize();
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
      static AttributeSet makeWallSet(const Config& cfg);

      bool isRoot() const;

      /// @brief Ubar, T and lambda from (Re, Wo, Wi or De); dt = T/steps.
      void deriveParameters();

      void setupSpaces();
      void setupFlow();
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
      void computeWallShear();
      void closeCycle(Real elapsed);
      void computeFluxes();
      Real boundaryMeasure(const AttributeSet& tags);

      /// @brief r(t) = (1 - cos(pi t/T_r))/2 for t < T_r, then 1.
      Real ramp(Real t) const;

      /// @brief r(t) q(t): the ramped normalised flow rate.
      Real inflowFactor(Real t) const;

      /// @brief h_K, the local element size.
      static Real cellSize(const Rodin::Geometry::Point& p);

      /// @brief Codina tau_1 on eta_s, as in the parent example.
      Real tau1At(const Rodin::Geometry::Point& p) const;
      /// @brief alpha3, Eq. (42), with the PTT factor at psi^n.
      Real alpha3At(const Rodin::Geometry::Point& p) const;

      /// @brief The alphas (P0) and every Pi_h, all at step n.
      void updateStabilization();

      /// @brief lambda0 = max(k lambda, lambda0_min).
      Real lambda0() const;

      /// @brief PTT factor f = 1 + epsilon (lambda/lambda0)(tr exp(psi) - d).
      Real pttFactor(Real traceExp) const;

      /// @brief sigma at the nodes, from m_psiIt.
      void updateConformation();

      void writeCSVHeader();
      void writeCSVRow(int cycle);

      Config m_cfg;

      PulsatileInflow m_inflow;
      /// @brief Mean inlet velocity Ubar (m/s) and period T (s).
      Real m_meanVelocity = 0.0;
      Real m_period = 0.0;
      Real m_dt = 0.0;

      MeshType m_mesh;
      Rodin::IO::XDMF m_xdmf;
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
      /// @brief sigma = (eta_p/lambda0)(exp(psi) - I) at the nodes (Voigt).
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
      Heart::OrthogonalProjection<VectorFESType> m_piConv;     ///< Pi[(grad u^n) u^n]
      Heart::OrthogonalProjection<VectorFESType> m_sub;        ///< dynamic subscale u'
      VectorGridFunctionType m_subOld;
      Heart::OrthogonalProjection<ScalarFESType> m_piDiv;      ///< Pi[div u^n]
      Heart::OrthogonalProjection<VectorFESType> m_piGradP;    ///< Pi[grad p^n]
      Heart::OrthogonalProjection<VectorFESType> m_piDivSigma; ///< Pi[div sigma^n]
      Heart::OrthogonalProjection<VectorFESType> m_piEps;      ///< Pi[eps(u^n)], Voigt
      Heart::OrthogonalProjection<VectorFESType> m_piAdvPsi;   ///< Pi[(u^n.grad) psi^n], Voigt

      // ---- Wall shear indices ----------------------------------------------
      /// @brief The two rows of 2 eps(u), recovered onto the nodes.
      VectorGridFunctionType m_wss;
      VectorGridFunctionType m_symRec0;
      VectorGridFunctionType m_symRec1;
      VectorGridFunctionType m_netShear;
      ScalarGridFunctionType m_absShear;
      ScalarGridFunctionType m_shearMagnitude;
      ScalarGridFunctionType m_tawss;
      ScalarGridFunctionType m_osi;

      // ---- Fluxes ----------------------------------------------------------
      ScalarTestFunctionType m_qFlux;
      ScalarGridFunctionType m_one;
      FluxFormType m_flux;

      // ---- Problems and solvers --------------------------------------------
      FlowProblemType m_flow;
      Rodin::Solver::KSP m_flowKSP;
      VectorProblemType m_vectorProjection;
      Rodin::Solver::KSP m_vectorProjectionKSP;
      VectorTrialFunctionType m_wssTrial;
      VectorTestFunctionType m_wssTest;
      VectorProblemType m_wssProjection;
      Rodin::Solver::KSP m_wssKSP;

      // ---- Time-dependent coefficients read by the frozen forms ------------
      Real m_t = 0.0;
      /// @brief r(t) q(t) at the current step, read by the inlet profile.
      Real m_inflowNow = 0.0;
      Real m_inletMeasure = 1.0;
      Real m_outletMeasure = 1.0;
      Real m_qIn = 0.0;
      Real m_qOut = 0.0;
      Real m_inletPressure = 0.0;
      Real m_outletMeanPressure = 0.0;
      Real m_speed = 0.0;
      /// @brief max |sigma_ij| (Pa), reported with the flow.
      Real m_stress = 0.0;
      /// @brief Newton iterations of the last step and their last increment.
      int m_conformationIts = 0;
      Real m_psiIncrement = 0.0;
      bool m_flowFieldSplitsSet = false;
      bool m_initialized = false;

      int m_step = 0;
      Timing m_timing;

      std::ofstream m_csv;
  };
}

#endif
