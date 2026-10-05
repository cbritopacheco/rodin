/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ArterialLesionAxi_viscous_logarithmic_implicit.h
 * @brief Pulsatile sPTT/Oldroyd-B flow through an idealised stenosis or
 *        fusiform aneurysm of a pipe, axisymmetric (meridian half-plane); the
 *        formulation of ArterialLesion2D_viscous_logarithmic_implicit with the
 *        changes the cylindrical geometry requires, and no others.
 *
 * Coordinates. x is the axis, y = r >= 0 the radius; there is no swirl,
 * u = u_x e_x + u_r e_r. The domain is the meridian half-plane of the pipe of
 * ArterialLesion3D (lesion_2D.geo with planar = 0): r(x)/R = 1 + a f(x),
 * a = sqrt(1 - S) - 1 (S the AREA reduction) or a = Gam - 1; tag 4 is the
 * axis.
 *
 * What changes with respect to the planar driver:
 *  - every volume and boundary integral carries the weight r (the 2 pi is
 *    dropped throughout, and restored only in the reported flow rates);
 *  - the velocity gradient is L = blockdiag(grad u, u_r/r), so
 *    div u = div_P u + u_r/r and eps(u) has the hoop entry eps_tt = u_r/r;
 *  - psi = log(conformation) is block-diagonal, psi = blockdiag(psi_P,
 *    psi_t): the in-plane block psi_P (Voigt xx, xr, rr) is the planar field,
 *    psi_t = psi_thetatheta is a new scalar unknown. exp, Dexp and Dexp^{-1}
 *    act block by block, so the psi_P equation is the planar one, except
 *    that the PTT factor f = 1 + epsilon (lambda/lambda0)(tr exp(psi_P) +
 *    exp(psi_t) - 3) couples it to psi_t; the psi_t equation is
 *      d_t psi_t + u . grad psi_t + a(psi_t) u_r/r + g(psi_P, psi_t) = 0,
 *      a = -2 + c e^{-psi_t},  g = (f/lambda)(1 - e^{-psi_t}),
 *    the hoop entry of N = -Dexp^{-1}[L T + T L^T - c eps] + (f/lambda)(I -
 *    exp(-psi)). Without swirl the advection of a tensor has no Christoffel
 *    terms, so (u . grad) acts componentwise as in the plane;
 *  - the momentum equation gains 2 eta_s (u_r/r)(v_r/r) and sigma_tt v_r/r,
 *    sigma_tt = s (e^{psi_t} - 1), and the pressure pairs with div v;
 *  - the axis carries u_r = 0, imposed by a penalty (Rodin has no
 *    component-wise value BC); psi needs no condition there (u . n = 0);
 *  - the inlet profile is the Womersley pipe flow of the 3D driver;
 *  - fluxes are 2 pi int u . n r ds, mean pressures int p r ds / int r ds.
 *
 * After multiplication by r the hoop terms read u_r v_r / r (bounded: u_r = 0
 * on the axis), u_r q, p v_r, s e^{psi_t} v_r and a u_r chi_t: only the
 * viscous one keeps a 1/r, and it is evaluated at interior quadrature points.
 *
 * The split-OSS stabilisation keeps the planar operators (weighted by r): it
 * is weakly consistent whatever the residual operator, so the hoop terms are
 * left out of it, except for the psi_t advection, which gets its own S3 term.
 *
 * Dimensionless groups, with eta_0 = eta_s + eta_p (as in the 3D driver):
 *   Re = rho Ubar D/eta_0,  Wo = (D/2) sqrt(omega rho/eta_0),
 *   Wi = lambda Ubar/D,     De = lambda/T.
 */
#ifndef EXAMPLES_VISCOELASTICFLUIDS_ARTERIALLESIONAXI_VISCOUS_LOGARITHMIC_IMPLICIT_H
#define EXAMPLES_VISCOELASTICFLUIDS_ARTERIALLESIONAXI_VISCOUS_LOGARITHMIC_IMPLICIT_H

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
   *        pipe of radius R, at a prescribed flow-rate waveform (copied from
   *        ArterialLesion3D_viscous_logarithmic_implicit).
   *
   * @details q(t) = sum_n c_n e^{i n omega t}, c_0 = 1, c_{-n} = conj(c_n).
   *          Each mode has the profile of unit cross-sectional mean
   *
   *            phi_0(r) = 2 (1 - r^2/R^2),
   *            phi_n(r) = [1 - I_0(k_n r)/I_0(k_n R)]
   *                       / [1 - 2 I_1(k_n R)/(k_n R I_0(k_n R))],
   *
   *          k_n^2 = i n omega rho / eta*(n omega).
   */
  class PulsatilePipeInflow
  {
    public:
      using Real = Rodin::Real;
      using Complex = std::complex<Real>;

      PulsatilePipeInflow() = default;

      /// @brief q(t) = 1 + A sin(omega t).
      PulsatilePipeInflow& setSinusoidal(Real amplitude);

      /// @brief q(t) from a two-column file (t, Q), normalised by its mean.
      PulsatilePipeInflow& load(const std::string& path, size_t harmonics);

      /// @brief Fixes the fluid and the pipe; computes every k_n.
      PulsatilePipeInflow& setFluid(Real rho, Real etaS, Real etaP, Real lambda,
        Real omega, Real radius);

      /// @brief q(t).
      Real flowRate(Real t) const;

      /// @brief U(r, t)/Ubar.
      Real velocity(Real r, Real t) const;

      /// @brief max_t q(t), sampled.
      Real getPeak() const;

      /// @brief e^{-z} I_nu(z), nu = 0, 1, Re z >= 0.
      static Complex besselIScaled(int nu, Complex z);

    private:
      /// @brief phi_n(r), n >= 1.
      Complex mode(size_t n, Real r) const;

      std::vector<Complex> m_c{ Complex(1.0) };
      std::vector<Complex> m_k;
      Real m_omega = 0.0;
      Real m_radius = 1.0;
  };

  class ArterialLesionAxiViscousLogImplicit
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

      /// @brief (u, p, psi_P, psi_t; v, q, chi_P, chi_t).
      using FlowProblemType =
        Rodin::Variational::Problem<LinearSystemType,
          VectorTrialFunctionType, ScalarTrialFunctionType,
          VectorTrialFunctionType, ScalarTrialFunctionType,
          VectorTestFunctionType, ScalarTestFunctionType,
          VectorTestFunctionType, ScalarTestFunctionType>;

      using VectorProblemType = Rodin::Variational::Problem<LinearSystemType,
        VectorTrialFunctionType, VectorTestFunctionType>;

      using FluxFormType = Rodin::Variational::LinearForm<ScalarFESType, ::Vec>;

      using CellFESType = Rodin::Variational::P0<Real, MeshType>;
      using CellGridFunctionType = Rodin::PETSc::Variational::GridFunction<CellFESType>;

      /// @brief sPTT/Oldroyd-B blood, as in the planar driver.
      struct OldroydB
      {
          Real etaS = 3.19e-3;  ///< solvent viscosity eta_s (Pa s)
          Real etaP = 4.0e-4;   ///< polymeric viscosity eta_p (Pa s)
          Real lambda = 0.06;   ///< relaxation time lambda (s)
          Real lambda0Factor = 1.0;
          Real lambda0Min = 0.0;
          /// @brief Linear sPTT extensibility parameter; 0 recovers Oldroyd-B.
          Real pttEpsilon = 0.1;
      };

      /**
       * @brief Boundary labels written by make_lesion_mesh.py --axi.
       *
       * @details 1 inlet, 2 outlet, 3 parent-vessel wall, 4 axis, 5 lesion
       *          wall (|x| <= ell).
       */
      struct Labels
      {
          Attribute inlet = 1;
          Attribute outlet = 2;
          std::array<Attribute, 2> wall{{3, 5}};
          Attribute axis = 4;
          Attribute lesion = 5;
      };

      struct Config
      {
          /// @brief Axisymmetric MEDIT "Dimension 2" mesh in units of D, y >= 0.
          std::string meshPath =
            "../resources/examples/viscoelastic_fluids/S75_axi_medium.mesh";

          std::string xdmfBasename = "ArterialLesionAxi_viscous_logarithmic_implicit";
          std::string csvPath = "ArterialLesionAxi_viscous_logarithmic_implicit.csv";

          Labels labels;

          /// @brief Parent-vessel diameter D (m); the mesh is scaled by it.
          Real diameter = 4.0e-3;

          Real rho = 1060.0;
          OldroydB oldroydB;
          Real pressurePenalty = 1.0e-12;

          /// @brief Re = rho Ubar D/eta_0.
          Real reynolds = 300.0;
          /// @brief Wo = (D/2) sqrt(omega rho/eta_0).
          Real womersley = 4.0;
          /// @brief q(t) = 1 + A sin(omega t); A = 0 is steady inflow.
          Real amplitude = 0.5;
          std::string flowWaveformPath;
          int harmonics = 10;

          /// @brief Wi = lambda Ubar/D. If > 0, it sets lambda.
          Real weissenberg = 0.0;
          /// @brief De = lambda/T. If > 0, it sets lambda (and wins over Wi).
          Real deborah = 0.0;

          Real rampCycles = 1.0;

          Real outletPressure = 0.0;
          Real outletBackflowStabilization = 1.0;

          /// @brief u_r = 0 on the axis: penalty gamma int u_r v_r ds with
          ///        gamma = axisPenalty eta_0/D.
          Real axisPenalty = 1.0e8;

          int stepsPerCycle = 400;
          int cycles = 5;

          /// @brief XDMF output period, in steps; every cycle boundary is
          ///        written regardless, and 0 keeps only those. The CSV is
          ///        written every step.
          int outputEvery = 20;

          Real vmsScale = 1.0;
          Real gradDivScale = 1.0;
          Real pressureScale = 1.0;
          Real stressDivScale = 1.0;
          Real stressScale = 1.0;
          bool useVMS = true;

          int conformationIterations = 6;
          Real conformationTolerance = 1.0e-6;
          Real newtonMaxStep = 2.0;
          bool inletConformation = true;

          Real maxVelocityFactor = 20.0;
      };

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

      ArterialLesionAxiViscousLogImplicit(const Rodin::Context::MPI& context, const Config& cfg);
      ~ArterialLesionAxiViscousLogImplicit();

      ArterialLesionAxiViscousLogImplicit(const ArterialLesionAxiViscousLogImplicit&) = delete;
      ArterialLesionAxiViscousLogImplicit& operator=(const ArterialLesionAxiViscousLogImplicit&) = delete;

      ArterialLesionAxiViscousLogImplicit& initialize();
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

      void deriveParameters();

      void setupSpaces();
      void setupFlow();
      void setupWallShear();

      template <class Expression>
      void projectVector(const Expression& expr, VectorGridFunctionType& out)
      {
        m_vectorProjection = Rodin::Variational::Integral(m_wTrial, m_wTest) -
          Rodin::Variational::Integral(expr, m_wTest);
        m_vectorProjection.solve(m_vectorProjectionKSP);
        out.setData(m_wTrial.getSolution().getData());
      }

      static void axpy(Real a, const ::Vec& x, ::Vec& y);

      bool solveFlow();
      void computeWallShear();
      void closeCycle(Real elapsed);
      void computeFluxes();
      /// @brief int r ds over the tagged facets.
      Real boundaryMeasure(const AttributeSet& tags);

      Real ramp(Real t) const;
      Real inflowFactor(Real t) const;

      static Real cellSize(const Rodin::Geometry::Point& p);

      Real tau1At(const Rodin::Geometry::Point& p) const;
      Real alpha3At(const Rodin::Geometry::Point& p) const;

      void updateStabilization();

      Real lambda0() const;

      /// @brief PTT factor f = 1 + epsilon (lambda/lambda0)(tr exp(psi) - 3),
      ///        the trace over the full 3x3 tensor.
      Real pttFactor(Real traceExp) const;

      /// @brief sigma_P and sigma_t at the nodes, from (psi_P^k, psi_t^k).
      void updateConformation();

      void writeCSVHeader();
      void writeCSVRow(int cycle);

      Config m_cfg;

      PulsatilePipeInflow m_inflow;
      Real m_meanVelocity = 0.0;
      Real m_period = 0.0;
      Real m_dt = 0.0;

      MeshType m_mesh;
      Rodin::IO::XDMF m_xdmf;
      AttributeSet m_wallSet;

      VectorFESType m_vh;
      ScalarFESType m_sh;
      /// @brief In-plane tensor space: P1, Voigt order (xx, xr, rr).
      VectorFESType m_tauh;
      CellFESType m_ch;

      /// @brief r = y, exact in P1: the weight of every integral.
      ScalarGridFunctionType m_radius;

      // ---- Flow ------------------------------------------------------------
      VectorTrialFunctionType m_u;
      ScalarTrialFunctionType m_p;
      VectorTestFunctionType m_v;
      ScalarTestFunctionType m_q;
      VectorGridFunctionType m_uOld;
      VectorGridFunctionType m_uIt;
      /// @brief psi_P (Voigt), its test function, psi_P^n and psi_P^k.
      VectorTrialFunctionType m_psi;
      VectorTestFunctionType m_chi;
      VectorGridFunctionType m_psiOld;
      VectorGridFunctionType m_psiIt;
      /// @brief psi_t = psi_thetatheta, its test function, psi_t^n and psi_t^k.
      ScalarTrialFunctionType m_psiT;
      ScalarTestFunctionType m_chiT;
      ScalarGridFunctionType m_psiTOld;
      ScalarGridFunctionType m_psiTIt;
      std::uint64_t m_psiRevision = 0;
      /// @brief sigma_P (Voigt) and sigma_t at the nodes.
      VectorGridFunctionType m_sigma;
      ScalarGridFunctionType m_sigmaT;

      VectorTrialFunctionType m_wTrial;
      VectorTestFunctionType m_wTest;

      // ---- VMS, split OSS --------------------------------------------------
      CellGridFunctionType m_tauK;
      CellGridFunctionType m_alpha1;
      CellGridFunctionType m_alpha2;
      CellGridFunctionType m_alpha3;
      CellGridFunctionType m_alphaPsi;
      Heart::OrthogonalProjection<VectorFESType> m_piConv;
      Heart::OrthogonalProjection<VectorFESType> m_sub;
      VectorGridFunctionType m_subOld;
      Heart::OrthogonalProjection<ScalarFESType> m_piDiv;
      Heart::OrthogonalProjection<VectorFESType> m_piGradP;
      Heart::OrthogonalProjection<VectorFESType> m_piDivSigma;
      Heart::OrthogonalProjection<VectorFESType> m_piEps;
      Heart::OrthogonalProjection<VectorFESType> m_piAdvPsi;
      /// @brief Pi[u^n . grad psi_t^n].
      Heart::OrthogonalProjection<ScalarFESType> m_piAdvPsiT;

      // ---- Wall shear indices ----------------------------------------------
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

      Real m_t = 0.0;
      Real m_inflowNow = 0.0;
      /// @brief int r ds over the inlet and the outlet (m^2): R^2/2 each.
      Real m_inletMeasure = 1.0;
      Real m_outletMeasure = 1.0;
      /// @brief Volumetric flow rates (m^3/s).
      Real m_qIn = 0.0;
      Real m_qOut = 0.0;
      Real m_inletPressure = 0.0;
      Real m_outletMeanPressure = 0.0;
      Real m_speed = 0.0;
      Real m_stress = 0.0;
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
