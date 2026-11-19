/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file AtriumViscoelastic.h
 * @brief Newtonian, generalised Newtonian, viscoelastic and thixotropic
 *        viscoelastic blood in the fluid region of the Atrium mesh, on a fixed
 *        domain: the log-conformation solver of
 *        examples/viscoelastic_fluids/ArterialLesion3D_viscous_logarithmic_implicit
 *        with the openings of the Atrium fluid.
 *
 * Formulation. Unchanged from the arterial driver: T = 2 eta_s eps(u) + sigma,
 * sigma = (eta_p/lambda0)(exp(psi) - I), (u, p, psi) by Newton on P1 spaces,
 * BDF1 or BDF2 with step control, term-by-term orthogonal subscales with
 * dynamic velocity subscales on momentum and continuity, and the
 * residual-based orthogonal subscale of the psi equation ("psi_rb"; "exp" and
 * "psi" are kept). The constitutive equation carries no artificial stress
 * diffusion: its stabilisation acts on the orthogonal projection of its
 * residual. README.md of this folder states what was checked.
 *
 * Domain. The atrial volume of the PTAG atrium (Coarse.mesh), or the cells of
 * material Config::Labels::fluidMaterial of a combined Atrium mesh (fluid 1,
 * solid 4), extracted on each rank from the partitioned parent as in
 * Atrium/Fluid. A mesh with no cell of that material is used whole.
 * Coordinates are scaled by Config::meshScale (mm to m).
 *
 * Boundary conditions, those of Atrium/Fluid on a fixed domain:
 *  - Pulmonary veins (labels 1, 2, 3, 5): traction -p_pv(t) n, p_pv = mean +
 *    amplitude sin(2 pi t/T), with the normal impedance and the tangential
 *    damping of Atrium/Fluid and the backflow term on the entering flux.
 *  - Mitral plane (label 6): traction -(p_out + Z u.n) n. With the lumped
 *    ventricle (Rodin::Heart::CCMLC2014) p_out = p_v and Z = |Gamma|/K, K the
 *    conductance of the valve law of that model, open while the mean pressure
 *    of the plane exceeds p_v; with feedback the ventricle reads that mean
 *    pressure as its atrial pressure. Without the ventricle, a prescribed
 *    pressure and resistance.
 *  - Wall (label 4): no slip.
 *  - psi, the aggregate size n and the age are advected fields: on every
 *    opening the upwind term |u.n|_- (x - x_in) imposes the entering value
 *    weakly, x_in = 0 (relaxed fluid), n_in and 0 on the veins and the value
 *    of the previous step on the mitral plane. Nothing is imposed where the
 *    flow leaves.
 *
 * Fluids (Config::fluid). "cy", "owens" with Owens::equilibrium and "owens"
 * share the steady shear viscosity eta_CY(gd) exactly, so that a difference
 * between them is one of memory and not of the flow curve:
 *  - "newtonian"  eta_ref (eta_inf by default), two-field (u, p);
 *  - "cy"         Carreau-Yasuda, two-field (u, p);
 *  - "owens" with equilibrium: the aggregate size follows the local shear
 *                 rate, n = n_st(gd), so lambda = (eta_CY(gd) - eta_plasma)/G:
 *                 elasticity without structural memory (a White-Metzner fluid);
 *  - "owens"      Owens (2006): n is transported and relaxes at the rates
 *                 a(gd), b(gd): elasticity and thixotropy;
 *  - "gsptt", "sptt": the single-mode sPTT fluids of the arterial driver.
 */
#ifndef RODIN_ATRIUM_VISCOELASTIC_ATRIUMVISCOELASTIC_H
#define RODIN_ATRIUM_VISCOELASTIC_ATRIUMVISCOELASTIC_H

#include <limits>
#include <array>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

#include <Rodin/Geometry.h>
#include <Rodin/Heart/CCMLC2014.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/MPI.h>
#include <Rodin/PETSc.h>
#include <Rodin/Solver.h>
#include <Rodin/Types.h>
#include <Rodin/Variational.h>

#include "OrthogonalProjection.h"

namespace Rodin::Examples::ViscoelasticFluids
{
  class AtriumViscoelastic
  {
    public:
      using Real = Rodin::Real;
      using Attribute = Rodin::Geometry::Attribute;
      using MeshType = Rodin::Geometry::Mesh<Rodin::Context::MPI>;

      /// @brief The fluid region, a submesh of the partitioned parent.
      using FluidMeshType = Rodin::Geometry::SubMesh<Rodin::Context::MPI>;

      /// @brief Lumped left ventricle with its mitral valve law.
      using VentricleType = Rodin::Heart::CCMLC2014T<>;
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

      /**
       * @brief (u, p, psi_D, psi_O; v, q, chi_D, chi_O). Each block of psi is
       *        a three-component P1 field, so it shares the vector types.
       */
      using FlowProblemType = Rodin::Variational::Problem<LinearSystemType,
        VectorTrialFunctionType, ScalarTrialFunctionType, VectorTrialFunctionType,
        VectorTrialFunctionType, VectorTestFunctionType, ScalarTestFunctionType,
        VectorTestFunctionType, VectorTestFunctionType>;

      using VectorProblemType = Rodin::Variational::Problem<LinearSystemType,
        VectorTrialFunctionType, VectorTestFunctionType>;

      using ScalarProblemType = Rodin::Variational::Problem<LinearSystemType,
        ScalarTrialFunctionType, ScalarTestFunctionType>;

      using FluxFormType = Rodin::Variational::LinearForm<ScalarFESType, ::Vec>;

      /// @brief Cellwise (P0) stabilisation parameters.
      using CellFESType = Rodin::Variational::P0<Real, MeshType>;
      using CellGridFunctionType = Rodin::PETSc::Variational::GridFunction<CellFESType>;

      /// @brief sPTT/Oldroyd-B blood, as in the 2D driver.
      struct OldroydB
      {
        Real etaS = 3.19e-3;  ///< solvent viscosity eta_s (Pa s)
        Real etaP = 4.0e-4;   ///< polymeric viscosity eta_p (Pa s)
        Real lambda = 0.1;   ///< relaxation time lambda (s)

        /// @brief lambda0 = max(k lambda, lambda0_min).
        Real lambda0Factor = 1.0;
        Real lambda0Min = 0.0;

        /// @brief Linear sPTT extensibility parameter; 0 recovers Oldroyd-B.
        Real pttEpsilon = 0.1;
      };

      /**
       * @brief Labels of the Atrium mesh (Atrium/Fluid): pulmonary veins 1, 2,
       *        3, 5; mitral plane 6; wall (fluid-solid interface) 4; fluid
       *        cells 1.
       */
      struct Labels
      {
        std::array<Attribute, 4> inlets{{1, 2, 3, 5}};
        Attribute outlet = 6;
        std::array<Attribute, 1> wall{{4}};
        Attribute fluidMaterial = 1;
      };

      /// @brief Piecewise-linear activation of the lumped ventricle (Atrium/Fluid).
      struct Activation
      {
        Real tRampStart = 0.15;
        Real tRampEnd = 0.21;
        Real tPlateauEnd = 0.36;
        Real tRelaxEnd = 0.45;
        Real tNegativeEnd = 0.6;
        Real positiveValue = 35.0;
        Real negativeValue = -20.0;
      };

      /**
       * @brief Atrial pressure driving the ventricle when it does not read the
       *        mitral plane (Atrium/Fluid).
       */
      struct AtrialPressure
      {
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

      /**
       * @brief Lumped left ventricle, the parameters and initial state of
       *        Atrium/Fluid.
       */
      struct Ventricle
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

      /**
       * @brief Ball that marks the left atrial appendage, in metres (the
       *        coordinates of the XDMF output). radius <= 0: no appendage
       *        diagnostics.
       */
      struct Appendage
      {
        std::array<Real, 3> center{{0.0, 0.0, 0.0}};
        Real radius = 0.0;
      };

      /**
       * @brief Carreau-Yasuda viscosity of whole blood, the constants of
       *        Atrium/Fluid. Cho & Kensey (1991), used by the arterial
       *        driver: 0.056, 0.00345, 1.902, 0.22, 1.25.
       */
      struct CarreauYasuda
      {
        Real etaZero = 0.0720; ///< zero-shear viscosity (Pa s)
        Real etaInf = 0.0049; ///< infinite-shear viscosity (Pa s)
        Real lambda = 6.1529; ///< time constant (s)
        Real n = 0.2081;
        Real a = 1.4173;
      };

      /**
       * @brief Owens (2006) blood model, JNNFM 140:57-70, single mode: red
       *        cells aggregate into rouleaux of mean size n (Hookean
       *        dumbbells, friction proportional to the size). Eqs. (23), (24),
       *        (10), (27), (28) of the paper:
       *          tau + mu tau^nabla = G mu gd,  G = N0 kB T = eta_p,inf/lambda_H,
       *          mu(n, gd) = n lambda_H / (1 + (b n (n - 1)/2 + aN0) lambda_H),
       *          D n/Dt = -(1/2) b (n - n_st)(n + n_st - 1),
       *          n_st = (eta(gd) - eta_s)/eta_p,inf (1 + (3/2) aN0 lambda_H),
       *          b = aN0/(n_st (n_st - 1)),
       *        so that G mu(n_st) = eta(gd) - eta_s exactly. eta(gd) is the
       *        Carreau-Yasuda of "cy" (Owens fits a Cross model: any steady
       *        fit will do, Eq. (29)); the plasma is the solvent, as in Fang &
       *        Owens (2006) and Moyers-Gonzalez et al. (JFM 2008), and eta_p,inf
       *        = eta_CY(inf) - eta_plasma.
       *        Aggregation rate aN0(gd): Moyers-Gonzalez, Owens & Fang (JFM
       *        2008, Part I, Appendix A, Eq. A1), the C1 cubics with zero slope
       *        at 0, gd_c = 5.78 1/s (Murata & Secomb) and gd_max = 900 1/s:
       *          aN0 = a0 + (ac - a0)(3 s^2 - 2 s^3),  s = gd/gd_c,          gd <= gd_c
       *          aN0 = ac (1 - 3 s^2 + 2 s^3),  s = (gd - gd_c)/(gd_max - gd_c), gd <= gd_max
       *          aN0 = 0,                                                     gd > gd_max,
       *        with a0 = a_{1,0} = 1 1/s and ac = 2 a0 (their coefficients
       *        a_{1,k}, a_{2,k} expand exactly to this). b -> 0 with aN0 near
       *        gd_max, which would freeze n in the jet; bFloor keeps b at its
       *        maximum above the shear rate where it peaks (steady viscosity
       *        changed by < 0.2 %, printed at start).
       */
      struct Owens
      {
        Real lambdaH = 0.005; ///< single-cell relaxation time (s)
        Real etaPlasma = 1.2e-3; ///< plasma (solvent) viscosity (Pa s)
        Real a0 = 1.0; ///< aN0 at rest, a_{1,0} (1/s)
        Real aCrit = -1.0; ///< aN0 at gd_c (1/s); < 0: 2 a0 (JFM 2008)
        Real gdCrit = 5.78; ///< shear rate of the aggregation peak (1/s)
        Real gdMax = 900.0; ///< aN0 = 0 above it (1/s)
        bool bFloor = true;

        /// @brief sPTT epsilon of the dumbbells (0: Hookean, Owens' model).
        Real epsilon = 0.0;

        /**
         * @brief n = n_st(gd) at every point and time: the limit of fast
         *        kinetics, with the same steady shear viscosity and the same
         *        modulus. lambda = (eta_CY(gd) - eta_plasma)/G follows the
         *        local shear rate at once, so the fluid has elasticity and
         *        no structural memory.
         */
        bool equilibrium = false;

        /**
         * @brief Initial aggregate size: < 0 the rest value n_st(0); 0 the
         *        inlet value; > 0 that value.
         */
        Real nInitial = -1.0;

        /**
         * @brief Aggregate size of the blood that enters through the veins.
         *        The pulmonary capillaries deliver single cells and the
         *        transit to the atrium (~1 s) is short against aggregation
         *        (n_st/a ~ 20 s at rest): n ~ 1 + a t/2 ~ 1.5. <= 0: the
         *        steady value at Config::inletShearRate.
         */
        Real nInlet = 1.5;

        /**
         * @brief Sets n to the periodic state of slow kinetics at the end of
         *        this cycle (0: never). Without advection, a point whose
         *        structure relaxes over many periods holds the n of zero
         *        mean rate over the cycle, <b (n - n_st)(n + n_st - 1)> = 0,
         *        that is n (n - 1) = <a>/<b> since b n_st (n_st - 1) = a.
         *        The means are taken over that cycle. The transport then
         *        continues from it: at rest n relaxes in ~20 s, which a run
         *        of a few cycles started from any uniform n does not reach.
         */
        int resetCycle = 2;
      };

      struct Config
      {
        /**
         * @brief Discretisation of the constitutive equation and of its
         *        stabilisation (file header). In all three, momentum and
         *        continuity carry the split OSS with dynamic subscales.
         *  - "psi_rb": the psi equation with the OSS of its full residual,
         *    a log-conformation subscale psi~ seen by the psi equation
         *    (streamline) and by the momentum (sigma~ = s Dexp[psi~]).
         *  - "exp": Moreno, Codina, Baiges, Castillo (CMAME 2019) and Moreno,
         *    Codina, Baiges (CMAME 2020). The constitutive equation is the
         *    sigma equation over 2 eta_p written in E = exp(psi) and tested
         *    with a stress-like chi; momentum and continuity carry the
         *    term-by-term (split) OSS, with dynamic subscales u1
         *    (convection), u2 (pressure gradient), u3 (div sigma) and the
         *    algebraic pressure subscale; the constitutive equation carries
         *    the OSS of its full residual, tested with the adjoint operator
         *    A(chi, v), with a dynamic stress subscale. Newton is exact in
         *    (u, p, psi): Dexp and D2exp by Daleckii-Krein.
         *  - "psi": the former scheme, the psi equation d_t psi + u.grad psi
         *    - Dexp^{-1}[L T + T L^T] + (f/lambda)(I - exp(-psi)) = 0 with the
         *    term-by-term OSS of eps(u) and (u.grad) psi.
         */
        std::string formulation = "psi_rb";

        /**
         * @brief Stabilisation constants of "psi_rb" and "exp" (tau_psi =
         *        (lambda/(2 eta_p)) alpha3), linear
         *        elements (Moreno et al. 2019): alpha1 = [c1 eta_0/h1^2 + c2
         *        rho |u|/h2]^-1, alpha2 = h1^2/(c1 alpha1), alpha3 = [c3 f/(2
         *        eta_p) + c4 (lambda |u|/(2 eta_p h2) + lambda |grad u|/eta_p)]^-1,
         *        h1 = |K|^(1/3), h2 the length of K along u.
         */
        Real c1 = 4.0;
        Real c2 = 2.0;
        Real c3 = 4.0;
        Real c4 = 0.25;

        /**
         * @brief Dynamic subscales (time-dependent, BDF1 in their own
         *        equations): alpha1 -> (rho/dt + 1/alpha1)^-1, alpha3 ->
         *        (lambda/(2 eta_p dt) + 1/alpha3)^-1, plus their memory terms.
         *        Off: quasi-static subscales ("psi_rb" and "exp").
         */
        bool dynamicSubscales = true;

        /**
         * @brief "psi_rb": test the log-conformation subscale with the full
         *        adjoint of the linearised psi operator, -(u^n.grad) chi +
         *        JN^*[chi] (the reaction part JN^* = (dN/dpsi)^* at step n,
         *        as in the stress adjoint of Moreno et al.); off (default):
         *        the convective part only. JN carries the deformation terms
         *        of N, of size |grad u| and of either sign, and tau |JN| =
         *        O(1) in the stenosis throat: with it, Newton failed at t =
         *        0.046 s in the S75 test (dt cut by 8), without it no step
         *        was rejected to t = 0.115 s. In the sPTT pipe both give the
         *        same errors.
         */
        bool adjointReaction = false;

        /**
         * @brief "psi_rb" and "exp": keep d_t of the FE field (BDF) in the
         *        residual that drives the stress subscale (on, default);
         *        off: without it, as in the OSS, where P^perp[d_t psi_h] = 0,
         *        and as in the term-by-term form of Castillo, Moreno, Baiges,
         *        Codina, ACME. Off was tried as the fix of the analysis
         *        below and refuted (S75 R Re 600, xcoarse, 200 steps/cycle):
         *        on ran the ramp cycle without a rejected step, off cut dt
         *        to 2e-5 s with max|sigma| ~ 2.5e4 Pa from t = 0.157 s. "exp": s~ = alpha3' (memS s~^n - R + Pi), R ~
         *        k_a (E - E_hist)/dtM, alpha3' -> 2 eta_p dt/lambda where
         *        alpha3 >> that: the stress s (E - I) + s~ seen by the
         *        momentum responds to E^{n+1} with s (1 - dt/dtM), -s/2 in
         *        BDF2. "psi_rb": With it, psi~^{n+1} contains
         *        -(tau'/dtM) psi_h^{n+1}, tau' = (1/dt + 1/tau)^-1 (BDF1 in
         *        the subscale) and dtM = 2 dt/3 (BDF2): the stress seen by
         *        the momentum, s Dexp[psi_h + psi~], responds to psi_h^{n+1}
         *        with the factor 1 - tau'/dtM, which is -1/2 wherever tau >>
         *        dt (slow flow, flow reversal, and everywhere as dt -> 0):
         *        the u-psi coupling changes sign. The S75 Re 600 test grew
         *        at ~8000 1/s from t = 0.135 s (deceleration) whatever dt.
         */
        bool subscaleTimeResidual = true;

        /**
         * @brief "sptt", "gsptt", "cy", "newtonian" or "owens" (file
         *        header, Owens). "newtonian": constant viscosity eta_ref
         *        (default eta_inf of the Carreau-Yasuda law), so that Re and
         *        Wo, hence Ubar and T, are those of "cy" and "gsptt".
         */
        std::string fluid = "cy";

        /**
         * @brief "cy" and "newtonian" solve (u, p) only (twoField()): the
         *        same momentum, continuity and split OSS as the three-field
         *        system with sigma = 0, without the psi unknowns, their
         *        equation and their projections. Off: the former passive
         *        polymer (psi kept, same cost as "gsptt").
         */
        bool twoField = true;

        /**
         * @brief Eulerian residence time (age) of the blood, solved after
         *        each accepted step: d_t a + u.grad a = 1, a = 0 in the blood
         *        that enters through the veins, upwind on the mitral
         *        backflow, SUPG; written as "residenceTime" and its volume
         *        and appendage means in the CSV.
         */
        bool residenceTime = true;

        /**
         * @brief Finite-difference check of the Newton Jacobian at this
         *        step (>= 1; 0 is off): F(x) = A(x) x - b(x) at x^k and at
         *        x^k + e d, d the first Newton increment, against A(x^k) d,
         *        for e = 1e-1 ... 1e-7; the run stops after the check.
         */
        int jacobianCheck = 0;
        CarreauYasuda carreauYasuda;
        Owens owens;

        /**
         * @brief Viscosity defining Re and Wo (Pa s); 0 means eta_s + eta_p
         *        for "sptt" and eta_inf for "cy" and "gsptt".
         */
        Real etaRef = 0.0;

        /**
         * @brief MEDIT mesh in mm: the atrial volume of the PTAG atrium
         *        (Coarse.mesh, fluid only, used whole), or a combined Atrium
         *        mesh (fluid 1, solid 4), with the same boundary labels.
         */
        std::string meshPath = "../resources/examples/Heart/Atrium/Coarse.mesh";

        /// @brief Output paths; missing directories are created.
        std::string xdmfBasename = "results/AtriumViscoelastic";
        std::string csvPath = "results/AtriumViscoelastic.csv";

        Labels labels;
        Appendage appendage;

        /// @brief Mesh coordinate scale applied after partitioning (mm to m).
        Real meshScale = 1.0e-3;

        /**
         * @brief Scales of the reported groups and of the divergence guard,
         *        not of the problem: a length (m) and a speed (m/s).
         */
        Real referenceLength = 0.03;
        Real referenceVelocity = 0.5;

        Real rho = 1060.0;
        OldroydB oldroydB;

        /// @brief Pressure penalty of the equal-order pair.
        Real pressurePenalty = 1.0e-12;

        /// @brief Cardiac period (s).
        Real period = 0.85;

        /// @brief Pulmonary venous pressure, p = mean + amplitude sin(2 pi t/T) (Pa).
        Real inletPressureMean = 1200.0;
        Real inletPressureAmplitude = 200.0;

        /**
         * @brief Normal impedance and tangential damping of the venous
         *        openings (Pa s/m), as in Atrium/Fluid: a pressure inlet
         *        without source impedance has no energy bound.
         */
        Real inletImpedance = 1.0e3;
        Real inletTangentialDamping = 1.0e3;
        Real inletBackflowStabilization = 1.0;
        Real outletBackflowStabilization = 1.0;

        /// @brief Shear rate that sets n_in when Owens::nInlet <= 0 (1/s).
        Real inletShearRate = 100.0;

        /**
         * @brief Lumped ventricle on the mitral plane (file header). Off: the
         *        prescribed pressure (Pa) and resistance (Pa s/m^3) below.
         */
        bool useVentricularModel = true;

        /**
         * @brief The ventricle reads the mean pressure of the mitral plane
         *        and the plane applies its valve law (conductances Kat, Kp).
         *        Off: the ventricle follows AtrialPressure and the plane the
         *        resistances mitralResistanceOpen/Closed.
         */
        bool ventricularFeedback = true;
        Real mitralResistanceOpen = 5.0e6;
        Real mitralResistanceClosed = 1.0e12;
        Real prescribedOutletPressure = 0.0;
        Real prescribedOutletResistance = 0.0;
        Activation activation;
        AtrialPressure atrialPressure;
        Ventricle lv;

        /// @brief De = lambda/T. If > 0, it sets lambda ("sptt", "gsptt").
        Real deborah = 0.0;

        /**
         * @brief The venous pressure, and the outlet pressure when the
         *        ventricle does not read the mitral plane, are scaled by
         *        r(t) = (1 - cos(pi t/T_r))/2 over this many cycles. 0
         *        applies them at once.
         */
        Real rampCycles = 0.25;

        int stepsPerCycle = 850;

        /**
         * @brief Time step (s). If > 0, it sets stepsPerCycle = round(T/dt),
         *        so that a cycle holds a whole number of steps, and the run
         *        uses dt = T/stepsPerCycle.
         */
        Real timeStep = 0.0;
        int cycles = 6;

        /**
         * @brief Time integrator: 1 = BDF1 (backward Euler), 2 = BDF2.
         *        BDF2 discretises d/dt by (3 x^{n+1} - 4 x^n + x^{n-1})/(2 dt)
         *        in the momentum and psi equations, the Galerkin convection
         *        advects with u* = 2 u^n - u^{n-1}, and Newton starts from
         *        (u*, p*, psi*). The orthogonal projections of the split-OSS
         *        stabilisation are explicit; see bdf2Extrapolate. The
         *        stabilisation parameters, the streamline directions,
         *        Pi[(grad u) u] and the dynamic subscale stay at step n. The
         *        first step takes x^{-1} = x^0, exact for a start from rest
         *        with the (1 - cos) ramp.
         */
        int timeOrder = 1;

        /**
         * @brief BDF2: the projections taken at the extrapolation (u*, p*,
         *        psi*) instead of step n, as a bit mask: Pi[div u] 1,
         *        Pi[grad p] 2, Pi[div sigma] 4, Pi[eps(u)] 8,
         *        Pi[(u^n . grad) psi] 16, Pi[u^n . grad psi_t] 32
         *        (axisymmetric). A projection lagged at step n leaves an
         *        O(dt) error times its weight. The default 8 extrapolates
         *        Pi[eps(u)], whose weight alpha3 ~ eta_p/2 does not vanish with
         *        h. Extrapolating the others diverged: Pi[grad p], Pi[div u]
         *        and Pi[div sigma] (masks 2, 3, 6) in the axisymmetric
         *        pulsatile Oldroyd-B test after 2 to 3 periods, and the psi
         *        advection projection (16) with the "gsptt" blood model within
         *        the first period. Lagged, Pi[grad p] and Pi[div sigma]
         *        (weights O(h^2)) leave an O(h^2 dt) error.
         */
        int bdf2Extrapolate = 8;

        /// @brief XDMF output period, in steps; cycle boundaries always.
        int outputEvery = 20;

        /// @brief Per-term stabilisation scales, as in the 2D driver.
        Real vmsScale = 1.0; ///< S1, dynamic convection
        Real gradDivScale = 1.0; ///< S2, alpha2
        Real pressureScale = 1.0; ///< S1, alpha1 <P^perp[grad p], grad q>
        Real stressDivScale = 1.0; ///< S1, alpha1 <P^perp[div sigma], div chi>
        Real stressScale = 1.0; ///< S3, alpha3
        bool useVMS = true;

        /// @brief Newton on the coupled system.
        int conformationIterations = 6;
        Real conformationTolerance = 1.0e-6;
        Real newtonMaxStep = 2.0;

        /**
         * @brief Relaxed fluid (psi = 0) enters through the veins: the
         *        upwind term |u*.n|_- psi : chi on their backflow-free side.
         */
        bool inletConformation = true;

        /// @brief Print the Newton history at every step, not only the first three.
        bool traceNewton = false;

        /**
         * @brief Keep the D2exp part of the S1 div sigma linearisation
         *        (quadratic Newton). Without it the converged solution is the
         *        same, Newton converges linearly, and each iteration costs
         *        less (four bilinear and two linear integrals fewer).
         */
        bool newtonD2 = true;

        /**
         * @brief Quadrature degree of the integrals that read the pointwise
         *        log-conformation data (tetrahedra: 2 -> 4 points, 4 -> 11).
         */
        int psiQuadratureOrder = 4;

        /**
         * @brief Divergence guard: max|u| above this multiple of
         *        referenceVelocity rejects the step (stops the run without
         *        adaptiveTimeStep).
         */
        Real maxVelocityFactor = 20.0;

        /**
         * @brief Volume diagnostics written to the CSV at every accepted
         *        step (README.md): energies, dissipations, Weissenberg
         *        number, subscale levels.
         */
        bool metrics = true;

        /**
         * @brief Probes: "x,y,z;x,y,z;..." in metres. Each records the
         *        Gaussian-weighted mean (standard deviation probeRadius, m)
         *        of u and tr(sigma) at every accepted step, in
         *        <csv stem>_probes.csv: the series for the spectra and the
         *        cycle-to-cycle statistics.
         */
        std::string probes;
        Real probeRadius = 1.5e-3;

        /**
         * @brief Diagnostics: XDMF output at every accepted step with
         *        max|sigma| > dumpStress (Pa; 0 off), at most dumpMax times;
         *        the last accepted state is also written when the run aborts.
         */
        Real dumpStress = 0.0;
        int dumpMax = 12;

        /**
         * @brief No-slip kinematics in the psi equation. At a no-slip wall
         *        grad u = a n^T with a . n = 0 (pure shear), but a P1 velocity
         *        on a faceted curved wall gives the first element layer a
         *        spurious normal extension n n^T, which stretched psi along
         *        the wall normal at wall nodes (u = 0, nothing advects it
         *        away): psi ~ 12, |sigma| ~ 2e4 Pa on the distal flank of S75,
         *        stretch within 3-7 deg of the normal (runs 13-14). With it,
         *        the gradient read by the psi equation (and its residual) is
         *        L' = (1 - w) L + w (I - n n^T) L n n^T, w the P1 interpolant
         *        of the wall-node indicator (1 on the wall, 0 from the first
         *        interior node), n the area-weighted nodal wall normal. L' is
         *        linear in L, so Newton stays exact; the momentum is untouched.
         */
        bool wallKinematics = true;

        /**
         * @brief Step control. A step is rejected, and retried from the same
         *        x^n with dt * dtShrink, when Newton does not converge in
         *        conformationIterations, when its undamped increment grows,
         *        when the linear solver fails or when the divergence guard
         *        trips. Nothing of x^n, x^{n-1} changes before a step is
         *        accepted, so a retry starts from the same state. After
         *        dtGrowAfter consecutive steps of at most newtonEasy
         *        iterations dt grows by dtGrow, never above the nominal
         *        T/stepsPerCycle; a step of newtonHard or more iterations
         *        shrinks the next one. Steps are clipped to end on the cycle
         *        boundaries; BDF2 is then variable-step. Iterations count
         *        undamped ones only. Without it, a step
         *        whose Newton does not converge is accepted with a warning and
         *        any other failure stops the run (the former behaviour).
         */
        bool adaptiveTimeStep = true;

        /// @brief Retries of one step; the smallest dt is dt0 dtShrink^maxStepRetries.
        int maxStepRetries = 6;
        Real dtShrink = 0.5;

        /**
         * @brief Growth factor; BDF2 with variable steps is zero-stable for
         *        dt_{n+1}/dt_n < 1 + sqrt(2), so it stays well below that.
         */
        Real dtGrow = 1.5;
        int dtGrowAfter = 10;
        int newtonEasy = 2;
        int newtonHard = 5;

        /**
         * @brief Newton acceptance by the quadratic prediction: after an
         *        undamped increment d_k, the next one is ~ C d_k^2, with C
         *        estimated from the steps that took two undamped iterations or
         *        more. The iterate is accepted when C d_k^2 is below
         *        newtonPredictSafety * conformationTolerance, which saves the
         *        confirming iteration (one assembly and one solve). Every
         *        newtonCalibrateEvery steps the prediction is not used, so
         *        that C keeps being measured.
         */
        bool newtonPredict = true;
        Real newtonPredictSafety = 0.1;
        int newtonCalibrateEvery = 20;

        /**
         * @brief Lagged preconditioner of the flow solve: the preconditioner
         *        (ILU blocks or LU factors) of an earlier matrix is reused for
         *        up to pcReuseSolves solves, across Newton iterations and
         *        steps. A solve with it that needs more than pcReuseGrowth
         *        times the iterations of the first Newton solve after the
         *        last build is kept, and the preconditioner is rebuilt for the
         *        next solve; one that fails is repeated, from the same guess,
         *        with a fresh preconditioner. A change of dt rebuilds it too.
         *        1 rebuilds it at every solve. Only meaningful with a Krylov
         *        method (-ksp_type gmres); with preonly it is ignored.
         */
        int pcReuseSolves = 1;
        Real pcReuseGrowth = 2.0;

        /**
         * @brief Weak inflow condition of the psi equation on the outlet:
         *        where u*.n < 0 the term |u*.n|_- (psi - psi^n) : chi makes
         *        the entering fluid carry the conformation found there at
         *        step n. The psi equation is hyperbolic, and without it a
         *        backflow on the outlet enters with no data.
         */
        bool outletConformationInflow = true;

        /**
         * @brief Skew-symmetric psi advection, + (1/2)(div u*) psi : chi: the
         *        discrete div u is not zero for P1/P1, and this restores the
         *        L2 stability of the advection of psi (Temam's term of the
         *        momentum equation, applied to psi). Weakly consistent.
         */
        bool psiSkewAdvection = false;
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

        /**
         * @brief Rejected attempts of the step, preconditioner builds and
         *        Krylov iterations over all its attempts.
         */
        int retries = 0;
        int pcBuilds = 0;
        int kspIterations = 0;
      };

      AtriumViscoelastic(const Rodin::Context::MPI& context, const Config& cfg);
      ~AtriumViscoelastic();

      AtriumViscoelastic(const AtriumViscoelastic&) = delete;
      AtriumViscoelastic& operator=(const AtriumViscoelastic&) = delete;

      AtriumViscoelastic& initialize();
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
      static MeshType makeParentMesh(
        const Rodin::Context::MPI& context, const Config& cfg);
      static FluidMeshType makeFluidMesh(const MeshType& parent, const Config& cfg);
      static VentricleType::Input makeVentricleInput(
        const Config& cfg, const Real* atrialPressure);
      static Real activationWave(const Activation& cfg, Real period, Real t);
      static Real atrialWave(const AtrialPressure& cfg, Real period, Real t);

      /**
       * @brief Opening data at time t: venous pressure, ventricular pressure
       *        (interpolated in the interval the ventricle has covered, and
       *        advanced beyond it), valve state and mitral impedance.
       */
      bool updateBoundaryData(Real t);
      Real mitralResistance(bool open) const;

      /// @brief Config::metrics: the volume integrals of the CSV.
      void computeMetrics();

      /// @brief Integral over the fluid of a pointwise scalar.
      template <class F>
      Real integrate(F&& f);
      bool inAppendage(const Rodin::Geometry::Point& p) const;
      void setupProbes();
      void writeProbes();

      /// @brief Owens::resetCycle: accumulates <a>, <b> over the step.
      void accumulateKinetics();
      void resetStructure(Real elapsed);
      static AttributeSet makeWallSet(const Config& cfg);

      bool isRoot() const;

      /// @brief Checks the fluid, sets its constants and dt = T/steps.
      void deriveParameters();

      void setupSpaces();
      void setupFlow();

      /// @brief The "exp" formulation (Config::formulation): its forms.
      void setupFlowExp();

      /// @brief "exp": the alphas at step n (dynamic or quasi-static).
      void updateStabilizationExp();

      /**
       * @brief "exp", after a step is accepted and before the history moves:
       *        the subscales of step n+1 at the quadrature points, then
       *        Pi[R^{n+1}], the projection of the constitutive residual read by
       *        the next step.
       */
      void updateSubscalesExp();

      /// @brief The "psi_rb" formulation (Config::formulation): its forms.
      void setupFlowPsiRB();

      /**
       * @brief "psi_rb", as updateSubscalesExp: the subscales of step n+1 and
       *        Pi[R_psi^{n+1}].
       */
      void updateSubscalesPsiRB();

      /**
       * @brief The length of the cell of p along the direction of u (h1 if
       *        |u| vanishes).
       */
      static Real streamlineLength(
        const Rodin::Geometry::Point& p, const Rodin::Math::SpatialVector<Real>& u);
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

      /// @brief z <- a x + b y, with the ghost values refreshed.
      static void combine(Real a, const ::Vec& x, Real b, const ::Vec& y, ::Vec& z);

      /**
       * @brief The time step of the BDF mass terms: dt (BDF1) or, for BDF2
       *        with w = dt_{n+1}/dt_n, dt (1 + w)/(1 + 2 w) (2 dt/3 at w = 1).
       */
      Real massStep() const;

      /**
       * @brief w = dt_{n+1}/dt_n of the BDF2 step about to be taken; 1 on the
       *        first step (x^{-1} = x^0) and for BDF1.
       */
      Real bdfRatio() const;

      /**
       * @brief Before each attempt of a step, from x^n, x^{n-1} and w: the BDF
       *        history x_hist, with d/dt x ~ (x - x_hist)/massStep(),
       *          x_hist = ((1 + w)^2 x^n - w^2 x^{n-1})/(1 + 2 w)  (BDF2),
       *        and the extrapolation u* = (1 + w) u^n - w u^{n-1}.
       */
      void prepareStep();

      /**
       * @brief After a step is accepted (m_uIt, the psi iterates and the
       *        pressure solution hold x^{n+1}): shifts x^n, x^{n-1}.
       */
      void advanceHistory();

      /// @brief After a rejected attempt: the solution fields back to x^n.
      void rejectStep();

      /**
       * @brief The orthogonal projections Pi[grad p], Pi[div sigma] and, if
       *        strain, Pi[div u], Pi[eps(u)] and Pi[(u^n . grad) psi], at step n
       *        (iterate = false) or at the Newton iterate (u^k, p^k, psi^k).
       */
      void projectConsistent(bool iterate, unsigned mask);

      bool solveFlow();
      void computeWallShear();
      void closeCycle(Real elapsed);
      void computeFluxes();
      Real boundaryMeasure(const AttributeSet& tags);

      /// @brief r(t) = (1 - cos(pi t/T_r))/2 for t < T_r, then 1.
      Real ramp(Real t) const;

      /// @brief h_K, the local element size.
      static Real cellSize(const Rodin::Geometry::Point& p);

      /// @brief Codina tau_1 on eta_s, as in the 2D driver.
      Real tau1At(const Rodin::Geometry::Point& p) const;

      // ---- Fluids "gsptt" and "cy" -----------------------------------------
      Real carreauYasuda(Real gammaDot) const;

      /// @brief Steady-shear sPTT factor: f^3 - f^2 - 2 eps lambda^2 gd^2 = 0.
      Real steadyPTTFactor(Real gammaDot) const;

      /// @brief eta_s(gd): eta_s ("sptt"), eta_CY ("cy"), eta_CY - eta_p/f ("gsptt").
      Real solventViscosity(Real gammaDot) const;

      /// @brief gd = sqrt(2 eps:eps) of u^n (BDF1) or u* (BDF2).
      Real shearRateAt(const Rodin::Geometry::Point& p) const;
      Real solventViscosityAt(const Rodin::Geometry::Point& p) const;
      Real referenceViscosity() const;

      /// @brief The solvent viscosity does not depend on the shear rate.
      bool constantSolvent() const
      {
        return m_cfg.fluid == "sptt" || m_cfg.fluid == "owens" ||
          m_cfg.fluid == "newtonian";
      }

      /// @brief (u, p) only: generalised Newtonian fluids (Config::twoField).
      bool twoField() const
      {
        return m_cfg.twoField && (m_cfg.fluid == "cy" || m_cfg.fluid == "newtonian");
      }

      /// @brief The (u, p) system of twoField().
      void setupFlowTwoField();

      /// @brief Eulerian residence time (Config::residenceTime).
      void setupResidence();
      void solveResidence();

      /**
       * @brief Config::jacobianCheck, on the assembled system of Newton
       *        iterate 0 (x0, A0, b0 copies, d the increment); returns after
       *        printing, the caller stops the run.
       */
      void checkJacobian(::Mat A0, ::Vec b0, ::Vec x0, ::Vec x1);
      bool isOwens() const
      {
        return m_cfg.fluid == "owens";
      }

      /**
       * @brief Relaxation time of the steady simple shear at this rate:
       *        lambda, or (eta_CY - eta_plasma)/G for Owens.
       */
      Real steadyLambda(Real gammaDot) const;

      /**
       * @brief Owens: aggregation rate aN0(gd), fragmentation rate b(gd),
       *        steady aggregate size n_st(gd) and mu = lambda(n, gd).
       */
      Real owensA(Real gammaDot) const;
      Real owensBRaw(Real gammaDot) const;
      Real owensB(Real gammaDot) const;
      Real owensNst(Real gammaDot) const;
      Real owensLambda(Real n, Real gammaDot) const;
      Real shearRateOf(
        const VectorGridFunctionType& u, const Rodin::Geometry::Point& p) const;

      /// @brief Owens: aggregate size of the entering blood.
      Real inletAggregateSize() const;

      /// @brief Owens: lambda by cell for the psi equation, or nullptr.
      const CellGridFunctionType* lambdaField() const
      {
        return isOwens() ? &m_lambdaCell : nullptr;
      }
      void setupStructure();

      /**
       * @brief Owens: n^{n+1} from the accepted (u^{n+1}, n^n), then
       *        lambda(n^{n+1}, gd^{n+1}) for the next step.
       */
      void solveStructure();
      void updateStructureFields();

      /// @brief alpha3, with the PTT factor at psi^n.
      Real alpha3At(const Rodin::Geometry::Point& p) const;

      /// @brief The alphas (P0) and every Pi_h, all at step n.
      void updateStabilization();

      /// @brief lambda0 = max(k lambda, lambda0_min).
      Real lambda0() const;

      /// @brief PTT factor f = 1 + epsilon (lambda/lambda0)(tr exp(psi) - 3).
      Real pttFactor(Real traceExp) const;

      /// @brief sigma at the nodes, from (m_psiItD, m_psiItO).
      void updateConformation();

      void writeCSVHeader();
      void writeCSVRow(int cycle);

      Config m_cfg;

      /// @brief Reference speed (m/s) and cardiac period (s).
      Real m_meanVelocity = 0.0;
      Real m_period = 0.0;

      /**
       * @brief Mean pressure of the mitral plane at the last accepted step:
       *        the atrial pressure read by the ventricle (declared before it).
       */
      Real m_outletMeanPressure = 0.0;
      VentricleType::Input m_input;
      VentricleType m_model;

      /**
       * @brief The ventricle is ahead of the flow by at most one step: its
       *        pressure at the two ends of the interval it last covered.
       */
      Real m_lvTime0 = 0.0;
      Real m_lvTime1 = 0.0;
      Real m_lvPressure0 = 0.0;
      Real m_lvPressure1 = 0.0;

      /// @brief Opening data read by the frozen forms.
      Real m_pIn = 0.0;
      Real m_pOut = 0.0;
      Real m_mitralZ = 0.0;
      bool m_valveOpen = true;

      /**
       * @brief dt of the step being taken, its nominal value T/stepsPerCycle,
       *        the dt of the last accepted step and the number of accepted steps.
       */
      Real m_dt = 0.0;
      Real m_dtNominal = 0.0;
      Real m_dtPrev = 0.0;
      int m_history = 0;

      /**
       * @brief (dt, massStep) the flow forms were built with: setupFlow bakes
       *        both into the forms, so a change of either rebuilds them.
       */
      Real m_formDt = -1.0;
      Real m_formMassStep = -1.0;

      MeshType m_parent;
      FluidMeshType m_mesh;
      Rodin::IO::XDMF m_xdmf;
      AttributeSet m_wallSet;
      AttributeSet m_inletSet;

      VectorFESType m_vh;
      ScalarFESType m_sh;

      /// @brief Space of one block of psi: P1, three components.
      VectorFESType m_tauh;
      CellFESType m_ch;

      // ---- Flow ------------------------------------------------------------
      VectorTrialFunctionType m_u;
      ScalarTrialFunctionType m_p;
      VectorTestFunctionType m_v;
      ScalarTestFunctionType m_q;

      /// @brief u^n: the stabilisation parameters are evaluated at it.
      VectorGridFunctionType m_uOld;

      /// @brief Newton velocity u^k of the constitutive equation.
      VectorGridFunctionType m_uIt;

      /**
       * @brief u^{n-1} (BDF2); the BDF history, u^n or (4 u^n - u^{n-1})/3, so
       *        that d/dt u ~ (u - m_uHist)/dt_eff, dt_eff = dt or 2 dt/3; and the
       *        advecting velocity of the Galerkin convection, u^n or u*.
       */
      VectorGridFunctionType m_uNm1;
      VectorGridFunctionType m_uHist;
      VectorGridFunctionType m_uAdv;

      /// @brief p^n, p^{n-1} (BDF2) and the pressure iterate p^k.
      ScalarGridFunctionType m_pN;
      ScalarGridFunctionType m_pNm1;
      ScalarGridFunctionType m_pIt;

      /**
       * @brief psi = log(tau): diagonal (xx, yy, zz) and off-diagonal
       *        (xy, xz, yz) blocks; the test tensor chi likewise.
       */
      VectorTrialFunctionType m_psiD;
      VectorTrialFunctionType m_psiO;
      VectorTestFunctionType m_chiD;
      VectorTestFunctionType m_chiO;
      VectorGridFunctionType m_psiOldD;
      VectorGridFunctionType m_psiOldO;

      /// @brief psi^{n-1} and the BDF history of psi, as for u.
      VectorGridFunctionType m_psiNm1D;
      VectorGridFunctionType m_psiNm1O;
      VectorGridFunctionType m_psiHistD;
      VectorGridFunctionType m_psiHistO;

      /// @brief Linearisation point psi^k of the Newton iterations.
      VectorGridFunctionType m_psiItD;
      VectorGridFunctionType m_psiItO;

      /**
       * @brief Bumped whenever (u^k, psi^k) or (u^n, psi^n) change; keys the
       *        cache of the pointwise log-conformation data.
       */
      std::uint64_t m_psiRevision = 0;

      /**
       * @brief Bumped when the step data (psi^n, u^n) change: the revision of
       *        the pointwise data read at step n ("psi_rb" adjoint).
       */
      std::uint64_t m_stepRevision = 0;

      /// @brief sigma = (eta_p/lambda0)(exp(psi) - I) at the nodes, by block.
      VectorGridFunctionType m_sigmaD;
      VectorGridFunctionType m_sigmaO;

      // ---- Reusable mass projection (wall-shear recovery) -------------------
      VectorTrialFunctionType m_wTrial;
      VectorTestFunctionType m_wTest;

      // ---- VMS, split OSS (updateStabilization) ----------------------------
      CellGridFunctionType m_tauK; ///< dynamic convective tau_K
      CellGridFunctionType m_alpha1; ///< tau_1/rho
      CellGridFunctionType m_alpha2; ///< rho h^2/(4 tau_1)
      CellGridFunctionType m_alpha3;
      CellGridFunctionType m_alphaPsi; ///< 2 k^2 s alpha3, the psi advection weight
      Heart::OrthogonalProjection<VectorFESType> m_piConv; ///< Pi[(grad u^n) u^n]
      Heart::OrthogonalProjection<VectorFESType> m_sub; ///< dynamic subscale u'
      VectorGridFunctionType m_subOld;
      Heart::OrthogonalProjection<ScalarFESType> m_piDiv; ///< Pi[div u^n]
      Heart::OrthogonalProjection<VectorFESType> m_piGradP; ///< Pi[grad p^n]
      Heart::OrthogonalProjection<VectorFESType> m_piDivSigma; ///< Pi[div sigma^n]
      Heart::OrthogonalProjection<VectorFESType> m_piEpsD; ///< Pi[eps(u^n)], D block
      Heart::OrthogonalProjection<VectorFESType> m_piEpsO; ///< Pi[eps(u^n)], O block
      Heart::OrthogonalProjection<VectorFESType> m_piAdvPsiD; ///< Pi[(u^n.grad) psi^n], D
      Heart::OrthogonalProjection<VectorFESType> m_piAdvPsiO; ///< Pi[(u^n.grad) psi^n], O

      /**
       * @brief "exp": Pi[R^n] of the constitutive residual (D and O blocks),
       *        the static and dynamic alpha3 and the PTT factor f(psi^n), by cell.
       */
      Heart::OrthogonalProjection<VectorFESType> m_piRD;
      Heart::OrthogonalProjection<VectorFESType> m_piRO;
      CellGridFunctionType m_alpha3Dyn;
      CellGridFunctionType m_pttN;

      /**
       * @brief Owens: the aggregate size n (P1; n^n, n^{n-1} and its BDF
       *        history) and lambda(n^n, gd(u^n)) by cell, read by the psi
       *        equation and its stabilisation during step n + 1.
       */
      ScalarTrialFunctionType m_n;
      ScalarTestFunctionType m_nTest;
      ScalarGridFunctionType m_nOld;
      ScalarGridFunctionType m_nNm1;
      ScalarGridFunctionType m_nHist;
      CellGridFunctionType m_lambdaCell;
      ScalarProblemType m_structure;
      Rodin::Solver::KSP m_structureKSP;
      Real m_nMin = 1.0;
      Real m_nMax = 1.0;
      Real m_lambdaMax = 0.0;
      int m_dumps = 0;

      /// @brief Owens: the maximum of b(gd) and where it is (bFloor).
      Real m_owensBPeak = 0.0;
      Real m_owensGdPeak = std::numeric_limits<Real>::infinity();

      // ---- Wall shear indices ----------------------------------------------
      /// @brief The three rows of 2 eps(u), recovered onto the nodes.
      VectorGridFunctionType m_wss;
      VectorGridFunctionType m_symRec0;
      VectorGridFunctionType m_symRec1;
      VectorGridFunctionType m_symRec2;
      VectorGridFunctionType m_netShear;
      ScalarGridFunctionType m_absShear;
      ScalarGridFunctionType m_shearMagnitude;
      ScalarGridFunctionType m_tawss;
      ScalarGridFunctionType m_osi;

      /// @brief TAWSS and OSI on the wall of the appendage only.
      ScalarGridFunctionType m_tawssAppendage;
      ScalarGridFunctionType m_osiAppendage;

      /**
       * @brief Nodal wall and appendage indicators, int_G phi_i dA (assembled,
       *        ADD across ranks): consistent on every rank, which a
       *        Region::Boundary projection is not (it visits the OWNED faces
       *        only, and flush() lets a stale ghost overwrite the owner:
       *        seams of zeros along the partition, as in AtriumPTAG).
       */
      ScalarGridFunctionType m_wallArea;
      ScalarGridFunctionType m_appendageArea;

      /**
       * @brief Config::wallKinematics: wall-node indicator (0/1) and the
       *        area-weighted nodal wall normal, by component.
       */
      ScalarGridFunctionType m_wallFlag;
      ScalarGridFunctionType m_wallN0;
      ScalarGridFunctionType m_wallN1;
      ScalarGridFunctionType m_wallN2;
      VectorGridFunctionType m_wssRec;
      Real m_wallAreaThreshold = 0.0;
      Real m_appendageAreaThreshold = 0.0;

      // ---- Fluxes ----------------------------------------------------------
      ScalarTestFunctionType m_qFlux;
      ScalarGridFunctionType m_one;
      FluxFormType m_flux;

      /// @brief "exp": the sweep that writes the subscales at the quadrature points.
      FluxFormType m_subscaleSweep;

      // ---- Problems and solvers --------------------------------------------
      FlowProblemType m_flow;
      Rodin::Solver::KSP m_flowKSP;

      /// @brief twoField(): (u, p; v, q), the same trial and test functions.
      using TwoFieldProblemType =
        Rodin::Variational::Problem<LinearSystemType, VectorTrialFunctionType,
          ScalarTrialFunctionType, VectorTestFunctionType, ScalarTestFunctionType>;
      TwoFieldProblemType m_flow2;
      Rodin::Solver::KSP m_flow2KSP;
      bool m_flow2FieldSplitsSet = false;
      bool m_jacobianChecked = false;

      /**
       * @brief Eulerian residence time a (P1), its last two steps, BDF
       *        history, and the volume and appendage means.
       */
      ScalarTrialFunctionType m_age;
      ScalarTestFunctionType m_ageTest;
      ScalarGridFunctionType m_ageOld;
      ScalarGridFunctionType m_ageNm1;
      ScalarGridFunctionType m_ageHist;
      ScalarProblemType m_ageProblem;
      Rodin::Solver::KSP m_ageKSP;
      Real m_ageMean = 0.0;
      Real m_ageAppendage = 0.0;
      Real m_appendageVolume = 0.0;
      VectorProblemType m_vectorProjection;
      Rodin::Solver::KSP m_vectorProjectionKSP;
      VectorTrialFunctionType m_wssTrial;
      VectorTestFunctionType m_wssTest;
      VectorProblemType m_wssProjection;
      Rodin::Solver::KSP m_wssKSP;

      // ---- Time-dependent coefficients read by the frozen forms ------------
      Real m_t = 0.0;
      Real m_inletMeasure = 1.0;
      Real m_outletMeasure = 1.0;
      Real m_qIn = 0.0;
      Real m_qOut = 0.0;
      Real m_inletPressure = 0.0;
      Real m_volume = 1.0;

      /// @brief Config::metrics (README.md), SI units.
      struct Metrics
      {
        Real kinetic = 0.0; ///< (rho/2) int |u|^2
        Real enstrophy = 0.0; ///< (1/2) int |curl u|^2
        Real viscousDissipation = 0.0; ///< int eta_s gd^2
        Real polymerPower = 0.0; ///< int sigma : eps(u)
        Real elasticEnergy = 0.0; ///< int tr(sigma)/2
        Real polymerDissipation = 0.0; ///< int f tr(sigma)/(2 lambda)
        Real convectiveSubscale = 0.0; ///< rho^2 int tau |P^perp[(grad u) u]|^2
        Real gradDivSubscale = 0.0; ///< int alpha2 (div u - Pi[div u]) div u
        Real velocitySubscale = 0.0; ///< (rho/2) int |u1 + u2 + u3|^2
        Real psiNorm = 0.0; ///< sqrt(int |psi|^2 / V)
        Real psiSubscale = 0.0; ///< sqrt(int |psi~|^2 / V)
        Real wiMax = 0.0; ///< max lambda gd
        Real wiVolume = 0.0; ///< volume fraction with lambda gd > 1
        Real wiAppendage = 0.0; ///< mean lambda gd in the appendage
        Real shearAppendage = 0.0; ///< mean gd in the appendage
        Real kineticAppendage = 0.0; ///< (rho/2) int |u|^2 over the appendage
      };
      Metrics m_metrics;
      CellGridFunctionType m_wiCell;

      /// @brief Probes: centres, normalisations int w, and the output file.
      std::vector<std::array<Real, 3>> m_probes;
      std::vector<Real> m_probeWeights;
      std::ofstream m_probeCSV;

      /// @brief Owens::resetCycle: int a dt and int b dt of the current cycle.
      ScalarGridFunctionType m_aSum;
      ScalarGridFunctionType m_bSum;
      ScalarGridFunctionType m_rateWork;

      /// @brief Cycles completed.
      int m_cyclesDone = 0;
      Real m_speed = 0.0;

      /// @brief max |sigma_ij| (Pa), reported with the flow.
      Real m_stress = 0.0;

      /// @brief Newton iterations of the last step and their last increment.
      int m_conformationIts = 0;
      Real m_psiIncrement = 0.0;

      /**
       * @brief Undamped Newton iterations of the last step: the measure of
       *        its difficulty for the step control (the damped ones of the
       *        start from rest are not).
       */
      int m_newtonUndamped = 0;

      /// @brief Why the last attempt was rejected (empty if accepted).
      std::string m_failure;

      /// @brief Newton accepted by the quadratic prediction (not by d_k < tol).
      bool m_newtonPredicted = false;

      /// @brief Recent estimates of C in d_{k+1} ~ C d_k^2 (ring buffer).
      std::vector<Real> m_contraction;
      size_t m_contractionNext = 0;

      /**
       * @brief Lagged preconditioner: built at least once, solves since its
       *        last build, and a forced rebuild (dt changed).
       */
      bool m_pcBuilt = false;
      int m_pcSolves = 0;

      /**
       * @brief Krylov iterations of the first Newton solve after the last
       *        build, the reference of pcReuseGrowth (< 0 until measured).
       */
      int m_pcReferenceIts = -1;
      bool m_pcStale = true;
      bool m_flowFieldSplitsSet = false;
      bool m_initialized = false;

      int m_step = 0;
      Timing m_timing;

      std::ofstream m_csv;
  };
}

#endif
