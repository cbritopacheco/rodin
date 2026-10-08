// ForcedConvection.cpp
//
// Flow boiling of water in a 3D channel with the Lee phase-change model.
// One-fluid formulation with vapour fraction phi in [0, 1]:
//
//   d phi/dt + div(phi u) = gamma div(eps grad phi - phi (1 - phi) n) + mdot / rho_v
//   div u                                         = mdot (1/rho_v - 1/rho_l)
//   rho (u - u^n)/dt + rho (grad u) u - div(2 mu e(u)) + grad p
//                                                 = rho g + sigma kappa grad phi
//   rho cp (dT/dt + u . grad T) - div(k grad T)   = - mdot h_lv
//
// with mixture properties phi-weighted (arithmetic for rho, mu, rho cp;
// harmonic for k) and the Lee source
//
//   mdot = r_l (1 - phi) rho_l [theta]_+  +  r_v phi rho_v [theta]_-,
//   theta = (T - T_sat)/T_sat,
//
// evaporation only above T_sat, condensation only where there is vapour, and
// mdot continuous at T_sat. The ramps are smoothed with a softplus of width
// dTsmooth ([theta]_+ >= 0, [theta]_- = theta - [theta]_+ <= 0), NOT by
// blending the two rate coefficients: a blended switch evaporates "negative
// vapour" from pure subcooled liquid and grows like exp(2 T/dTsmooth).
// The source is linearised (Newton) in (phi, T) about the previous Picard
// iterate: dmdot/dphi <= 0 and dmdot/dT >= 0, so both reaction terms are
// positive on the left-hand side.
//
// The conservative Allen-Cahn term (Chiu & Lin 2011) keeps the interface at a
// fixed thickness eps (tanh profile) and phi in [0, 1]; n = Pi_h[grad phi /
// |grad phi|] is shared with the curvature. The compression is semi-implicit:
// gamma ((1 - phi^k) n^k phi, grad a). gamma is a velocity scale (default: the
// peak inlet velocity) and eps defaults to 1.5 mean cell sizes.
//
// Discretisation: backward Euler. (phi, T) are solved monolithically in P1/P1 (P2 with -DFC_PHASE_ORDER=2),
// (u, p) monolithically in equal-order P1/P1. Both blocks are stabilised with
// split orthogonal subgrid scales (OSS): tau (A(u_h) - Pi_h[A(u_h^n)], B(v_h))
// with Pi_h the L2 projection from OrthogonalProjection.h, lagged at t^n. The
// two blocks are coupled by Picard iterations within each time step
// (properties and mdot refreshed every iteration; normal and curvature at t^n
// unless -fc_lag_interface false).
// Curvature is kappa = -div Pi_h[grad phi / |grad phi|].
//
// Run from the build directory:
//   mpirun -n 4 ./examples/Heat/ForcedConvection -fc_dt 1e-5 -fc_tend 0.2
//   ./examples/Heat/ForcedConvection -fc_2d -fc_tend 0.2   (2D x-z section)
//   ./examples/Heat/ForcedConvection -fc_topview -fc_fluid r134a   (2D top view)
//
// Top view (-fc_topview): the diverging channel seen from above, plane x-y,
// width W(x) linear from inletWidth to outletWidth, depth d, every field
// averaged over the depth (Hele-Shaw / Brinkman):
//   momentum  ... + (C_f mu/d^2) u            friction of floor and ceiling
//   fluid T   ... + (H_w/d)(T - T_w)          exchange with the heated walls
//   wall T_w  C_w dT_w/dt - div(K_w grad T_w) + H_w (T_w - T) + H_f (T_w - T_sat) = q''
//   film      mdot_f = H_f (T_w - T_sat)/(h_lv d), in continuity and the phi equation
// with H = h P_h/W per unit plan area (P_h the heated perimeter: floor and
// side walls by default, -fc_heated_walls), h_w = Nu k/D_h(x) (k
// phi-weighted) and h_f = F(phi) k_l/delta_f the
// evaporating liquid film left on the floor under vapour (thin-film
// evaporation; delta_f ~ Taylor-Bretherton, a few microns). The heat flux q''
// of the floor (W/m^2) is no longer a boundary condition but the source of
// the wall equation; its latent part h_f (T_w - T_sat) leaves the wall as
// vapour without heating the depth-averaged fluid. Nucleation is triggered by
// T_w. Options: -fc_wi, -fc_wo, -fc_depth, -fc_darcy, -fc_nu, -fc_heated_walls, -fc_film,
// -fc_fluid (water-paper, water, ethanol, r245fa, r134a, r1234ze: saturated at
// T_sat = 310 K, p_out = p_sat, T_in = T_sat - 10 K unless -fc_tin is given).
//
// Coefficients (mixture properties, Lee terms in P2; tau in P0) are finite
// element functions refreshed every Picard iteration, so every integrand has a
// known polynomial order and the implicit and lagged halves of each term are
// integrated with the same rule (see updateCoefficients()).
//
// Options: -fc_mesh, -fc_dt, -fc_tend, -fc_G, -fc_q, -fc_tin, -fc_r, -fc_g,
// -fc_eps, -fc_gamma (0 disables Allen-Cahn), -fc_maxit, -fc_tol, -fc_output_every,
// -fc_print_every (progress line on screen, default every step), -fc_picard_monitor,
// -fc_neps (normal regularisation, 1/m; default 0.05/eps), -fc_dtsmooth (K, Lee switch width),
// -fc_lag_interface (n, kappa at t^n; default true), -fc_2d, -fc_nx, -fc_ny, -fc_length,
// -fc_clapeyron (T_sat(p), default true), -fc_dc (discontinuity capturing on T, 0 = off), -fc_wall_capacity (J/(m^2 K), substrate under the heated bottom), -fc_wall_conductance (W/K,
// k_s t_s of that substrate, lateral conduction), -fc_pout (Pa), -fc_tsat (K),
// -fc_interface_transfer (default true; false = bulk Lee), -fc_accommodation (s_e), -fc_area_cutoff,
// -fc_site_spacing (m, 0 = no
// nucleation sites), -fc_onset (K), -fc_seed_radius (m), -fc_site_wait (s), -fc_nucleation_every,
// -fc_dt_boil (s, dt cap after the first nucleation), -fc_adaptive_dt (Picard-driven dt after
// nucleation, default true), -fc_picard_easy, -fc_picard_hard, -fc_dt_grow, -fc_dt_shrink,
// -fc_dt_min, -fc_max_retries (rejected steps are retried with dt/2)
// (2D channel); solvers under the prefixes -fc_flow_
// (MUMPS), -fc_phase_ (GMRES) and -fc_proj_* (CG + Jacobi mass solves).
#include "Rodin/Variational/ForwardDecls.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

#include <petscksp.h>

#include <Rodin/Alert.h>
#include <Rodin/Configure.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/MPI.h>
#include <Rodin/PETSc.h>
#include <Rodin/Solver.h>
#include <Rodin/Types.h>
#include <Rodin/Variational.h>

#ifdef RODIN_USE_SCOTCH
#include <Rodin/Scotch/MeshPartitioner.h>
#endif

#include "OrthogonalProjection.h"

/// Polynomial order of (phi, T, T_w). P2 Galerkin has no discrete maximum
/// principle at all (its mass and stiffness matrices are not M-matrices):
/// across an interface resolved by 1-2 cells the edge-midpoint DOFs of T
/// overshoot by tens of K in the vapour (rho cp 5e4 times below the liquid)
/// while the vertex values stay bounded (seen: 326 K at midpoints, 313 K at
/// vertices, wall at 316 K; the case of the 540 K blow-up). P1 with the
/// discontinuity capturing on T is close to monotone. Compile with
/// -DFC_PHASE_ORDER=2 for the previous P2 discretisation.
#ifndef FC_PHASE_ORDER
#define FC_PHASE_ORDER 1
#endif

namespace Rodin::Examples::Heat
{
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;

  class ForcedConvection
  {
    public:
      using MeshType = Mesh<Context::MPI>;
      using VectorFES = H1<1, Math::SpatialVector<Real>, MeshType>;
      using ScalarFES = H1<1, Real, MeshType>;
      /// Space of (phi, T, T_w) and of the coefficients: P1 by default (see
      /// FC_PHASE_ORDER at the top).
      using ScalarFES2 = H1<FC_PHASE_ORDER, Real, MeshType>;
      using VectorGF = PETSc::Variational::GridFunction<VectorFES>;
      using ScalarGF = PETSc::Variational::GridFunction<ScalarFES>;
      using ScalarGF2 = PETSc::Variational::GridFunction<ScalarFES2>;
      using VectorTrial = PETSc::Variational::TrialFunction<VectorGF, VectorFES>;
      using ScalarTrial = PETSc::Variational::TrialFunction<ScalarGF, ScalarFES>;
      using ScalarTrial2 = PETSc::Variational::TrialFunction<ScalarGF2, ScalarFES2>;
      using VectorTest = PETSc::Variational::TestFunction<VectorFES>;
      using ScalarTest = PETSc::Variational::TestFunction<ScalarFES>;
      using ScalarTest2 = PETSc::Variational::TestFunction<ScalarFES2>;
      using System = PETSc::Math::LinearSystem;
      using FlowProblem = Problem<System, VectorTrial, ScalarTrial, VectorTest, ScalarTest>;
      using PhaseProblem = Problem<System, ScalarTrial2, ScalarTrial2, ScalarTest2, ScalarTest2>;
      /// Top view: (phi, T, T_w) monolithic.
      using PhaseProblemTop = Problem<System, ScalarTrial2, ScalarTrial2, ScalarTrial2,
                                      ScalarTest2, ScalarTest2, ScalarTest2>;
      using CellFES = P0<Real, MeshType>;
      using CellGF = PETSc::Variational::GridFunction<CellFES>;

      struct Labels
      {
        Attribute inlet = 1, outlet = 2, bottom = 3, top = 4, sides = 5;
      };

      struct Config
      {
        std::string meshPath = "../resources/examples/Heat/DivergingMicrochannel3D.mesh";
        std::string output = "ForcedConvection";
        Labels labels;

        // Water as in Alugoju, Dubey & Javed, IJHMT 160 (2020) 120212, Sec. 2.4:
        // outlet at P_sat = 6500 Pa (absolute), T_sat = 310 K. The paper does
        // not give h_lv and prints sigma = 0.662 N/m (a typo): both are taken
        // from saturated water at 6500 Pa (IAPWS, CoolProp).
        Real rhoL = 996.5;     ///< kg/m^3 (paper)
        Real rhoV = 0.045;     ///< kg/m^3 (paper; rho_l/rho_v ~ 22000)
        Real muL = 2.82e-4;    ///< Pa s (paper)
        Real muV = 3.0e-5;     ///< Pa s (paper)
        Real kL = 0.606;       ///< W/(m K) (paper)
        Real kV = 0.0202;      ///< W/(m K) (paper)
        Real cpL = 4179.0;     ///< J/(kg K) (paper)
        Real cpV = 1882.0;     ///< J/(kg K) (paper)
        Real hLV = 2.412e6;    ///< J/kg (saturation at 6500 Pa)
        Real sigma = 0.0701;   ///< N/m (saturation at 6500 Pa)
        Real Tsat = 310.0;     ///< K (paper), at p = outletPressure
        /// Saturation temperature from the local pressure (Clausius-Clapeyron,
        /// through (outletPressure, Tsat)): 1/T_sat(p) = 1/T_sat - R/(M h_lv) ln(p/p_out).
        /// At 6500 Pa, dT_sat/dp = 2.8e-3 K/Pa: 350 Pa of overpressure in a
        /// growing bubble cancel 1 K of wall superheat. With a fixed T_sat the
        /// evaporation ignores the pressure it creates and a confined bubble
        /// in front of the fixed-velocity inlet drives MPa spikes. The
        /// dependence is implicit in the continuity equation,
        ///   div u + beta (p - p^k) = (1/rho_v - 1/rho_l) mdot(p^k),
        ///   beta = (1/rho_v - 1/rho_l) (dmdot/dT) dT_sat/dp >= 0,
        /// (a compressibility-like term), and lagged (p^k) in the phase block.
        bool clapeyron = true;
        /// Discontinuity capturing on T (Codina 1993, orthogonal form):
        ///   k_dc = C_dc (h/2) |u.grad T - Pi[u.grad T]| / |grad T|  (<= h|u|/2),
        /// added as rho cp k_dc (grad T, grad w), lagged at the Picard iterate.
        /// Galerkin + OSS is not monotone; the vapour has rho cp ~5e4 times
        /// below the liquid, so a small overshoot there is tens of K and
        /// grows with the vapour velocity (seen: T 310 -> 540 K at 30 m/s
        /// with the wall at 316 K). Vanishes where the solution is resolved.
        /// 0 disables.
        Real dcT = 0.7;

        /// Thin silicon substrate under the heated bottom (conjugate effect of
        /// the paper, lumped across its 0.4 mm thickness): heat capacity per
        /// unit area rho_s c_s t_s = 2330 * 712 * 4e-4 J/(m^2 K), added as
        /// C_w dT/dt on the bottom boundary. Lateral conduction in the silicon
        /// is not represented. 0 disables (flux straight into the fluid).
        Real wallCapacity = 2330.0 * 712.0 * 4.0e-4;
        /// Lateral conduction in the same substrate, k_s t_s = 148 * 4e-4 W/K:
        /// thin-wall model on the heated bottom,
        ///   C_w dT/dt - div_G(k_s t_s grad_G T) + q_fluid = q'',
        /// with grad_G T = grad T - (grad T . n) n (the bottom is flat). It
        /// carries the heat entering under a vapour patch to the wetted wall
        /// around it: ~6e4 W/(m^2 K) over 1 mm. Without it a dry patch takes
        /// the full 150 kW/m^2 into vapour (rho c_p ~5e4 times smaller than
        /// the liquid), T runs away and the Hertz-Knudsen source with it.
        /// 0 disables.
        Real wallConductance = 148.0 * 4.0e-4;

        // Lee model. r is not physical: calibrate on the 1D Stefan problem.
        Real rL = 100.0;       ///< 1/s, evaporation
        Real rV = 100.0;       ///< 1/s, condensation
        Real dTsmooth = 0.05;  ///< K, width of the softplus ramps at T_sat

        /// Interfacial mass transfer (default) instead of the bulk Lee model:
        /// the Hertz-Knudsen flux of the paper, Eq. (12),
        ///   F = 2 s_e/(2 - s_e) sqrt(M/(2 pi R T_sat)) h_lv rho_l rho_v/(rho_l - rho_v) (T - T_sat)/T_sat,
        /// times the interfacial area density of the phase field,
        ///   A_i = |grad phi| = phi (1 - phi)/eps   (conservative Allen-Cahn profile),
        /// in place of 6 alpha_v/d (Eq. 13), which needs an arbitrary d. Pure
        /// liquid (phi = 0) does not evaporate: vapour appears at nucleation
        /// sites only. Bulk Lee evaporates any superheated liquid and fills the
        /// superheated layer with a diffuse vapour "fog" instead of bubbles.
        bool interfaceTransfer = true;
        /// s_e. The paper takes 1, but at 6500 Pa (rho_v = 0.045) that gives
        /// F/(T - T_sat) = 0.74 kg/(m^2 s K), an interface speed F/rho_v of
        /// ~16 m/s per K of superheat: a seed then drives MPa pressure
        /// oscillations against the fixed-velocity inlet and Picard stalls.
        /// 0.01 (within the measured 0.01-1 range for water) gives ~0.08 m/s
        /// per K and inlet pressures of ~10 kPa, the range of the paper's
        /// Fig. 10. Calibrate it (e.g. on the 1D Stefan problem).
        Real accommodation = 0.01;
        Real molarMass = 0.018015;   ///< kg/mol
        /// A_i is ramped to zero (C^1) for phi < 2 c: round-off vapour would
        /// otherwise grow at ~F/(eps rho_v) ~ 1e6 1/s and nucleate everywhere.
        Real areaCutoff = 0.01;

        /// Nucleation sites on the heated bottom, every siteSpacing along x
        /// (first at siteSpacing/2; on y = 0 in 3D). A site seeds a vapour
        /// nucleus (tanh profile of radius seedRadius, hemisphere/half disc
        /// on the wall) when the mean wall temperature within seedRadius
        /// exceeds T_sat + onsetSuperheat, the mean temperature over the seed
        /// itself is >= T_sat, no vapour is within 2 seedRadius,
        /// and siteWaiting has elapsed since its last seed. 0 disables.
        Real siteSpacing = 2.0e-3;   ///< m
        Real onsetSuperheat = 1.0;   ///< K
        Real seedRadius = -1.0;      ///< m (auto: 2.5 eps)
        Real siteWaiting = 5.0e-4;   ///< s
        int nucleationEvery = 1;     ///< steps between site checks
        /// Time step from the first nucleation on (<= 0: keep dt). Lets the
        /// single-phase heating stage (~45 ms in the paper) run with a larger dt.
        Real dtBoiling = -1.0;
        /// After the first nucleation dt is controlled by the Picard coupling
        /// between the (phi, T) and (u, p) blocks, not by a Courant number:
        /// each BDF1 block is solved implicitly, but the fixed point that
        /// couples them (phi transported by u^k; density and expansion source
        /// from phi^k) contracts worse as the interface moves a larger part of
        /// its thickness per step. dt grows by dtGrow after a step converged in
        /// at most picardEasy iterations, shrinks by dtShrink after one that
        /// needed more than picardHard, always within [dtMin, dtBoiling (or dt)].
        bool adaptiveDt = true;
        int picardEasy = 3;
        int picardHard = 6;
        Real dtGrow = 1.5;
        Real dtShrink = 0.5;
        /// A step whose Picard iterations do not converge (or whose linear
        /// solve diverges) is rejected and retried with dt/2, at most
        /// maxRetries times and never below dtMin.
        int maxRetries = 6;
        Real dtMin = 1.0e-8;
        Real gravity = 0.0;    ///< m/s^2 along -z; 0 keeps the microchannel horizontal
        /// 1/m, regularisation of n = grad phi / sqrt(|grad phi|^2 + neps^2).
        /// Must be comparable to interface gradients (~1/(4 eps)), not ~0: with
        /// neps -> 0, n is a unit vector even where phi ~ 1e-5, and the
        /// compression gamma div(phi n) turns into an anti-diffusive reaction of
        /// rate ~gamma |div n| ~ gamma/h. Negative = automatic: 0.05/eps.
        Real normalEps = -1.0;

        // Conservative Allen-Cahn interface regularisation. Negative = automatic.
        // Defaults are mild: eps = h_mean (interface ~4 eps wide, the channel
        // is 0.4 mm deep) and gamma = mean inlet velocity. With interfacial
        // transfer there is no vapour fog, so both terms act on real
        // interfaces only (they vanish identically where phi = 0 or 1).
        Real epsilon = -1.0;   ///< m, interface half-thickness (auto: h_mean)
        Real gamma = -1.0;     ///< m/s, mobility (auto: mean inlet velocity G/rho_l); 0 disables

        Real massFlux = 240.0;          ///< G (kg/m^2 s); mean inlet velocity G/rho_l
        Real inletWidth = 2.0e-4;       ///< m, inlet centred on y = 0
        Real depth = 4.0e-4;            ///< m, z in [0, depth]
        Real inletTemperature = 300.0;  ///< K
        Real heatFlux = 1.5e5;          ///< W/m^2 on the bottom
        Real outletPressure = 6500.0;   ///< Pa (absolute, paper: P_sat at the outlet)

        /// 2D variant: the x-z section of the 3D channel, [0, length] x [0, depth],
        /// built in code (inlet x = 0, outlet x = length, heated bottom y = 0,
        /// adiabatic no-slip top y = depth); the inflow is a single parabola.
        bool twoD = false;
        Real length = 3.0e-2;  ///< m (as the 3D channel)
        size_t nx = 262;       ///< streamwise cells (dx = 0.115 mm, as the .geo)
        size_t ny = 10;        ///< cells across the depth (nD in the .geo)

        /// Top view (see the header): plane x-y, y in [-W(x)/2, W(x)/2],
        /// nx x ny cells (defaults 600 x 12 with -fc_topview).
        bool topView = false;
        Real outletWidth = 4.0e-4;     ///< m (inletWidth at x = 0)
        /// C_f of the depth friction C_f mu/d^2 u. Hele-Shaw (W >> d) gives 12;
        /// with the in-plane viscous term (Brinkman) 14 reproduces the square
        /// duct (fRe = 56.9 on D_h) and is 8% low at W = d/2 (fRe = 62.2).
        Real darcy = 14.0;
        /// Heated walls in contact with the fluid: 3 = floor and both side
        /// walls (channel etched in silicon, glass cover as in the paper), 1 =
        /// floor only. The exchange and film areas per unit plan area are
        /// P_h/W, P_h = W + 2d (3) or W (1).
        int heatedWalls = 3;
        Real nusselt = 3.5;            ///< wall-fluid Nu on D_h (3 heated walls, W/d = 0.5-1)
        Real filmThickness = 3.0e-6;   ///< m, liquid film under vapour
        std::string fluid = "water-paper";

        Real dt = 1.0e-5;
        Real tEnd = 0.02;
        Real ramp = 5.0e-3;  ///< s, smooth start of the inflow
        int maxIterations = 10;   ///< Picard iterations per step
        Real tolerance = 1.0e-4;  ///< relative increment of (u, phi, T)
        int outputEvery = 20;
        int outputEveryBoiling = -1;  ///< output interval after the first nucleation (<= 0: outputEvery)
        int printEvery = 1;       ///< steps between progress lines on screen
        bool picardMonitor = false; ///< print (dPhi, dT, dU) at every Picard iteration
        /// Normal n and curvature kappa at t^n (true) or at the Picard iterate
        /// (false). At the iterate, the compression gamma (1 - phi^k) n^k phi
        /// is a lagged term quadratic in phi when |grad phi| << neps, and
        /// Picard oscillates (even/odd) with a gain that grows with phi.
        bool lagInterface = true;
      };

      ForcedConvection(const Context::MPI& context, const Config& cfg)
        : m_cfg(cfg),
          m_mesh(makeMesh(context, m_cfg)),
          m_xdmf(context.getCommunicator(), m_cfg.output),
          m_vh(std::integral_constant<size_t, 1>{}, m_mesh, m_mesh.getSpaceDimension()),
          m_sh(std::integral_constant<size_t, 1>{}, m_mesh),
          m_sh2(std::integral_constant<size_t, FC_PHASE_ORDER>{}, m_mesh),
          m_ch(m_mesh),
          m_u(m_vh), m_p(m_sh), m_v(m_vh), m_q(m_sh),
          m_phi(m_sh2), m_T(m_sh2), m_a(m_sh2), m_w(m_sh2), m_Tw(m_sh2), m_b(m_sh2),
          m_uOld(m_vh), m_uIt(m_vh), m_pIt(m_sh),
          m_phiOld(m_sh2), m_TOld(m_sh2), m_phiIt(m_sh2), m_TIt(m_sh2), m_TwOld(m_sh2), m_TwIt(m_sh2),
          m_phiOut(m_sh), m_TOut(m_sh), m_TwOut(m_sh), m_TsatOut(m_sh),
          m_g(m_vh),
          m_rhoGF(m_sh2), m_muGF(m_sh2), m_kGF(m_sh2), m_rhoCpGF(m_sh2), m_phiClipGF(m_sh2),
          m_mdotGF(m_sh2), m_cTGF(m_sh2), m_cPhiGF(m_sh2),
          m_srcGF(m_sh2), m_srcTGF(m_sh2), m_srcPhiGF(m_sh2),
          m_mdotTotGF(m_sh2), m_srcTwGF(m_sh2), m_hwGF(m_sh2), m_hfGF(m_sh2),
          m_betaGF(m_sh2), m_TsatGF(m_sh2),
          m_tauMGF(m_ch), m_tauMRhoGF(m_ch), m_tauTGF(m_ch), m_tauPhiGF(m_ch), m_tauCGF(m_ch), m_dcTGF(m_ch),
          m_piConv(m_vh, "fc_proj_conv_"), m_piGradP(m_vh, "fc_proj_gradp_"),
          m_piNormal(m_vh, "fc_proj_normal_"), m_piDiv(m_sh, "fc_proj_div_"),
          m_piKappa(m_sh, "fc_proj_kappa_"),
          m_piConvPhi(m_sh2, "fc_proj_phi_"), m_piConvT(m_sh2, "fc_proj_T_"),
          m_flow(m_u, m_p, m_v, m_q), m_flowKSP(m_flow),
          m_phase(m_phi, m_T, m_a, m_w), m_phaseKSP(m_phase),
          m_phaseTop(m_phi, m_T, m_Tw, m_a, m_w, m_b), m_phaseTopKSP(m_phaseTop),
          m_z(m_sh), m_one(m_sh), m_flux(m_z)
      {
        const size_t dim = m_mesh.getSpaceDimension();
        Math::SpatialVector<Real> zero(dim);
        for (size_t i = 0; i < dim; ++i)
          zero(i) = 0.0;
        m_uOld = zero;
        m_uIt = zero;
        m_pIt = m_cfg.outletPressure;  // T_sat(p^k) is evaluated from the first iteration
        m_TOld = m_cfg.inletTemperature;
        m_TIt = m_cfg.inletTemperature;
        m_phiOld = 0.0;  // all liquid
        m_phiIt = 0.0;
        m_TwOld = m_cfg.inletTemperature;
        m_TwIt = m_cfg.inletTemperature;
        m_one = Real(1);
        {
          Math::SpatialVector<Real> g(dim);
          for (size_t i = 0; i < dim; ++i)
            g(i) = 0.0;
          g(dim - 1) = m_cfg.topView ? 0.0 : -m_cfg.gravity;  // top view: g normal to the plane
          m_g = g;
        }

        m_flowKSP.setPrefix("fc_flow_");
        m_phaseKSP.setPrefix("fc_phase_");
        m_phaseTopKSP.setPrefix("fc_phase_");

        const auto& L = m_cfg.labels;
        m_bottomArea = m_cfg.topView ? volumeIntegral(m_one) : boundaryIntegral(m_one, L.bottom);
        m_inletArea = boundaryIntegral(m_one, L.inlet);
        m_outletArea = boundaryIntegral(m_one, L.outlet);

        // Allen-Cahn parameters from the mesh and the inflow when not given.
        const Real hMean = std::pow(volumeIntegral(m_one) / static_cast<Real>(m_mesh.getCellCount()),
                                    1.0 / static_cast<Real>(m_mesh.getDimension()));
        m_tauDt = (m_cfg.dtBoiling > 0.0) ? std::min(m_cfg.dtBoiling, m_cfg.dt) : m_cfg.dt;
        m_eps = (m_cfg.epsilon > 0.0) ? m_cfg.epsilon : hMean;
        m_normalEps = (m_cfg.normalEps > 0.0) ? m_cfg.normalEps : 0.05 / m_eps;
        m_gamma = (m_cfg.gamma >= 0.0) ? m_cfg.gamma : m_cfg.massFlux / m_cfg.rhoL;

        // Hertz-Knudsen coefficient (kg/(m^2 s) per unit (T - T_sat)/T_sat).
        {
          const Real se = m_cfg.accommodation, Rg = 8.314462618;
          m_hk = 2.0 * se / (2.0 - se) * std::sqrt(m_cfg.molarMass / (2.0 * M_PI * Rg * m_cfg.Tsat))
               * m_cfg.hLV * m_cfg.rhoL * m_cfg.rhoV / (m_cfg.rhoL - m_cfg.rhoV);
        }
        m_seedRadius = (m_cfg.seedRadius > 0.0) ? m_cfg.seedRadius : 2.5 * m_eps;
        if (m_cfg.siteSpacing > 0.0)
        {
          const Real length = (m_cfg.twoD || m_cfg.topView) ? m_cfg.length : 3.0e-2;
          for (Real x = 0.5 * m_cfg.siteSpacing; x < length; x += m_cfg.siteSpacing)
            m_sites.push_back(x);
          m_lastSeed.assign(m_sites.size(), -1.0e30);
        }

        // Projections must hold consistent data before the forms are assembled.
        updateStabilisation();
        updateInterface();
        updateCoefficients();

        setupFlow();
        setupPhase();

        m_xdmf.setMesh(m_mesh);
        m_xdmf.add("velocity", m_u.getSolution());
        m_xdmf.add("pressure", m_p.getSolution());
        m_xdmf.add("temperature", m_TOut);  // phase fields interpolated to P1 for output
        m_xdmf.add("fraction", m_phiOut);
        if (m_cfg.topView)
          m_xdmf.add("wall_temperature", m_TwOut);
        m_xdmf.add("saturation_temperature", m_TsatOut);  // T_sat(p)
        m_xdmf.add("curvature", m_piKappa.get());

        const auto cells = m_mesh.getCellCount(), vertices = m_mesh.getVertexCount();  // collective
        if (isRoot())
        {
          m_csv.open(m_cfg.output + ".csv");
          m_csv << "t,iterations,maxU,maxT,minPhi,maxPhi,vaporVolume,bulkTout,dp,pIn,maxTw,meanTw\n";
          Alert::Info() << "[mesh] cells=" << cells << " vertices=" << vertices
                        << "  bottom area=" << m_bottomArea << " m^2"
                        << "  Allen-Cahn eps=" << m_eps << " m gamma=" << m_gamma << " m/s"
                        << " neps=" << m_normalEps << " 1/m"
                        << Alert::Raise;
          Alert::Info() << "[phase change] "
                        << (m_cfg.interfaceTransfer ? "interfacial Hertz-Knudsen, F/dT = " : "bulk Lee")
                        << (m_cfg.interfaceTransfer ? std::to_string(m_hk / m_cfg.Tsat) + " kg/(m^2 s K)" : "")
                        << "  sites=" << m_sites.size() << " seed radius=" << m_seedRadius << " m"
                        << "  T_sat=" << m_cfg.Tsat << " K  p_out=" << m_cfg.outletPressure << " Pa"
                        << Alert::Raise;
          if (m_cfg.topView)
            Alert::Info() << "[top view] fluid=" << m_cfg.fluid << "  W=" << m_cfg.inletWidth << " -> "
                          << m_cfg.outletWidth << " m  d=" << m_cfg.depth << " m  C_f=" << m_cfg.darcy
                          << "  Nu=" << m_cfg.nusselt << " heated walls=" << m_cfg.heatedWalls
                          << "  film=" << m_cfg.filmThickness << " m"
                          << "  T_in=" << m_cfg.inletTemperature << " K" << Alert::Raise;
        }
      }

      int run()
      {
        // dt may change at the first nucleation (dtBoiling), so the loop runs on t.
        int n = 0;
        while (m_t < m_cfg.tEnd - 0.5 * m_cfg.dt)
        {
          ++n;
          int iterations = 0;
          bool retried = false;
          for (int attempt = 0;; ++attempt)
          {
            m_t += m_cfg.dt;
            bool ok = false, threw = false;
            try
            {
              iterations = step();
              ok = m_converged;
            }
            catch (const std::runtime_error& e)
            {
              threw = true;
              if (isRoot())
                Alert::Warning() << e.what() << Alert::Raise;
            }
            if (ok || attempt >= m_cfg.maxRetries || 0.5 * m_cfg.dt < m_cfg.dtMin)
            {
              if (threw)  // out of retries with a diverged solve: stop
                throw std::runtime_error("step failed at t = " + std::to_string(m_t)
                                         + " after " + std::to_string(attempt) + " retries");
              break;      // unconverged Picard after all retries: accept with a warning
            }
            // Reject: back to t^n, half the step.
            m_t -= m_cfg.dt;
            m_uIt.setData(m_uOld.getData());
            m_phiIt.setData(m_phiOld.getData());
            m_TIt.setData(m_TOld.getData());
            m_TwIt.setData(m_TwOld.getData());
            setTimeStep(0.5 * m_cfg.dt, "retry");
            retried = true;
          }

          m_uOld.setData(m_u.getSolution().getData());
          m_phiOld.setData(m_phi.getSolution().getData());
          m_TOld.setData(m_T.getSolution().getData());
          if (m_cfg.topView)
            m_TwOld.setData(m_Tw.getSolution().getData());
          if (!m_sites.empty() && n % m_cfg.nucleationEvery == 0)
            nucleate();           // seeds go into phi^n before n, kappa, Pi_h
          updateStabilisation();  // Pi_h at t^n for the next step
          if (m_cfg.lagInterface)
            updateInterface();    // n, kappa at t^n for the next step

          report(n, iterations);
          if (m_boiling && m_cfg.adaptiveDt)
            adaptTimeStep(iterations, retried);
          const bool last = !(m_t < m_cfg.tEnd - 0.5 * m_cfg.dt);
          const int every = (m_boiling && m_cfg.outputEveryBoiling > 0) ? m_cfg.outputEveryBoiling : m_cfg.outputEvery;
          if (n % every == 0 || last)
          {
            m_phiOut = m_phi.getSolution();
            m_TOut = m_T.getSolution();
            m_TsatOut.project(RealFunction([this](const Point& p) { return tsatAt(p); }));
            if (m_cfg.topView)
              m_TwOut = m_Tw.getSolution();
            m_xdmf.write(m_t).flush();
          }
        }
        m_xdmf.close();
        return 0;
      }

    private:
      /// The forms hold 1/dt: changing dt rebuilds them.
      void setTimeStep(Real dt, const char* why)
      {
        m_cfg.dt = dt;
        setupFlow();
        setupPhase();
        if (isRoot())
          Alert::Info() << "[dt] " << why << ": dt -> " << dt << " s" << Alert::Raise;
      }

      /// Picard-driven controller: grow after an easy step, shrink after a hard
      /// one, never grow right after a rejected attempt.
      void adaptTimeStep(int iterations, bool retried)
      {
        Real target = m_cfg.dt;
        if (iterations > m_cfg.picardHard || !m_converged)
          target = m_cfg.dtShrink * m_cfg.dt;
        else if (iterations <= m_cfg.picardEasy && !retried)
          target = m_cfg.dtGrow * m_cfg.dt;
        target = std::clamp(target, m_cfg.dtMin, m_dtMax);
        if (target != m_cfg.dt)
          setTimeStep(target, target < m_cfg.dt ? "Picard hard" : "Picard easy");
      }

      /// Nucleation sites: seed phi = max(phi, phi_seed) at every site whose
      /// wall superheat exceeds the onset value with no vapour nearby. The
      /// seed is the equilibrium profile of the Allen-Cahn term,
      /// phi_seed = (1 - tanh((r - R)/(2 eps)))/2, r = distance to the site.
      void nucleate()
      {
        const size_t dim = m_mesh.getSpaceDimension();
        const Real R = m_seedRadius, eps = m_eps;
        const auto& L = m_cfg.labels;
        const auto distance = [dim](const Point& p, Real xs) {
          Real r2 = (p.x() - xs) * (p.x() - xs);
          for (size_t i = 1; i < dim; ++i)
            r2 += p(i) * p(i);  // site at y = 0 (2D) or (y, z) = (0, 0) (3D)
          return std::sqrt(r2);
        };
        const Real seedVolume = (dim == 2) ? (m_cfg.topView ? 1.0 : 0.5) * M_PI * R * R : 2.0 / 3.0 * M_PI * R * R * R;

        std::vector<Real> active;
        for (size_t i = 0; i < m_sites.size(); ++i)
        {
          const Real xs = m_sites[i];
          if (m_t - m_lastSeed[i] < m_cfg.siteWaiting)
            continue;
          const RealFunction window([=](const Point& p) { return distance(p, xs) <= R ? 1.0 : 0.0; });
          // Wall temperature: T on the heated bottom (side view, 3D) or the
          // floor temperature T_w around the site (top view, sites on y = 0).
          const Real area = m_cfg.topView ? volumeIntegral(window) : boundaryIntegral(window, L.bottom);
          if (area <= 0.0)
            continue;
          const Real Tw = (m_cfg.topView ? volumeIntegral(m_TwIt * window)
                                         : boundaryIntegral(m_TIt * window, L.bottom)) / area;
          // The nucleus must sit in superheated liquid: a seed poking into the
          // subcooled core condenses at once and the implosion drives large
          // negative pressures (seen at 1 K wall superheat, 300 K inlet).
          // T_sat at the local pressure.
          const RealFunction cap([=](const Point& p) { return distance(p, xs) <= R ? 1.0 : 0.0; });
          const Real capVolume = volumeIntegral(cap);
          if (capVolume <= 0.0)
            continue;
          const Real Ts = saturation(volumeIntegral(m_pIt * cap) / capVolume).first;
          if (Tw - Ts < m_cfg.onsetSuperheat || volumeIntegral(m_TIt * cap) / capVolume < Ts)
            continue;
          const RealFunction ball([=](const Point& p) { return distance(p, xs) <= 2.0 * R ? 1.0 : 0.0; });
          if (volumeIntegral(m_phiIt * ball) > 0.05 * seedVolume)
            continue;
          active.push_back(xs);
          m_lastSeed[i] = m_t;
          if (isRoot())
            Alert::Info() << "[nucleation] t=" << m_t << " s  x=" << xs << " m  T_wall=" << Tw
                          << " K  T_sat(p)=" << Ts << " K" << Alert::Raise;
        }
        if (active.empty())
          return;

        ScalarGF2 seeded(m_sh2);
        seeded.project(RealFunction([&](const Point& p) {
          Real f = m_phiIt.getValue(p);
          for (const Real xs : active)
            f = std::max(f, 0.5 * (1.0 - std::tanh((distance(p, xs) - R) / (2.0 * eps))));
          return f; }));
        m_phiIt.setData(seeded.getData());
        m_phiOld.setData(seeded.getData());
        m_phi.getSolution().setData(seeded.getData());

        if (!m_boiling)
        {
          m_dtMax = (m_cfg.dtBoiling > 0.0) ? std::min(m_cfg.dtBoiling, m_cfg.dt) : m_cfg.dt;
          if (m_dtMax < m_cfg.dt)
            setTimeStep(m_dtMax, "nucleation");
        }
        m_boiling = true;
      }

      /// One time step: Picard iterations phase -> flow until the relative
      /// increments of (phi, T, u) fall below the tolerance.
      int step()
      {
        int k = 0;
        for (; k < m_cfg.maxIterations; ++k)
        {
          // (phi, T) monolithic, with properties, mdot and Pi_h at the iterate.
          updateCoefficients();
          Real dTw = 0.0;
          if (m_cfg.topView)
          {
            solve(m_phaseTop, m_phaseTopKSP);
            dTw = increment(m_Tw.getSolution(), m_TwIt);
            m_TwIt.setData(m_Tw.getSolution().getData());
          }
          else
            solve(m_phase, m_phaseKSP);
          const Real dPhi = increment(m_phi.getSolution(), m_phiIt);
          const Real dT = std::max(increment(m_T.getSolution(), m_TIt), dTw);
          m_phiIt.setData(m_phi.getSolution().getData());
          m_TIt.setData(m_T.getSolution().getData());

          // (u, p) with the fresh fraction: density, viscosity, capillary force
          // and continuity source see (phi, T)^{k+1}. The OSS projections stay
          // at t^n: lagged at the iterate they make Picard contract only like
          // tau_M/(tau_M + dt) on smooth pressure modes (~0.37 measured).
          if (!m_cfg.lagInterface)
            updateInterface();
          updateCoefficients();
          solve(m_flow, m_flowKSP);
          const Real dU = increment(m_u.getSolution(), m_uIt);
          m_uIt.setData(m_u.getSolution().getData());
          m_pIt.setData(m_p.getSolution().getData());

          if (m_cfg.picardMonitor && isRoot())
            Alert::Info() << "  Picard " << k + 1 << ": dPhi=" << dPhi << " dT=" << dT
                          << " dU=" << dU << Alert::Raise;

          if (std::max({ dPhi, dT, dU }) < m_cfg.tolerance)
          {
            m_converged = true;
            return k + 1;
          }
        }
        m_converged = false;
        if (isRoot())
          Alert::Warning() << "Picard did not converge at t = " << m_t << Alert::Raise;
        return k;
      }

      /// Relative L2 increment ||new - old|| / ||new|| (collective).
      template <class GF>
      static Real increment(const GF& current, const GF& previous)
      {
        GF diff(current.getFiniteElementSpace());
        diff.setData(current.getData());
        PetscReal scale = 0.0, delta = 0.0;
        PetscErrorCode ierr = VecAXPY(diff.getData(), -1.0, previous.getData());
        assert(ierr == PETSC_SUCCESS);
        ierr = VecNorm(current.getData(), NORM_2, &scale);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecNorm(diff.getData(), NORM_2, &delta);
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
        return delta / (scale > 0.0 ? scale : 1.0);
      }

      /// OSS projections Pi_h of the residuals. Called once per time step, at
      /// the converged state, so every Picard iteration of the next step sees
      /// Pi_h at t^n.
      void updateStabilisation()
      {
        m_piConv.project(Mult(Jacobian(m_uIt), m_uIt));
        m_piGradP.project(Grad(m_pIt));
        m_piDiv.project(Div(m_uIt));
        m_piConvPhi.project(Dot(m_uIt, Grad(m_phiIt)));
        m_piConvT.project(Dot(m_uIt, Grad(m_TIt)));
      }

      /// Coefficients of the forms as finite element functions of known
      /// polynomial order, refreshed from (phi^k, T^k, u^k): mixture properties
      /// and Lee terms in P2 (the space of phi, T), tau in P0. Rodin integrates
      /// an integrand of unknown order (a lambda) with the rule of the test
      /// space alone in a linear form but of trial + test in a bilinear one;
      /// on tetrahedra a P1 linear form then gets the one-point centroid rule.
      /// The lagged halves of the time derivatives and of the OSS terms,
      /// (rho u^n, v) and tau (Pi_h, B(v)), were thus integrated differently
      /// from their implicit halves (rho u, v) and tau (A(u), B(v)): every
      /// step applied u <- M_exact^{-1} M_centroid u^n, independent of dt.
      void updateCoefficients()
      {
        const auto project = [](auto& gf, auto&& f) { gf.project(RealFunction(f)); };
        project(m_rhoGF, [this](const Point& p) { return rhoAt(p); });
        project(m_muGF, [this](const Point& p) { return muAt(p); });
        project(m_kGF, [this](const Point& p) { return kAt(p); });
        project(m_rhoCpGF, [this](const Point& p) { return rhoCpAt(p); });
        project(m_phiClipGF, [this](const Point& p) { return phiAt(p); });
        project(m_mdotGF, [this](const Point& p) { return mdotAt(p); });
        project(m_cTGF, [this](const Point& p) { return dmdotdTAt(p); });
        project(m_cPhiGF, [this](const Point& p) { return dmdotdphiAt(p); });
        // Flow tau with its time scale floored at m_tauDt (the boiling dt):
        // with 2/dt, tau_M -> dt/2 as dt -> 0 removes the pressure
        // stabilisation tau_M/rho (grad p, grad q) that equal-order P1/P1
        // needs, and the pressure grows like 1/dt (every step retried with
        // dt/2 made it worse: 5.8e5 -> 3.8e8 Pa between dt = 1.2e-7 and 1e-8 s).
        // Without any time term (tau ~ h^2/(4 nu) ~ 1 ms in the liquid) the
        // lagged projection Pi[grad p^n] penalises every change of grad p with
        // weight tau/dt >> 1 and suppresses the inertial pressure.
        project(m_tauMGF, [this](const Point& p) { return tau(p, muAt(p) / rhoAt(p), darcyRateAt(p), true); });
        project(m_tauMRhoGF, [this](const Point& p) {
          return tau(p, muAt(p) / rhoAt(p), darcyRateAt(p), true) / rhoAt(p); });
        project(m_tauTGF, [this](const Point& p) {
          const Real C = rhoCpAt(p);
          return tau(p, kAt(p) / C, (m_cfg.hLV * dmdotdTAt(p) + wallExchangeAt(p) / m_cfg.depth) / C); });
        project(m_mdotTotGF, [this](const Point& p) { return mdotAt(p) + filmAt(p); });
        if (m_cfg.dcT > 0.0)
        {
          const auto gradT = Grad(m_TIt);
          project(m_dcTGF, [this, &gradT](const Point& p) {
            const auto g = gradT.getValue(p);
            const auto u = m_uIt.getValue(p);
            const Real gn = std::sqrt(Math::dot(g, g)), un = std::sqrt(Math::dot(u, u));
            if (gn <= 0.0 || un <= 0.0)
              return 0.0;
            const Real h = cellSize(p);
            const Real r = std::abs(Math::dot(u, g) - m_piConvT.get().getValue(p));
            return rhoCpAt(p) * std::min(m_cfg.dcT * 0.5 * h * r / gn, 0.5 * h * un); });
        }
        project(m_betaGF, [this](const Point& p) { return betaAt(p); });
        project(m_srcGF, [this](const Point& p) { return phiSourceGain(p) * (mdotAt(p) + filmAt(p)); });
        if (m_cfg.topView)
        {
          project(m_hwGF, [this](const Point& p) { return wallExchangeAt(p) / m_cfg.depth; });
          project(m_hfGF, [this](const Point& p) { return filmConductanceAt(p) / m_cfg.depth; });
          project(m_srcTwGF, [this](const Point& p) { return phiSourceGain(p) * dfilmdTwAt(p); });
          project(m_TsatGF, [this](const Point& p) { return tsatAt(p); });
        }
        project(m_srcTGF, [this](const Point& p) { return phiSourceGain(p) * dmdotdTAt(p); });
        project(m_srcPhiGF, [this](const Point& p) { return phiSourcePhi(p); });
        project(m_tauPhiGF, [this](const Point& p) {
          return tau(p, m_gamma * m_eps, -phiSourcePhi(p)); });
        project(m_tauCGF, [this](const Point& p) {
          return muAt(p) + 0.5 * rhoAt(p) * speed(p) * cellSize(p); });
      }

      /// Interface normal and curvature at the current iterate phi^k (physics,
      /// not stabilisation: refreshed at every Picard iteration).
      void updateInterface()
      {
        const auto g = Grad(m_phiIt);
        const Real eps = m_normalEps;
        const auto normal = (1.0 / Sqrt(Dot(g, g) + eps * eps)) * g;

        m_piNormal.project(normal);
        m_piKappa.project(-1.0 * Div(m_piNormal.get()));
      }

      static MeshType makeMesh(const Context::MPI& context, const Config& cfg)
      {
        const auto& comm = context.getCommunicator();
        Rodin::MPI::Sharder sharder(context);
        if (comm.rank() == 0)
        {
          Mesh<Context::Local> mesh;
          if (cfg.topView)
            mesh = makeChannelTop(cfg);
          else if (cfg.twoD)
            mesh = makeChannel2D(cfg);
          else
          {
            mesh.load(cfg.meshPath, IO::FileFormat::MEDIT);
            if (mesh.getSpaceDimension() != 3 || mesh.getDimension() != 3)
              throw std::runtime_error("ForcedConvection expects a tetrahedral mesh (or -fc_2d).");
            connect(mesh);
          }
#ifdef RODIN_USE_SCOTCH
          Scotch::Partitioner partitioner(mesh);
#else
          BalancedCompactPartitioner partitioner(mesh);
#endif
          partitioner.partition(static_cast<size_t>(comm.size()));
          sharder.shard(partitioner);
          sharder.scatter(0);
        }
        MeshType mesh = sharder.gather(0);
        connect(mesh);
        if (mesh.getDimension() == 3)
          mesh.reconcile(2);
        mesh.reconcile(1);
        return mesh;
      }

      /// [0, length] x [0, depth] triangulated, facets tagged by position.
      static Mesh<Context::Local> makeChannel2D(const Config& cfg)
      {
        auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, { cfg.nx + 1, cfg.ny + 1 });
        const Real hx = cfg.length / static_cast<Real>(cfg.nx);
        const Real hy = cfg.depth / static_cast<Real>(cfg.ny);
        for (auto it = mesh.getVertex(); it; ++it)
        {
          Math::SpatialPoint x = mesh.getVertexCoordinates(it->getIndex());
          x(0) *= hx;
          x(1) *= hy;
          mesh.setVertexCoordinates(it->getIndex(), x);
        }
        mesh.flush();
        connect(mesh);

        const Real tol = 1e-3 * std::min(hx, hy);
        const auto& L = cfg.labels;
        std::vector<std::pair<Index, Attribute>> tags;
        for (auto it = mesh.getBoundary(); it; ++it)
        {
          Math::SpatialPoint c(2);
          c.setZero();
          for (const auto& v : it->getVertices())
            c += mesh.getVertexCoordinates(v);
          c /= static_cast<Real>(it->getVertices().size());
          Attribute a = L.top;
          if (c(0) < tol)
            a = L.inlet;
          else if (c(0) > cfg.length - tol)
            a = L.outlet;
          else if (c(1) < tol)
            a = L.bottom;
          tags.emplace_back(it->getIndex(), a);
        }
        for (const auto& [f, a] : tags)
          mesh.setAttribute({ 1, f }, a);
        return mesh;
      }

      /// Top view: [0, length] x [-W(x)/2, W(x)/2], W linear from inletWidth
      /// to outletWidth; facets tagged inlet (x = 0), outlet (x = length), sides.
      static Mesh<Context::Local> makeChannelTop(const Config& cfg)
      {
        auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, { cfg.nx + 1, cfg.ny + 1 });
        const Real hx = cfg.length / static_cast<Real>(cfg.nx);
        for (auto it = mesh.getVertex(); it; ++it)
        {
          Math::SpatialPoint x = mesh.getVertexCoordinates(it->getIndex());
          const Real s = x(0) / static_cast<Real>(cfg.nx), eta = x(1) / static_cast<Real>(cfg.ny);
          const Real W = cfg.inletWidth + (cfg.outletWidth - cfg.inletWidth) * s;
          x(0) *= hx;
          x(1) = (eta - 0.5) * W;
          mesh.setVertexCoordinates(it->getIndex(), x);
        }
        mesh.flush();
        connect(mesh);

        const Real tol = 1e-3 * std::min(hx, cfg.inletWidth / static_cast<Real>(cfg.ny));
        const auto& L = cfg.labels;
        std::vector<std::pair<Index, Attribute>> tags;
        for (auto it = mesh.getBoundary(); it; ++it)
        {
          Math::SpatialPoint c(2);
          c.setZero();
          for (const auto& v : it->getVertices())
            c += mesh.getVertexCoordinates(v);
          c /= static_cast<Real>(it->getVertices().size());
          Attribute a = L.sides;
          if (c(0) < tol)
            a = L.inlet;
          else if (c(0) > cfg.length - tol)
            a = L.outlet;
          tags.emplace_back(it->getIndex(), a);
        }
        for (const auto& [f, a] : tags)
          mesh.setAttribute({ 1, f }, a);
        return mesh;
      }

      template <class M>
      static void connect(M& mesh)
      {
        const size_t D = mesh.getDimension();
        mesh.getConnectivity().compute(D, D);
        mesh.getConnectivity().compute(D, 0);
        mesh.getConnectivity().compute(D, D - 1);
        mesh.getConnectivity().compute(D - 1, D);
        mesh.getConnectivity().compute(D - 1, 0);
        if (D == 3)
          mesh.getConnectivity().compute(D - 1, 1);
        mesh.getConnectivity().compute(1, 0);
      }

      static Real cellSize(const Point& p)
      {
        return std::pow(p.getPolytope().getMeasure(), 1.0 / p.getPolytope().getDimension());
      }

      Real speed(const Point& p) const
      {
        const auto u = m_uIt.getValue(p);
        return std::sqrt(Math::dot(u, u));
      }

      // ---- Mixture properties and Lee coefficients at the last iterate ----
      // Evaluated at quadrature points from phi^k, T^k (P2), never from nodal
      // averages. phi is clipped only for property evaluation.

      Real phiAt(const Point& p) const
      {
        return std::min(std::max(m_phiIt.getValue(p), Real(0)), Real(1));
      }

      Real rhoAt(const Point& p) const
      {
        const Real f = phiAt(p);
        return f * m_cfg.rhoV + (1.0 - f) * m_cfg.rhoL;
      }

      Real muAt(const Point& p) const
      {
        const Real f = phiAt(p);
        return f * m_cfg.muV + (1.0 - f) * m_cfg.muL;
      }

      /// Harmonic mean: correct normal heat flux across the interface.
      Real kAt(const Point& p) const
      {
        const Real f = phiAt(p);
        return 1.0 / (f / m_cfg.kV + (1.0 - f) / m_cfg.kL);
      }

      Real rhoCpAt(const Point& p) const
      {
        const Real f = phiAt(p);
        return f * m_cfg.rhoV * m_cfg.cpV + (1.0 - f) * m_cfg.rhoL * m_cfg.cpL;
      }

      // ---- Lee source at the iterate (phi^k clipped, T^k) ----
      // x = T^k - T_sat, d = dTsmooth:
      //   [x]_+ = max(x, 0) + d log(1 + exp(-|x|/d))  (softplus, >= 0)
      //   [x]_- = x - [x]_+                            (<= 0)
      //   d[x]_+/dx = logistic(x/d)

      /// {T_sat(p), dT_sat/dp} (Clausius-Clapeyron through (p_out, Tsat)); p
      /// clamped to [0.05, 100] p_out, with zero slope outside.
      std::pair<Real, Real> saturation(Real p) const
      {
        if (!m_cfg.clapeyron)
          return { m_cfg.Tsat, 0.0 };
        const Real p0 = m_cfg.outletPressure, Rs = 8.314462618 / m_cfg.molarMass;
        const Real pc = std::clamp(p, 0.05 * p0, 100.0 * p0);
        const Real Ts = 1.0 / (1.0 / m_cfg.Tsat - Rs / m_cfg.hLV * std::log(pc / p0));
        return { Ts, (pc == p) ? Ts * Ts * Rs / (m_cfg.hLV * pc) : 0.0 };
      }

      Real tsatAt(const Point& p) const { return saturation(m_pIt.getValue(p)).first; }
      Real dtsatdpAt(const Point& p) const { return saturation(m_pIt.getValue(p)).second; }

      Real superheatAt(const Point& p) const  ///< [theta]_+
      {
        const Real x = m_TIt.getValue(p) - tsatAt(p), d = m_cfg.dTsmooth;
        return (std::max(x, Real(0)) + d * std::log1p(std::exp(-std::abs(x) / d))) / m_cfg.Tsat;
      }

      Real subcoolingAt(const Point& p) const  ///< [theta]_-
      {
        return (m_TIt.getValue(p) - tsatAt(p)) / m_cfg.Tsat - superheatAt(p);
      }

      Real logisticAt(const Point& p) const  ///< d[x]_+/dx in [0, 1]
      {
        const Real z = (m_TIt.getValue(p) - tsatAt(p)) / m_cfg.dTsmooth;
        if (z >= 0.0)
          return 1.0 / (1.0 + std::exp(-z));
        const Real e = std::exp(z);
        return e / (1.0 + e);
      }

      /// Interfacial area density A_i(phi) = s(phi) phi (1 - phi)/eps and its
      /// derivative; s ramps (C^1 smoothstep) from 0 at phi = c to 1 at 2c.
      std::pair<Real, Real> areaDensity(Real f) const
      {
        const auto [s, ds] = cutoff(f);
        const Real a = f * (1.0 - f) / m_eps, da = (1.0 - 2.0 * f) / m_eps;
        return { s * a, ds * a + s * da };
      }

      /// C^1 smoothstep s(phi): 0 for phi <= c, 1 for phi >= 2c, and s'.
      std::pair<Real, Real> cutoff(Real f) const
      {
        const Real c = m_cfg.areaCutoff;
        if (c <= 0.0)
          return { 1.0, 0.0 };
        const Real x = std::clamp((f - c) / c, Real(0), Real(1));
        return { x * x * (3.0 - 2.0 * x), (f > c && f < 2.0 * c) ? 6.0 * x * (1.0 - x) / c : 0.0 };
      }

      // ---- Top view: depth friction, floor exchange, evaporating film ----
      // All zero outside the top view.

      Real widthAt(Real x) const
      {
        const Real s = std::clamp(x / m_cfg.length, Real(0), Real(1));
        return m_cfg.inletWidth + (m_cfg.outletWidth - m_cfg.inletWidth) * s;
      }

      /// C_f nu/d^2 (1/s), reaction rate of the depth friction (enters tau_M).
      Real darcyRateAt(const Point& p) const
      {
        if (!m_cfg.topView)
          return 0.0;
        return m_cfg.darcy * muAt(p) / (rhoAt(p) * m_cfg.depth * m_cfg.depth);
      }

      /// Heated perimeter over width, P_h/W (exchange area per plan area).
      Real heatedPerimeterRatio(Real W) const
      {
        return (m_cfg.heatedWalls >= 3) ? (W + 2.0 * m_cfg.depth) / W : 1.0;
      }

      /// h_w P_h/W, h_w = Nu k/D_h (W/(m^2 K) per unit plan area), k
      /// phi-weighted: the liquid exchanges with the walls, a vapour slug
      /// barely (its film is h_f).
      Real wallExchangeAt(const Point& p) const
      {
        if (!m_cfg.topView)
          return 0.0;
        const Real W = widthAt(p.x()), d = m_cfg.depth, f = phiAt(p);
        return m_cfg.nusselt * ((1.0 - f) * m_cfg.kL + f * m_cfg.kV) * (W + d) / (2.0 * W * d)
             * heatedPerimeterRatio(W);
      }

      /// Film fraction F(phi) = s(phi) phi of the floor and dF/dphi: the cut-off
      /// keeps round-off vapour (phi ~ 1e-6) from growing at h_f/(h_lv d rho_v)
      /// ~ 5e3 1/s per K of wall superheat.
      std::pair<Real, Real> filmFraction(Real f) const
      {
        const auto [s, ds] = cutoff(f);
        return { s * f, ds * f + s };
      }

      /// h_f P_h/W, h_f = F(phi) k_l/delta_f (film on every heated wall).
      Real filmConductanceAt(const Point& p) const
      {
        if (!m_cfg.topView)
          return 0.0;
        return filmFraction(phiAt(p)).first * m_cfg.kL / m_cfg.filmThickness
             * heatedPerimeterRatio(widthAt(p.x()));
      }

      /// mdot_f = h_f (T_w - T_sat)/(h_lv d), kg/(m^3 s), at the iterate.
      Real filmAt(const Point& p) const
      {
        if (!m_cfg.topView)
          return 0.0;
        return filmConductanceAt(p) * (m_TwIt.getValue(p) - tsatAt(p)) / (m_cfg.hLV * m_cfg.depth);
      }

      Real dfilmdTwAt(const Point& p) const  ///< >= 0
      {
        if (!m_cfg.topView)
          return 0.0;
        return filmConductanceAt(p) / (m_cfg.hLV * m_cfg.depth);
      }

      Real dfilmdphiAt(const Point& p) const
      {
        if (!m_cfg.topView)
          return 0.0;
        return filmFraction(phiAt(p)).second * m_cfg.kL / m_cfg.filmThickness
             * heatedPerimeterRatio(widthAt(p.x()))
             * (m_TwIt.getValue(p) - tsatAt(p)) / (m_cfg.hLV * m_cfg.depth);
      }

      Real thetaAt(const Point& p) const  ///< (T - T_sat(p^k))/T_sat
      {
        return (m_TIt.getValue(p) - tsatAt(p)) / m_cfg.Tsat;
      }

      /// beta = (1/rho_v - 1/rho_l) d(mdot + mdot_f)/d(-p) >= 0: every source
      /// depends on T - T_sat(p) (or T_w - T_sat(p)).
      Real betaAt(const Point& p) const
      {
        const Real expansion = 1.0 / m_cfg.rhoV - 1.0 / m_cfg.rhoL;
        return expansion * (dmdotdTAt(p) + dfilmdTwAt(p)) * dtsatdpAt(p);
      }

      /// Interfacial: mdot = F_HK A_i(phi) theta (evaporation and condensation).
      /// Bulk Lee:    mdot = r_l (1 - phi) rho_l [theta]_+ + r_v phi rho_v [theta]_-
      Real mdotAt(const Point& p) const
      {
        const Real f = phiAt(p);
        if (m_cfg.interfaceTransfer)
          return m_hk * areaDensity(f).first * thetaAt(p);
        return m_cfg.rL * (1.0 - f) * m_cfg.rhoL * superheatAt(p)
             + m_cfg.rV * f * m_cfg.rhoV * subcoolingAt(p);
      }

      /// d mdot / d phi, implicit part only (<= 0, keeps the phi equation
      /// M-matrix-like); the rest is lagged through the residual form.
      Real dmdotdphiAt(const Point& p) const
      {
        if (m_cfg.interfaceTransfer)
          return std::min(m_hk * areaDensity(phiAt(p)).second * thetaAt(p), Real(0));
        return -m_cfg.rL * m_cfg.rhoL * superheatAt(p) + m_cfg.rV * m_cfg.rhoV * subcoolingAt(p);
      }

      /// d mdot / d T >= 0
      Real dmdotdTAt(const Point& p) const
      {
        const Real f = phiAt(p);
        if (m_cfg.interfaceTransfer)
          return m_hk * areaDensity(f).first / m_cfg.Tsat;
        const Real s = logisticAt(p);
        return (m_cfg.rL * (1.0 - f) * m_cfg.rhoL * s + m_cfg.rV * f * m_cfg.rhoV * (1.0 - s))
             / m_cfg.Tsat;
      }

      /// g(phi) = (1 - phi)/rho_v + phi/rho_l: volume created per kg evaporated
      /// in the phi equation consistent with continuity.
      Real phiSourceGain(const Point& p) const
      {
        const Real f = phiAt(p);
        return (1.0 - f) / m_cfg.rhoV + f / m_cfg.rhoL;
      }

      /// d(g mdot)/d phi, implicit part only (<= 0).
      Real phiSourcePhi(const Point& p) const
      {
        const Real expansion = 1.0 / m_cfg.rhoV - 1.0 / m_cfg.rhoL;
        return std::min(phiSourceGain(p) * (dmdotdphiAt(p) + dfilmdphiAt(p))
                        - expansion * (mdotAt(p) + filmAt(p)), Real(0));
      }

      /// tau = [(2/dt)^2 + (2|u^k|/h)^2 + (4 D/h^2)^2 + s^2]^{-1/2}
      /// with D a diffusivity and s a (nonnegative) reaction rate. For the flow
      /// the time scale is floored at m_tauDt (see updateCoefficients).
      Real tau(const Point& p, Real diffusivity, Real reaction, bool flow = false) const
      {
        const Real h = cellSize(p);
        const Real a = 2.0 / (flow ? std::max(m_cfg.dt, m_tauDt) : m_cfg.dt);
        const Real b = 2.0 * speed(p) / h, c = 4.0 * diffusivity / (h * h);
        return 1.0 / std::sqrt(a * a + b * b + c * c + reaction * reaction);
      }

      void setupFlow()
      {
        const size_t dim = m_mesh.getSpaceDimension();
        const Real dt = m_cfg.dt;
        const auto& L = m_cfg.labels;
        const auto n = BoundaryNormal(m_mesh);

        const auto convU = Mult(Jacobian(m_u), m_uIt);    // (grad u) u^k
        const auto streamV = Mult(Jacobian(m_v), m_uIt);  // (grad v) u^k
        const auto symU = 0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));
        const auto symV = 0.5 * (Jacobian(m_v) + Transpose(Jacobian(m_v)));
        const auto backflow = 0.5 * m_rhoGF * Max(-Dot(m_uIt, n), 0.0);
        const RealFunction pOut = m_cfg.outletPressure;

        // CSF: sigma kappa^k grad phi^k, kappa^k = -div Pi_h[n^k].
        const auto capillary = m_cfg.sigma * m_piKappa.get() * Grad(m_phiIt);

        // Mean G/rho_l, smooth ramp in time. 3D: product of parabolas in (y, z)
        // (peak 9/4 of the mean); 2D: parabola in y across the depth (peak 3/2).
        const auto inflow = VectorFunction(dim, [this, dim](const Point& p) {
          const Real W = m_cfg.inletWidth, d = m_cfg.depth;
          const Real s = std::min(m_t / m_cfg.ramp, 1.0);
          const Real ramp = 0.5 * (1.0 - std::cos(M_PI * s));
          const Real G = m_cfg.massFlux / m_cfg.rhoL;
          Math::SpatialVector<Real> u(dim);
          for (size_t i = 0; i < dim; ++i)
            u(i) = 0.0;
          if (m_cfg.topView)  // parabola across the inlet width (depth-averaged)
          {
            const Real eta = 2.0 * p.y() / W;
            u(0) = std::max(1.5 * G * ramp * (1.0 - eta * eta), 0.0);
          }
          else if (dim == 2)
          {
            const Real zeta = 2.0 * p.y() / d - 1.0;
            u(0) = std::max(1.5 * G * ramp * (1.0 - zeta * zeta), 0.0);
          }
          else
          {
            const Real eta = 2.0 * p.y() / W, zeta = 2.0 * p.z() / d - 1.0;
            u(0) = std::max(2.25 * G * ramp * (1.0 - eta * eta) * (1.0 - zeta * zeta), 0.0);
          }
          return u;
        });

        // Galerkin with rho^k, mu^k; the convective form is not
        // skew-symmetrised because div u = mdot (1/rho_v - 1/rho_l) != 0.
        // Split OSS: tau_M rho ((grad u) u^k - Pi[(grad u^k) u^k], (grad v) u^k)
        //          + tau_M/rho (grad p - Pi[grad p^k], grad q)
        //          + tau_C (div u - Pi[div u^k], div v).
        const Real expansion = 1.0 / m_cfg.rhoV - 1.0 / m_cfg.rhoL;
        // Top view: depth friction C_f mu/d^2 u (0 otherwise).
        const Real darcy = m_cfg.topView ? m_cfg.darcy / (m_cfg.depth * m_cfg.depth) : 0.0;
        m_flow = (1.0 / dt) * Integral(m_rhoGF * m_u, m_v) - (1.0 / dt) * Integral(m_rhoGF * m_uOld, m_v)
               + Integral(m_rhoGF * convU, m_v)
               + 2.0 * Integral(m_muGF * symU, symV)
               + darcy * Integral(m_muGF * m_u, m_v)
               - Integral(m_p, Div(m_v))
               + Integral(Div(m_u), m_q) - expansion * Integral(m_mdotTotGF, m_q)
               + Integral(m_betaGF * m_p, m_q) - Integral(m_betaGF * m_pIt, m_q)
               - Integral(m_rhoGF * m_g, m_v)
               - Integral(capillary, m_v)

               + Integral(m_tauMGF * m_rhoGF * convU, streamV)
               - Integral(m_tauMGF * m_rhoGF * m_piConv.get(), streamV)

               + Integral(m_tauMRhoGF * Grad(m_p), Grad(m_q))
               - Integral(m_tauMRhoGF * m_piGradP.get(), Grad(m_q))

               + Integral(m_tauCGF * Div(m_u), Div(m_v))
               - Integral(m_tauCGF * m_piDiv.get(), Div(m_v))

               + BoundaryIntegral(pOut * Dot(m_v, n)).over(L.outlet)
               + BoundaryIntegral(backflow * Dot(m_u, m_v)).over(L.outlet)

               + DirichletBC(m_u, Zero(dim)).on(FlatSet<Attribute>{ L.bottom, L.top, L.sides })
               + DirichletBC(m_u, inflow).on(L.inlet);
      }

      void setupPhase()
      {
        const Real dt = m_cfg.dt, rhoV = m_cfg.rhoV, hLV = m_cfg.hLV, Cw = m_cfg.wallCapacity;
        const Real Kw = m_cfg.wallConductance;
        const auto nb = BoundaryNormal(m_mesh);
        const auto& L = m_cfg.labels;
        const RealFunction q = m_cfg.heatFlux;
        const RealFunction Tin = m_cfg.inletTemperature;
        const RealFunction liquid = 0.0;
        const auto streamA = Dot(m_uIt, Grad(m_a));
        const auto streamW = Dot(m_uIt, Grad(m_w));
        const Real gammaEps = m_gamma * m_eps;
        const auto compression = m_gamma * ((1.0 - m_phiClipGF) * m_piNormal.get());  // (1 - phi^k) n^k

        // Linearised Lee source about (phi^k, T^k), in residual form so that at
        // convergence the source is exactly mdot^k (same coefficient, same rule):
        //   mdot ~= mdot^k + cT (T - T^k) + cPhi (phi - phi^k_clipped),
        //   cT = dmdot/dT >= 0,  cPhi = dmdot/dphi <= 0.
        // phi-equation in the form consistent with continuity: phi div u is
        // replaced by its exact value phi mdot (1/rho_v - 1/rho_l), so
        //   d phi/dt + u^k.grad phi = mdot g(phi),  g = (1 - phi)/rho_v + phi/rho_l.
        // The source vanishes (to O(1/rho_l)) at phi = 1 and phi stays bounded
        // whatever div u_h of the previous flow iterate is; with
        // phi div u^k the mismatch between div u^k (old mdot) and the new
        // mdot is amplified by 1/rho_v ~ 22 m^3/kg and phi overshot to ~2.
        // Linearised in residual form: S = S^k + g cT (T - T^k) + a (phi - phi^k),
        // a = min(g cPhi - (1/rho_v - 1/rho_l) mdot^k, 0) (implicit part only).
        // T-equation:   rho cp (T - T^n)/dt + rho cp u^k.grad T - div(k grad T)
        //               = -mdot h_lv.
        // Allen-Cahn (weak): gamma eps (grad phi, grad a) - gamma ((1 - phi^k) n^k phi, grad a).
        // Split OSS on the convective operators only.
        if (m_cfg.topView)
        {
          // Same phi and T equations plus the film source g mdot_f (implicit in
          // T_w, residual form), the floor exchange (h_w/d)(T - T_w), and the
          // wall equation divided by d (all three rows per unit volume):
          //   C_w/d dT_w/dt - div(K_w/d grad T_w) + h_w/d (T_w - T)
          //     + h_f/d (T_w - T_sat) = q''/d.
          const Real d = m_cfg.depth;
          const RealFunction qd = m_cfg.heatFlux / d;
          m_phaseTop = (1.0 / dt) * Integral(m_phi, m_a) - (1.0 / dt) * Integral(m_phiOld, m_a)
                     + Integral(Dot(m_uIt, Grad(m_phi)), m_a)
                     - Integral(m_srcTGF * m_T, m_a)
                     - Integral(m_srcTwGF * m_Tw, m_a)
                     - Integral(m_srcPhiGF * m_phi, m_a)
                     - Integral(m_srcGF - m_srcTGF * m_TIt - m_srcTwGF * m_TwIt - m_srcPhiGF * m_phiClipGF, m_a)
                     + gammaEps * Integral(Grad(m_phi), Grad(m_a))
                     - Integral(compression * m_phi, Grad(m_a))

                     + Integral(m_tauPhiGF * Dot(m_uIt, Grad(m_phi)), streamA)
                     - Integral(m_tauPhiGF * m_piConvPhi.get(), streamA)

                     + (1.0 / dt) * Integral(m_rhoCpGF * m_T, m_w) - (1.0 / dt) * Integral(m_rhoCpGF * m_TOld, m_w)
                     + Integral(m_rhoCpGF * Dot(m_uIt, Grad(m_T)), m_w)
                     + Integral(m_kGF * Grad(m_T), Grad(m_w))
                     + Integral(m_dcTGF * Grad(m_T), Grad(m_w))
                     + Integral(m_hwGF * m_T, m_w)
                     - Integral(m_hwGF * m_Tw, m_w)
                     + hLV * Integral(m_cTGF * m_T, m_w)
                     + hLV * Integral(m_cPhiGF * m_phi, m_w)
                     + hLV * Integral(m_mdotGF - m_cTGF * m_TIt - m_cPhiGF * m_phiClipGF, m_w)

                     + Integral(m_tauTGF * m_rhoCpGF * Dot(m_uIt, Grad(m_T)), streamW)
                     - Integral(m_tauTGF * m_rhoCpGF * m_piConvT.get(), streamW)

                     + (Cw / (d * dt)) * Integral(m_Tw, m_b) - (Cw / (d * dt)) * Integral(m_TwOld, m_b)
                     + (Kw / d) * Integral(Grad(m_Tw), Grad(m_b))
                     + Integral(m_hwGF * m_Tw, m_b)
                     - Integral(m_hwGF * m_T, m_b)
                     + Integral(m_hfGF * m_Tw, m_b)
                     - Integral(m_hfGF * m_TsatGF, m_b)
                     - Integral(qd, m_b)

                     + DirichletBC(m_T, Tin).on(L.inlet)
                     + DirichletBC(m_phi, liquid).on(L.inlet);
          return;
        }
        m_phase = (1.0 / dt) * Integral(m_phi, m_a) - (1.0 / dt) * Integral(m_phiOld, m_a)
                + Integral(Dot(m_uIt, Grad(m_phi)), m_a)
                - Integral(m_srcTGF * m_T, m_a)
                - Integral(m_srcPhiGF * m_phi, m_a)
                - Integral(m_srcGF - m_srcTGF * m_TIt - m_srcPhiGF * m_phiClipGF, m_a)
                + gammaEps * Integral(Grad(m_phi), Grad(m_a))
                - Integral(compression * m_phi, Grad(m_a))

                + Integral(m_tauPhiGF * Dot(m_uIt, Grad(m_phi)), streamA)
                - Integral(m_tauPhiGF * m_piConvPhi.get(), streamA)

                + (1.0 / dt) * Integral(m_rhoCpGF * m_T, m_w) - (1.0 / dt) * Integral(m_rhoCpGF * m_TOld, m_w)
                + Integral(m_rhoCpGF * Dot(m_uIt, Grad(m_T)), m_w)
                + Integral(m_kGF * Grad(m_T), Grad(m_w))
                + Integral(m_dcTGF * Grad(m_T), Grad(m_w))
                - BoundaryIntegral(q * m_w).over(L.bottom)
                + (Cw / dt) * BoundaryIntegral(m_T, m_w).over(L.bottom)
                - (Cw / dt) * BoundaryIntegral(m_TOld * m_w).over(L.bottom)
                + Kw * BoundaryIntegral(Grad(m_T), Grad(m_w)).over(L.bottom)
                - Kw * BoundaryIntegral(Dot(Grad(m_T), nb), Dot(Grad(m_w), nb)).over(L.bottom)
                + hLV * Integral(m_cTGF * m_T, m_w)
                + hLV * Integral(m_cPhiGF * m_phi, m_w)
                + hLV * Integral(m_mdotGF - m_cTGF * m_TIt - m_cPhiGF * m_phiClipGF, m_w)

                + Integral(m_tauTGF * m_rhoCpGF * Dot(m_uIt, Grad(m_T)), streamW)
                - Integral(m_tauTGF * m_rhoCpGF * m_piConvT.get(), streamW)

                + DirichletBC(m_T, Tin).on(L.inlet)
                + DirichletBC(m_phi, liquid).on(L.inlet);
      }

      template <class ProblemType>
      void solve(ProblemType& problem, Solver::KSP& ksp)
      {
        problem.assemble();
        problem.solve(ksp);
        ::KSPConvergedReason reason;
        PetscErrorCode ierr = KSPGetConvergedReason(ksp.getHandle(), &reason);
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
        if (reason < 0)
          throw std::runtime_error("KSP diverged at t = " + std::to_string(m_t));
      }

      template <class Expression>
      Real boundaryIntegral(const Expression& f, Attribute a)
      {
        m_flux = BoundaryIntegral(f, m_z).over(a);
        m_flux.assemble();
        return m_flux(m_one);
      }

      template <class Expression>
      Real volumeIntegral(const Expression& f)
      {
        m_flux = Integral(f, m_z);
        m_flux.assemble();
        return m_flux(m_one);
      }

      /// Outlet bulk temperature, inlet-outlet mean pressure drop, vapour
      /// volume and bounds of phi.
      void report(int step, int iterations)
      {
        const auto& L = m_cfg.labels;
        const auto n = BoundaryNormal(m_mesh);
        const auto& u = m_u.getSolution();
        const auto& T = m_T.getSolution();
        const auto& p = m_p.getSolution();
        const auto& phi = m_phi.getSolution();

        const Real flow = boundaryIntegral(Dot(u, n), L.outlet);
        const Real enthalpy = boundaryIntegral(T * Dot(u, n), L.outlet);
        const Real pIn = boundaryIntegral(p, L.inlet) / m_inletArea;
        const Real dp = pIn - boundaryIntegral(p, L.outlet) / m_outletArea;
        const Real bulk = (flow > 0.0) ? enthalpy / flow : m_cfg.inletTemperature;
        const Real vapour = volumeIntegral(phi);
        const Real umax = std::max(std::abs(u.max()), std::abs(u.min()));
        const Real Tmax = T.max();
        m_TsatOut.project(RealFunction([this](const Point& q) { return tsatAt(q); }));
        const Real TsatMax = m_TsatOut.max(), pMax = p.max(), pMin = p.min();
        const Real phiMin = phi.min(), phiMax = phi.max();
        // Floor temperature: T_w (top view) or T on the heated bottom.
        Real TwMax = Tmax, TwMean = 0.0;
        if (m_cfg.topView)
        {
          const auto& Tw = m_Tw.getSolution();
          TwMax = Tw.max();
          TwMean = volumeIntegral(Tw) / m_bottomArea;
        }
        else if (m_bottomArea > 0.0)
          TwMean = boundaryIntegral(T, L.bottom) / m_bottomArea;

        if (!isRoot())
          return;
        m_csv << m_t << ',' << iterations << ',' << umax << ',' << Tmax << ',' << phiMin << ','
              << phiMax << ',' << vapour << ',' << bulk << ',' << dp << ',' << pIn << ','
              << TwMax << ',' << TwMean << '\n';
        m_csv.flush();
        if (step % m_cfg.printEvery == 0)
          Alert::Info() << "t=" << m_t << " s  dt=" << m_cfg.dt << "  it=" << iterations << "  max|u|=" << umax
                        << " m/s  maxT=" << Tmax << " K  phi in [" << phiMin << ", " << phiMax
                        << "]  Vv=" << vapour << " m^3  Tb,out=" << bulk << " K  p_in=" << pIn << " Pa  dp=" << dp
                        << " Pa  Tw max/mean=" << TwMax << "/" << TwMean << " K  p in [" << pMin << ", " << pMax
                        << "] Pa  Tsat max=" << TsatMax << " K" << Alert::Raise;
      }

      bool isRoot() const { return m_mesh.getContext().getCommunicator().rank() == 0; }

      Config m_cfg;
      MeshType m_mesh;
      IO::XDMF m_xdmf;
      VectorFES m_vh;
      ScalarFES m_sh;
      ScalarFES2 m_sh2;
      CellFES m_ch;

      // Flow unknowns (P1/P1) and phase unknowns (P2/P2).
      VectorTrial m_u;
      ScalarTrial m_p;
      VectorTest m_v;
      ScalarTest m_q;
      ScalarTrial2 m_phi;
      ScalarTrial2 m_T;
      ScalarTest2 m_a;
      ScalarTest2 m_w;
      ScalarTrial2 m_Tw;  ///< top view: floor (wall) temperature
      ScalarTest2 m_b;

      // Time level n and Picard iterate k.
      VectorGF m_uOld, m_uIt;
      ScalarGF m_pIt;
      ScalarGF2 m_phiOld, m_TOld, m_phiIt, m_TIt, m_TwOld, m_TwIt;
      ScalarGF m_phiOut, m_TOut, m_TwOut, m_TsatOut;  ///< P1 copies for XDMF
      VectorGF m_g;               ///< gravity, constant

      // Coefficients actually used in the forms (known polynomial order),
      // projected from the point evaluators below by updateCoefficients().
      ScalarGF2 m_rhoGF, m_muGF, m_kGF, m_rhoCpGF, m_phiClipGF, m_mdotGF, m_cTGF, m_cPhiGF;
      ScalarGF2 m_srcGF, m_srcTGF, m_srcPhiGF;  ///< phi source g mdot and its linearisation
      /// Top view: total mdot (interface + film), g dmdot_f/dT_w, h_w/d, h_f/d.
      ScalarGF2 m_mdotTotGF, m_srcTwGF, m_hwGF, m_hfGF;
      ScalarGF2 m_betaGF, m_TsatGF;  ///< Clapeyron: continuity beta, T_sat(p^k)
      CellGF m_tauMGF, m_tauMRhoGF, m_tauTGF, m_tauPhiGF, m_tauCGF;
      CellGF m_dcTGF;  ///< rho cp k_dc, discontinuity capturing on T


      // Orthogonal projections Pi_h at the iterate.
      OrthogonalProjection<VectorFES> m_piConv, m_piGradP, m_piNormal;
      OrthogonalProjection<ScalarFES> m_piDiv, m_piKappa;
      OrthogonalProjection<ScalarFES2> m_piConvPhi, m_piConvT;

      FlowProblem m_flow;
      Solver::KSP m_flowKSP;
      PhaseProblem m_phase;
      Solver::KSP m_phaseKSP;
      PhaseProblemTop m_phaseTop;
      Solver::KSP m_phaseTopKSP;

      ScalarTest m_z;
      ScalarGF m_one;
      LinearForm<ScalarFES, ::Vec> m_flux;

      Real m_t = 0.0;
      Real m_eps = 0.0, m_gamma = 0.0, m_normalEps = 0.0;
      Real m_hk = 0.0, m_seedRadius = 0.0;  ///< Hertz-Knudsen coefficient, seed radius
      std::vector<Real> m_sites, m_lastSeed;  ///< nucleation sites (x) and last seed time
      bool m_boiling = false;
      bool m_converged = true;
      Real m_dtMax = 0.0;
      Real m_tauDt = 0.0;  ///< floor of the flow-tau time scale
      Real m_bottomArea = 0.0, m_inletArea = 0.0, m_outletArea = 0.0;
      std::ofstream m_csv;
  };
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  const auto setDefault = [](const char* name, const char* value) {
    PetscBool set = PETSC_FALSE;
    PetscErrorCode ierr = PetscOptionsHasName(PETSC_NULLPTR, PETSC_NULLPTR, name, &set);
    if (ierr == PETSC_SUCCESS && !set)
      ierr = PetscOptionsSetValue(PETSC_NULLPTR, name, value);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
  };
  // Flow: direct (MUMPS). Phase block: GMRES + block Jacobi/ILU.
  // Projections: CG + Jacobi on the mass matrices.
  setDefault("-fc_flow_ksp_type", "preonly");
  setDefault("-fc_flow_pc_type", "lu");
  setDefault("-fc_flow_pc_factor_mat_solver_type", "mumps");
  // This MUMPS build faults in its distributed RHS scatter: keep the
  // right-hand side and solution centralized (as in the Heart examples).
  for (const char* prefix : { "-fc_flow_", "-fc_phase_" })
  {
    setDefault((std::string(prefix) + "mat_mumps_icntl_20").c_str(), "0");
    setDefault((std::string(prefix) + "mat_mumps_icntl_21").c_str(), "0");
  }
  setDefault("-fc_phase_ksp_type", "gmres");
  setDefault("-fc_phase_pc_type", "bjacobi");
  setDefault("-fc_phase_ksp_rtol", "1e-8");
  setDefault("-fc_phase_ksp_gmres_restart", "100");
  for (const char* prefix : { "-fc_proj_conv_", "-fc_proj_gradp_", "-fc_proj_normal_",
                              "-fc_proj_div_", "-fc_proj_kappa_", "-fc_proj_phi_", "-fc_proj_T_" })
  {
    setDefault((std::string(prefix) + "ksp_type").c_str(), "cg");
    setDefault((std::string(prefix) + "pc_type").c_str(), "jacobi");
    setDefault((std::string(prefix) + "ksp_rtol").c_str(), "1e-10");
  }

  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world(PETSC_COMM_WORLD, boost::mpi::comm_attach);
  Rodin::Context::MPI context(env, world);

  int status = 0;
  try
  {
    using Rodin::Examples::Heat::ForcedConvection;
    ForcedConvection::Config cfg;

    char buffer[512];
    PetscBool got = PETSC_FALSE;
    PetscOptionsGetString(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_mesh", buffer, sizeof(buffer), &got);
    if (got)
      cfg.meshPath = buffer;
    PetscBool twoD = PETSC_FALSE;
    PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_2d", &twoD, PETSC_NULLPTR);
    cfg.twoD = (twoD == PETSC_TRUE);
    PetscBool top = PETSC_FALSE;
    PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_topview", &top, PETSC_NULLPTR);
    cfg.topView = (top == PETSC_TRUE);
    if (cfg.topView)
    {
      cfg.nx = 600;  // dx = 50 um, 12 cells across W = 0.2-0.4 mm
      cfg.ny = 12;
    }
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_wi", &cfg.inletWidth, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_wo", &cfg.outletWidth, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_depth", &cfg.depth, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_darcy", &cfg.darcy, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_nu", &cfg.nusselt, PETSC_NULLPTR);
    {
      PetscInt walls = cfg.heatedWalls;
      PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_heated_walls", &walls, PETSC_NULLPTR);
      cfg.heatedWalls = static_cast<int>(walls);
    }
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_film", &cfg.filmThickness, PETSC_NULLPTR);
    PetscOptionsGetString(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_fluid", buffer, sizeof(buffer), &got);
    if (got)
    {
      // Saturated at T_sat = 310 K (CoolProp 8): p_sat, rho_l, rho_v, mu_l,
      // mu_v, k_l, k_v, cp_l, cp_v, h_lv, sigma, M. "water-paper" keeps the
      // defaults (Alugoju et al. 2020).
      struct Fluid { const char* name; double p, rl, rv, ml, mv, kl, kv, cl, cv, h, s, M; };
      static const Fluid fluids[] = {
        { "water",   6231.0,   993.3, 0.04366, 6.933e-4, 1.008e-5, 0.6242,  0.01928, 4179.0, 1927.0, 2.414e6, 0.07019,  0.018015 },
        { "ethanol", 15170.0,  774.8, 0.2736,  8.666e-4, 9.135e-6, 0.1612,  0.01635, 2532.0, 1487.0, 9.073e5, 0.02074,  0.04607 },
        { "r245fa",  225700.0, 1306.0, 12.67,  3.439e-4, 1.232e-5, 0.08844, 0.01682, 1347.0, 941.9,  1.842e5, 0.01212,  0.13404 },
        { "r134a",   933400.0, 1160.0, 45.79,  1.680e-4, 1.222e-5, 0.07607, 0.01508, 1481.0, 1118.0, 1.663e5, 0.006509, 0.10203 },
        { "r1234ze", 702900.0, 1123.0, 37.18,  1.624e-4, 1.301e-5, 0.07025, 0.01465, 1430.0, 1033.0, 1.575e5, 0.007332, 0.11404 },
      };
      cfg.fluid = buffer;
      if (cfg.fluid != "water-paper")
      {
        const Fluid* f = nullptr;
        for (const auto& candidate : fluids)
          if (cfg.fluid == candidate.name)
            f = &candidate;
        if (!f)
          throw std::runtime_error("unknown -fc_fluid " + cfg.fluid
                                   + " (water-paper, water, ethanol, r245fa, r134a, r1234ze)");
        cfg.Tsat = 310.0;
        cfg.outletPressure = f->p;
        cfg.rhoL = f->rl; cfg.rhoV = f->rv; cfg.muL = f->ml; cfg.muV = f->mv;
        cfg.kL = f->kl; cfg.kV = f->kv; cfg.cpL = f->cl; cfg.cpV = f->cv;
        cfg.hLV = f->h; cfg.sigma = f->s; cfg.molarMass = f->M;
        cfg.inletTemperature = cfg.Tsat - 10.0;
      }
    }
    PetscInt nx = static_cast<PetscInt>(cfg.nx), ny = static_cast<PetscInt>(cfg.ny);
    PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_nx", &nx, PETSC_NULLPTR);
    PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_ny", &ny, PETSC_NULLPTR);
    cfg.nx = static_cast<size_t>(nx);
    cfg.ny = static_cast<size_t>(ny);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_length", &cfg.length, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_dt", &cfg.dt, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_tend", &cfg.tEnd, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_G", &cfg.massFlux, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_q", &cfg.heatFlux, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_tin", &cfg.inletTemperature, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_g", &cfg.gravity, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_tol", &cfg.tolerance, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_eps", &cfg.epsilon, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_gamma", &cfg.gamma, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_neps", &cfg.normalEps, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_dtsmooth", &cfg.dTsmooth, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_wall_capacity", &cfg.wallCapacity, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_wall_conductance", &cfg.wallConductance, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_pout", &cfg.outletPressure, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_tsat", &cfg.Tsat, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_dc", &cfg.dcT, PETSC_NULLPTR);
    {
      PetscBool cc = cfg.clapeyron ? PETSC_TRUE : PETSC_FALSE;
      PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_clapeyron", &cc, PETSC_NULLPTR);
      cfg.clapeyron = (cc == PETSC_TRUE);
    }
    PetscBool interfacial = cfg.interfaceTransfer ? PETSC_TRUE : PETSC_FALSE;
    PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_interface_transfer", &interfacial, PETSC_NULLPTR);
    cfg.interfaceTransfer = (interfacial == PETSC_TRUE);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_area_cutoff", &cfg.areaCutoff, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_accommodation", &cfg.accommodation, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_site_spacing", &cfg.siteSpacing, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_onset", &cfg.onsetSuperheat, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_seed_radius", &cfg.seedRadius, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_site_wait", &cfg.siteWaiting, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_dt_boil", &cfg.dtBoiling, PETSC_NULLPTR);
    {
      PetscBool adaptive = cfg.adaptiveDt ? PETSC_TRUE : PETSC_FALSE;
      PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_adaptive_dt", &adaptive, PETSC_NULLPTR);
      cfg.adaptiveDt = (adaptive == PETSC_TRUE);
      PetscInt easy = cfg.picardEasy, hard = cfg.picardHard;
      PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_picard_easy", &easy, PETSC_NULLPTR);
      PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_picard_hard", &hard, PETSC_NULLPTR);
      cfg.picardEasy = static_cast<int>(easy);
      cfg.picardHard = static_cast<int>(hard);
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_dt_grow", &cfg.dtGrow, PETSC_NULLPTR);
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_dt_shrink", &cfg.dtShrink, PETSC_NULLPTR);
    }
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_dt_min", &cfg.dtMin, PETSC_NULLPTR);
    {
      PetscInt retries = cfg.maxRetries;
      PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_max_retries", &retries, PETSC_NULLPTR);
      cfg.maxRetries = static_cast<int>(retries);
    }
    {
      PetscInt every = cfg.nucleationEvery;
      PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_nucleation_every", &every, PETSC_NULLPTR);
      cfg.nucleationEvery = std::max(static_cast<int>(every), 1);
    }
    Rodin::Real r = cfg.rL;
    PetscBool gotR = PETSC_FALSE;
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_r", &r, &gotR);
    if (gotR)
      cfg.rL = cfg.rV = r;
    PetscInt every = cfg.outputEvery, maxit = cfg.maxIterations;
    PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_output_every", &every, PETSC_NULLPTR);
    {
      PetscInt everyBoil = cfg.outputEveryBoiling;
      PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_output_every_boil", &everyBoil, PETSC_NULLPTR);
      cfg.outputEveryBoiling = static_cast<int>(everyBoil);
    }
    PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_maxit", &maxit, PETSC_NULLPTR);
    PetscInt print = cfg.printEvery;
    PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_print_every", &print, PETSC_NULLPTR);
    cfg.printEvery = std::max(static_cast<int>(print), 1);
    PetscBool monitor = PETSC_FALSE;
    PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_picard_monitor", &monitor, PETSC_NULLPTR);
    cfg.picardMonitor = (monitor == PETSC_TRUE);
    PetscBool lag = cfg.lagInterface ? PETSC_TRUE : PETSC_FALSE;
    PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_lag_interface", &lag, PETSC_NULLPTR);
    cfg.lagInterface = (lag == PETSC_TRUE);
    cfg.outputEvery = static_cast<int>(every);
    cfg.maxIterations = static_cast<int>(maxit);

    ForcedConvection simulation(context, cfg);
    status = simulation.run();
  }
  catch (const std::exception& e)
  {
    std::cerr << "Fatal error: " << e.what() << "\n";
    status = 1;
  }
  PetscFinalize();
  return status;
}
