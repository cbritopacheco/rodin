/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file CoupledLV0DCoronary3D.h
 * @brief Coupled 0D left-ventricle / 3D coronary flow solver.
 */
#ifndef EXAMPLES_HEART_CORONARYARTERY_COUPLEDLV0DCORONARY3D_H
#define EXAMPLES_HEART_CORONARYARTERY_COUPLEDLV0DCORONARY3D_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "Rodin/Heart/CCMLC2014.h"
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/MPI.h>
#include <Rodin/PETSc.h>
#include <Rodin/Solver.h>
#include <Rodin/Types.h>
#include <Rodin/Variational.h>

#include "CoronaryArteryTiming.h"
#include "VMSConvectionIntegrator.h"

namespace Rodin::Examples::Heart
{
  class CoupledLV0DCoronary3D
  {
    public:
      using Real = Rodin::Real;
      using Model = Rodin::Heart::CCMLC2014T<>;
      using Attribute = Rodin::Geometry::Attribute;
      using MeshType = Rodin::Geometry::Mesh<Rodin::Context::MPI>;

      using VelocityFESType =
        Rodin::Variational::H1<2, Rodin::Math::SpatialVector<Real>, MeshType>;

      using PressureFESType = Rodin::Variational::H1<1, Real, MeshType>;

      using VelocityGridFunctionType =
        Rodin::PETSc::Variational::GridFunction<VelocityFESType>;

      using PressureGridFunctionType =
        Rodin::PETSc::Variational::GridFunction<PressureFESType>;

      using VelocityTrialFunctionType =
        Rodin::PETSc::Variational::TrialFunction<VelocityGridFunctionType,
          VelocityFESType>;

      using VelocityTestFunctionType =
        Rodin::PETSc::Variational::TestFunction<VelocityFESType>;

      using PressureTrialFunctionType =
        Rodin::PETSc::Variational::TrialFunction<PressureGridFunctionType,
          PressureFESType>;

      using PressureTestFunctionType =
        Rodin::PETSc::Variational::TestFunction<PressureFESType>;

      using FluxLinearFormType = Rodin::Variational::LinearForm<PressureFESType, ::Vec>;

      using LinearSystemType = Rodin::PETSc::Math::LinearSystem;

      using FlowProblemType =
        Rodin::Variational::Problem<LinearSystemType, VelocityTrialFunctionType,
          PressureTrialFunctionType, VelocityTestFunctionType, PressureTestFunctionType>;

      using ScalarProjectionProblemType = Rodin::Variational::Problem<LinearSystemType,
        PressureTrialFunctionType, PressureTestFunctionType>;

      using VelocityProjectionProblemType = Rodin::Variational::Problem<LinearSystemType,
        VelocityTrialFunctionType, VelocityTestFunctionType>;

      using TauFESType = Rodin::Variational::P1<Real, MeshType>;

      using TauGridFunctionType = Rodin::PETSc::Variational::GridFunction<TauFESType>;

      using TauTrialFunctionType =
        Rodin::PETSc::Variational::TrialFunction<TauGridFunctionType, TauFESType>;

      using TauTestFunctionType = Rodin::PETSc::Variational::TestFunction<TauFESType>;

      using TauProblemType = Rodin::Variational::Problem<LinearSystemType,
        TauTrialFunctionType, TauTestFunctionType>;

      using VMSFESType =
        Rodin::Variational::H1<2, Rodin::Math::SpatialVector<Real>, MeshType>;

      using VMSGridFunctionType = Rodin::PETSc::Variational::GridFunction<VMSFESType>;

      using VMSTrialFunctionType =
        Rodin::PETSc::Variational::TrialFunction<VMSGridFunctionType, VMSFESType>;

      using VMSTestFunctionType = Rodin::PETSc::Variational::TestFunction<VMSFESType>;

      using VMSProblemType = Rodin::Variational::Problem<LinearSystemType,
        VMSTrialFunctionType, VMSTestFunctionType>;

      /**
       * @brief State and calibrated constants of one coronary outlet.
       *
       * @details One state per outlet: the microvascular transmural pressure
       *          p_tm = p_c - p_im. The compartment is intramyocardial, so
       *          p_im is not a source added to the balance but the reference
       *          of the whole compartment; the venous limb is a Starling
       *          resistor whose throat closes at p_im (Permutt-Riley;
       *          Downey & Kirk 1975). That single constitutive statement makes
       *          p_tm a stable relaxation variable, bounded below by zero and
       *          above by R_v max(q_a), which is why no collapsible-tube law,
       *          no unstressed areas, no parallel multiplicities and no
       *          non-negativity penalty are needed: they all existed to bound
       *          a variable that is now bounded by construction.
       */
      struct RCR
      {
    /// @brief Microvascular transmural pressure p_tm = p_c - p_im. THE state.
          Real ptm = 0.0;
    /// @brief Diagnostic: microvascular pressure, p_c = p_tm + p_im.
          Real pc = 0.0;
    /// @brief Outlet pressure applied to the 3D model (the p_c + p_im part;
    ///        the R_a Q part is assembled implicitly, see Ra below).
          Real pout = 0.0;
    /// @brief Flow leaving the microvascular compartment towards the right
    ///        atrium through the Starling throat.
          Real qd = 0.0;
    /// @brief Intramyocardial (tissue) pressure, p_im = alpha p_LV.
          Real pim = 0.0;
    /// @brief Arteriolar (dominant) lumped resistance at the reference
    ///        viscosity. Assembled implicitly on the outlet boundary, so it is
    ///        never exposed to the one-step lag.
          Real Ra = 0.0;
    /// @brief Venular lumped resistance at the reference viscosity.
          Real Rv = 0.0;
    /// @brief Microvascular compliance, C_tot weighted by the Murray split.
          Real C = 0.0;
    /// @brief Measured outlet area. Only used to scale the implicit
    ///        resistance term, R_a A (u.n)(v.n).
          Real area = 0.0;
    /// @brief Calibrated branch flow, reference point of the rheological
    ///        modulation Phi(q).
          Real q0 = 0.0;
    /// @brief Derived resting wall shear rate of the arteriolar limb (1/s).
    /// @details g_a = r_a dP_a / (2 mu_N L_a). Not an input: it follows from
    ///          the calibre, the path length and the pressure budget.
          Real gammaA = 0.0;
    /// @brief Derived resting wall shear rate of the venular limb (1/s).
          Real gammaV = 0.0;
    /// @brief Predicted arteriolar bed multiplicity, N_a. Diagnostic.
          Real Na = 0.0;
    /// @brief Predicted venular bed multiplicity, N_v. Diagnostic.
          Real Nv = 0.0;
          /// @brief Diagnostic: apparent-viscosity multiplier of the arteriolar limb.
          Real muA = 1.0;
          /// @brief Diagnostic: apparent-viscosity multiplier of the venular limb.
          Real muV = 1.0;
          /// @brief Diagnostic: stored microvascular volume, V = C p_tm.
          Real vol = 0.0;
          /// @brief Diagnostic: resting stored volume.
          Real vol0 = 0.0;
      };

      /**
   * @brief Rheological model used by the 0D outlet closure.
   */
      enum class RheologyModel
      {
        /// @brief Carreau-Yasuda, mu = mu_inf + (mu_0-mu_inf)[1+(lambda g)^a]^((n-1)/a).
        CarreauYasuda,
        /// @brief Quemada, parameterized by haematocrit and plasma viscosity.
        Quemada
      };

      /**
   * @brief Quemada blood viscosity parameters.
   *
   * @details mu = mu_F [1 - k(g) phi / 2]^-2, with
   *          k(g) = (k_0 + k_inf sqrt(g/g_c)) / (1 + sqrt(g/g_c)).
   *
   *          Its value here is not a better fit: it is that it separates the
   *          two mechanisms that Carreau-Yasuda entangles. The haematocrit phi
   *          and the plasma viscosity mu_F set the *high-shear* level, where
   *          the coronary bed actually operates; k_0 sets the low-shear
   *          aggregation rise. A Carreau-Yasuda pair fitted to a healthy and a
   *          hyperviscous condition can share mu_inf -- and then differ only
   *          below ~50 1/s, which is a regime a normally perfused bed never
   *          visits (see the gamma_rest identity in the calibration). Quemada
   *          makes the haematocrit axis explicit, and along that axis the
   *          effect is first order: phi from 0.45 to 0.55 raises the apparent
   *          viscosity 1.9-fold at 400 1/s and costs 47 per cent of resting
   *          coronary flow.
   */
      struct Quemada
      {
          /// @brief Plasma viscosity (Pa s).
          /// @details The second pathological axis: paraproteinaemia raises it
          ///          without touching the haematocrit.
          Real plasmaViscosity = 0.0017963;
          /// @brief Haematocrit (volume fraction).
          Real hematocrit = 0.45;

          /// @brief Zero-shear intrinsic viscosity. Derived from phi if <= 0.
          /// @details k_0, k_inf and gamma_c are *not* constants: they are functions
          ///          of the haematocrit, and treating them as fixed breaks the
          ///          model. Quemada's law has a packing limit phi_max = 2/k, so a
          ///          k_0 frozen at its phi = 0.45 value makes the viscosity diverge
          ///          above phi = 0.46 -- exactly the range a polycythaemia study
          ///          needs. The Cokelet correlations are used when these are left
          ///          non-positive:
          ///
          ///            k_0     = exp( 3.874 - 10.41 phi + 13.8 phi^2 - 6.738 phi^3)
          ///            k_inf   = exp( 1.3435 - 2.803 phi + 2.711 phi^2 - 0.6479 phi^3)
          ///            gamma_c = exp(-6.1508 + 27.923 phi - 25.6 phi^2 + 3.697 phi^3)
          Real k0 = 3.7;
          /// @brief Infinite-shear intrinsic viscosity. Derived from phi if <= 0.
          Real kInf = 1.66;
          /// @brief Critical shear rate of the aggregation transition (1/s).
          ///        Derived from phi if <= 0.
          Real gammaC = 2.29;
      };

      /**
   * @brief Carreau-Yasuda blood viscosity parameters.
   */
      struct CarreauYasuda
      {
          /// @brief Low-shear viscosity.
          Real mu0 = 0.301;
          /// @brief Infinite-shear viscosity.
          Real muInf = 0.0055;
          /// @brief Relaxation time.
          Real lambda = 16.15;
          /// @brief Power-law index.
          Real n = 0.21;
          /// @brief Yasuda transition exponent.
          Real yasuda = 0.77;
          /// @brief Shear-rate regularization used in the 3D viscosity.
          Real gammaRegularization = 1.0e-3;
      };

      /**
   * @brief Tabulated WRMS closure and scalar-solve tolerances for the 0D
   *        outlet.
   *
   * @details The Weissenberg-Rabinowitsch-Mooney-Schofield closure is retained
   *          exactly. Writing the flow of a generalized-Newtonian fluid in a
   *          tube as Hagen-Poiseuille with an apparent viscosity,
   *
   *            Q = pi R^4 dp / (8 mu_ap L)   and   Q = pi R^3 I(tau_w)/tau_w^3
   *            =>  mu_ap(tau_w) = tau_w^4 / (4 I(tau_w)),
   *            with I(tau_w) = int_0^{tau_w} tau^2 gammadot(tau) dtau,
   *
   *          shows that mu_ap is a *universal* function of the wall shear
   *          stress for a given rheology: it does not depend on R, on L or on
   *          the branch. It is therefore tabulated once, in the nominal shear
   *          rate gammadot = 4Q/(pi R^3) = tau_w / mu_ap, and the 0D element
   *          only interpolates. With 241 log-spaced nodes the maximum log-log
   *          interpolation error is 0.09 per cent, against a per-call Newton
   *          plus a 100-step RK4 quadrature in the previous implementation.
   */
      struct OutletFlowLaw
      {
          /// @brief Lower bound of the tabulated wall shear stress (Pa).
          Real tableTauMin = 1.0e-6;
          /// @brief Upper bound of the tabulated wall shear stress (Pa).
          Real tableTauMax = 1.0e4;
          /// @brief Number of log-spaced table nodes.
          int tableNodes = 241;
          /// @brief Quadrature nodes used once per table entry to evaluate I(tau_w).
          int integralSteps = 2000;
          /// @brief Wall shear root solver step tolerance, in log(gammadot).
          Real shearStepTolerance = 1.0e-12;
          /// @brief Wall shear root solver maximum iterations.
          int shearMaxIterations = 200;
          /// @brief Newton tolerance on the outlet transmural pressure (Pa).
          Real outletStepTolerance = 1.0e-9;
          /// @brief Maximum Newton iterations of the outlet update.
          int outletMaxIterations = 50;
          /// @brief Flow magnitude below which the rheological modulation is frozen
          ///        at its low-shear plateau.
          Real zeroFlowTolerance = 1.0e-16;
          /// @brief Residual venous conductance kept when the waterfall is shut, as a
          ///        fraction of the open value.
          /// @details With a hard on/off switch the Newton Jacobian loses its only
          ///          restoring term the instant p_tm crosses zero, leaving a pure
          ///          integrator of the 3D flux: a transient backflow then drives
          ///          p_tm to large negative values with nothing to stop it. A small
          ///          leak keeps J bounded away from C/dt and makes the switch
          ///          differentiable enough for Newton. It is a regularisation, not a
          ///          physical shunt: at 1e-3 it moves the open-lumen solution by less
          ///          than 0.1 per cent.
          Real closedLumenLeak = 1.0e-3;
          /// @brief Enforce the lower bound p_tm >= 0 that the formulation assumes.
          /// @details The RCR closure is derived assuming p_tm is bounded below by
          ///          zero. That bound holds only while the venous lumen is open; it
          ///          is imposed here so it holds unconditionally.
          bool clampTransmuralPressure = true;
      };

      /**
   * @brief Universal WRMS apparent-viscosity table, log-log interpolated.
   *
   * @details Nodes are stored in the nominal shear rate
   *          gammadot = 4Q/(pi R^3), which is monotone in tau_w, so the lookup
   *          needed by the 0D element (given a flow, return the apparent
   *          viscosity) is a direct interpolation with no root finding.
   */
      struct WRMSTable
      {
          /// @brief log(gammadot) nodes, strictly increasing.
          std::vector<Real> logGamma;
          /// @brief log(mu_ap) at each node.
          std::vector<Real> logMu;

          /// @brief Apparent viscosity at a nominal shear rate. Clamped, not
          ///        extrapolated: outside the table the Carreau-Yasuda plateaus are
          ///        the correct continuation.
          Real operator()(Real gamma) const
          {
            if (logGamma.size() < 2)
              return std::exp(logMu.empty() ? 0.0 : logMu.front());

            const Real lg = std::log(std::max(gamma, 1e-300));

            if (lg <= logGamma.front())
              return std::exp(logMu.front());
            if (lg >= logGamma.back())
              return std::exp(logMu.back());

            const auto it = std::upper_bound(logGamma.begin(), logGamma.end(), lg);
            const std::size_t i = static_cast<std::size_t>(it - logGamma.begin()) - 1;
            const Real w = (lg - logGamma[i]) / (logGamma[i + 1] - logGamma[i]);

            return std::exp(logMu[i] + w * (logMu[i + 1] - logMu[i]));
          }
      };

      /**
   * @brief Piecewise-linear LV activation waveform parameters.
   */
      struct Activation
      {
          /// @brief Period of the prescribed cardiac cycle.
          Real period = 0.85;
          /// @brief Activation ramp start.
          Real tRampStart = 0.15;
          /// @brief Activation ramp end.
          Real tRampEnd = 0.21;
          /// @brief Activation plateau end.
          Real tPlateauEnd = 0.36;
          /// @brief Relaxation ramp end.
          Real tRelaxEnd = 0.45;
          /// @brief Negative activation plateau end.
          Real tNegativeEnd = 0.6;
          /// @brief Positive activation plateau value.
          Real positiveValue = 35.0;
          /// @brief Negative activation plateau value.
          Real negativeValue = -20.0;
      };

      /**
   * @brief Piecewise-linear atrial pressure waveform parameters.
   */
      struct AtrialPressure
      {
          /// @brief Period of the prescribed cardiac cycle.
          Real period = 0.85;
          /// @brief Baseline atrial pressure.
          Real minValue = 500.0;
          /// @brief First plateau atrial pressure.
          Real maxValue = 1000.0;
          /// @brief Second plateau atrial pressure.
          Real secondThreshold = 1250.0;
          /// @brief End of first down-ramp.
          Real t1 = 0.02;
          /// @brief End of first plateau.
          Real t2 = 0.15;
          /// @brief End of first return ramp.
          Real t3 = 0.17;
          /// @brief End of second up-ramp.
          Real t4 = 0.56;
          /// @brief End of second plateau.
          Real t5 = 0.62;
          /// @brief End of cycle return ramp.
          Real t6 = 0.85;
      };

      /**
   * @brief Parameters passed to the 0D LV model.
   */
      struct LVModel
      {
          /// @brief 0D fluid density.
          Real rho = 1.0e3;
          /// @brief Reference radius.
          Real R0 = 2.45e-2;
          /// @brief Reference wall thickness.
          Real d0 = 1.4e-2;
          /// @brief Passive elastic stiffness.

          Real Es = 3.0e6;
          /// @brief Viscous parameter.
          Real mu = 70.0;
          /// @brief Viscous parameter.
          Real eta = 70.0;
          /// @brief Active stress gain.
          Real alpha = 1.5;
          /// @brief Load-dependent relaxation time scale.
          Real alphaR = 0.12;
          /// @brief Active stiffness scale.
          Real k0 = 1.0e5;
          /// @brief Active stress scale.
          Real sigma0 = 1.2e5;
          /// @brief Proximal arterial resistance.
          Real Rp = 5e7;
          /// @brief Proximal arterial compliance.
          Real Cp = 5e-9;
          /// @brief Distal arterial resistance.
          Real Rd = 1.0e8;
          /// @brief Distal arterial compliance.
          Real Cd = 5.0e-10;
          /// @brief Atrial valve coefficient.
          Real Kat = 1.0e-6;
          /// @brief Peripheral valve coefficient.
          Real Kp = 5.0e-10;
          /// @brief Arterial valve coefficient.
          Real Kar = 2.e-7;
          /// @brief LV cavity capacity.
          Real cavityCapacity = 5.0e-12;
          /// @brief Local 0D Newton absolute tolerance.
          Real localTolerance = 1.0e-12;
          /// @brief Local 0D Newton maximum iterations.
          int localMaxIterations = 50;
          /// @brief Local 0D Newton damping.
          Real localDamping = 1.0;
          /// @brief Absolute-value regularization.
          Real absRegularization = 1.0e-14;
          /// @brief Initial fiber deformation.
          Real initFibDef = 0.0;
          /// @brief Initial active stiffness.
          Real initActiveStiffness = 0.0;
          /// @brief Initial active stress.
          Real initActiveStress = 0.0;
          /// @brief Low-fiber-deformation target for load-dependent relaxation.
          Real relaxationM0Low = 1.6;
          /// @brief High-fiber-deformation target for load-dependent relaxation.
          Real relaxationM0High = 1.0;
          /// @brief Fiber deformation at which m0 reaches relaxationM0Low.
          Real relaxationM0LowEc = 0.0;
          /// @brief Fiber deformation at which m0 reaches relaxationM0High.
          Real relaxationM0HighEc = 2.0;
          /// @brief Systemic venous pressure callback value.
          Real systemicVenousPressure = 1.0e3;
          /// @brief Passive energy parameter mu1.
          Real passiveMu1 = 0.0;
          /// @brief Passive energy parameter mu2.
          Real passiveMu2 = 0.0;
          /// @brief Passive energy parameter C0.
          Real passiveC0 = 1.9e3;
          /// @brief Passive energy parameter C1.
          Real passiveC1 = 1.1e-1;
          /// @brief Passive energy parameter C2.
          Real passiveC2 = 1.9e3;
          /// @brief Passive energy parameter C3.
          Real passiveC3 = 1.1e-1;
          /// @brief 0D model maximum Newton iterations.
          int maxIterations = 200;
          /// @brief 0D model absolute tolerance.
          Real absoluteTolerance = 1.0e-8;
          /// @brief 0D model relative tolerance.
          Real relativeTolerance = 1.0e-8;
          /// @brief 0D model step tolerance.
          Real stepTolerance = 1.0e-10;
          /// @brief 0D model damping factor.
          Real dampingFactor = 1.0;
          /// @brief Initial LV fiber deformation state.
          Real initialY = 0.0;
          /// @brief Initial LV velocity state.
          Real initialV = 0.0;
          /// @brief Offset applied to atrial pressure to initialize pv.
          Real initialPvOffset = -100.0;
          /// @brief Initial arterial pressure.
          Real initialPar = 11000.0;
          /// @brief Initial distal pressure.
          Real initialPd = 10000.0;
          /// @brief Low-shear viscosity.
          Real mu_0 = 0.301;
          /// @brief Infinite-shear viscosity.
          Real mu_Inf = 0.0055;
          /// @brief Relaxation time.
          Real lambda = 16.152;
          /// @brief Power-law index.
          Real n = 0.21;
          /// @brief Yasuda transition exponent.
          Real yasuda = 0.77;
          /// @brief Proximal surrogate vessel radius.
          Real proximalRadius = 0.0125;
          /// @brief Proximal surrogate vessel length.
          Real proximalLength = 0.4;
          /// @brief Distal surrogate vessel radius.
          Real distalRadius = 0.00175;
          /// @brief Distal surrogate vessel length.
          Real distalLength = 0.2;
      };

      /**
   * @brief Linearization strategy for the 3D coronary flow solve.
   */
      enum class FlowMode
      {
        /**
     * @brief Full Newton linearization of the nonlinear 3D flow residual.
     *
     * The convective term is differentiated with respect to the current
     * velocity iterate, and the Carreau-Yasuda viscosity includes its
     * directional derivative in the Jacobian.
     */
        Newton,

        /**
     * @brief Oseen/Picard linearization with lagged transport coefficients.
     *
     * The advecting velocity, Temam divergence factor, and
     * Carreau-Yasuda viscosity are evaluated using the previous time-step
     * velocity. This gives a lagged linearized 3D flow system that is
     * solved directly with PETSc KSP, without PETSc SNES nonlinear
     * iterations.
     */
        Oseen
      };

      struct Config
      {
          /// @brief Input coronary fluid mesh path.
          std::string meshPath = "CoronaryArtery.mesh";
          /// @brief Basename for XDMF and related output files.
          std::string xdmfBasename = "CoronaryArtery";
          /// @brief CSV diagnostics output path.
          std::string csvPath = "CoronaryArtery.csv";

          /// @brief No-slip wall boundary attribute.
          Attribute wall = 2;
          /// @brief Inlet boundary attribute.
          Attribute inlet = 4;

          /// @brief Outlet boundary attributes, in the same order used by RCR data.
          std::array<Attribute, 6> outlets{{7, 8, 9, 10, 14, 15}};

          /// @brief Mesh coordinate scale applied after partitioning.
          Real meshScale = 1.0e-3;
          /// @brief Pressure stabilization parameter.
          Real eps = 1.0e-12;
          /// @brief 3D blood density.
          Real rho = 1060.0;
          /// @brief Inlet reversed-flow damping multiplier. Set to 0 to disable.
          Real inletBackflowStabilization = 1.0;
          /// @brief Inlet normal impedance coefficient in Pa s / m. Set to 0 to
          /// disable.
          /// @details Defaults to defaultRCR.Rp times the scaled inlet area.
          Real inletImpedance = 5.e2;
          /// @brief Outlet backflow damping multiplier. Set to 0 to disable.
          Real outletBackflowStabilization = 1.0;

          /// @brief Time-step size.
          Real dt = 1.0e-3;
          /// @brief Number of time steps.
          size_t nsteps = 3 * static_cast<int>(0.85 / 1.0e-3);
          /// @brief Factor applied to dt when the 3D KSP/SNES solve fails.
          Real timeAdaptivityReductionFactor = 0.5;
          /// @brief Maximum number of successive dt reductions per accepted step.
          int timeAdaptivityMaxLevels = 8;

          /// @brief 3D coronary flow linearization mode. Defaults to Oseen/Picard.
          FlowMode flowMode = FlowMode::Oseen;
          /// @brief Blood viscosity model shared by 3D flow and outlet laws.
          CarreauYasuda viscosity;

          Real newtonianCalibrationViscosity = 0.0035;
          /// @brief Non-Newtonian outlet flow-law parameters.
          OutletFlowLaw outletFlowLaw;
          /// @brief Prescribed LV activation waveform parameters.
          Activation activation;
          /// @brief Prescribed atrial pressure waveform parameters.
          AtrialPressure atrialPressure;
          /// @brief 0D LV model parameters and initial conditions.
          LVModel lv;
          /// @brief Default RCR parameters copied to every outlet at startup.
          /// @details Designated initializers: robust against member reordering.
          /// @details p_c is the microvascular pressure, so it initialises near the
          ///          right atrial level. The calibration overwrites these with the
          ///          steady state of the calibrated network.
          RCR defaultRCR{.ptm = 1400.0, .pc = 2000.0, .pout = 2000.0};

          /// @brief Enable the anatomical outlet calibration at startup.
          /// @details Outlet areas are measured on the mesh and branch flows split as
          ///          Q_i proportional to r_i^3 (Murray's law) summing to
          ///          lcaTargetFlow. Disabling it leaves a single vessel per bed
          ///          and the configured geometry, for diagnostic use only.
          bool autoCalibrateOutlets = true;
          /// @brief Total target coronary inflow distributed across outlets (m^3/s).
          /// @details ~1.0e-6 m^3/s is about 60 mL/min; LCA rest flow ~150-250.
          Real lcaTargetFlow = 1.5e-6;

          /// @brief Derive the bed resistance from morphometry instead of
          ///        calibrating it to lcaTargetFlow.
          /// @details When false (the default) R_a = dP_a / Q_i with
          ///          Q_i = lcaTargetFlow * w: the resting flow is an input and
          ///          the resistance follows, so nothing upstream -- a stenosis
          ///          included -- can move the calibration. When true the
          ///          arteriolar multiplicity N_a is prescribed instead
          ///          (arteriolarCount, split by the same Murray weight) and the
          ///          flow becomes a prediction, Q_i = pi r_a^2 v_a N_a. The two
          ///          closures are algebraically identical when
          ///          arteriolarCount = lcaTargetFlow / (pi r_a^2 v_a), which is
          ///          the default below: switching this on changes nothing until
          ///          arteriolarCount is edited away from that value.
          bool morphometricResistance = false;

          /// @brief Terminal arteriolar count of the whole bed, split by Murray.
          /// @details Used only when morphometricResistance is true. The default
          ///          reproduces lcaTargetFlow at the default (r_a, v_a).
          ///          Coronary morphometry puts it at 1e5-1e6.
          Real arteriolarCount = 1.5279e5;

          /// @brief Microvascular vasodilation factor. 1 = rest.
          /// @details R_a and R_v are divided by it, so 4-5 reproduces the
          ///          adenosine hyperaemia under which FFR is measured. At rest
          ///          an epicardial stenosis below ~90 per cent area is not flow
          ///          limiting, because the microvascular bed carries almost all
          ///          the resistance; it only becomes limiting once the bed is
          ///          dilated. A stenosis study therefore needs this above one.
          Real vasodilationFactor = 1.0;

          /// @brief Total microvascular compliance (m^3/Pa), split across
          ///        outlets by the Murray weight.
          /// @details The identified intramyocardial compliance. Together with
          ///          R_v it predicts the run-off time constant tau = C R_v,
          ///          which is reported at calibration and is the quantity to
          ///          check against the diastolic decay. The model is
          ///          remarkably insensitive to it: over 1e-10 to 3e-9 the
          ///          systolic/diastolic inflow ratio moves only from 0.38 to
          ///          0.30, because the systolic impediment comes from the
          ///          tissue pressure entering p_out, not from the storage.
          Real coronaryComplianceTotal = 4e-10;

          /// @brief Rheological model of the 0D outlet closure.
          RheologyModel rheologyModel = RheologyModel::CarreauYasuda;
          /// @brief Quemada parameters, used when rheologyModel == Quemada.
          Quemada quemada;

          /// @brief Morphometric radius of the arteriolar limb (m).
          /// @details The prescribed pair is (r, v): both are read directly off
          ///          intravital microscopy, which measures a diameter and a
          ///          red-cell velocity. Everything else is derived,
          ///
          ///            g_0 = 4 v / r                       (Poiseuille wall shear)
          ///            N   = Q_i / (pi r^2 v)              (bed multiplicity)
          ///            L   = r dP_share / (2 mu_N g_0)     (pressure budget)
          ///            T   = L / v                         (transit time)
          ///
          ///          so the effective path length is an *output*. That is the
          ///          right way round: L is a lumped path through several
          ///          generations in series and is not a measurable quantity,
          ///          whereas r and v are.
          ///
          ///          Three measurables (r, v, T) for two degrees of freedom
          ///          leaves one consistency check: the derived transit time
          ///          must match the indicator-dilution value (~1-2 s at rest).
          ///          It is reported at calibration and warned about when it
          ///          falls outside 0.5-3 s. With r = 25 um and v = 5 mm/s the
          ///          arteriolar limb needs T = 8 s, four times the physiological
          ///          value: a single (r, L, N) cannot carry 87 per cent of the
          ///          pressure drop at that calibre and that velocity. The check
          ///          is a statement about the lumping, not a numerical fault.
          Real arteriolarRadius = 25.0e-6;
          /// @brief Mean arteriolar velocity (m/s). Physiological range 2-10 mm/s.
          Real arteriolarVelocity = 5.0e-3;
          /// @brief Morphometric radius of the venular limb (m).
          Real venularRadius = 30.0e-6;
          /// @brief Mean venular velocity (m/s). Physiological range 1-5 mm/s.
          Real venularVelocity = 3.0e-3;
          /// @brief Reference mean transit time used only by the consistency
          ///        check (s). Indicator dilution, ~1-2 s at rest.
          Real referenceTransitTime = 1.5;

          /// @brief Venular share of the microvascular pressure budget.
          /// @details The measured coronary distribution of head loss is about
          ///          75 per cent arteriolar/capillary and 10-15 per cent
          ///          venular. It is the only quantity that splits the pressure
          ///          budget: R_v = f_v dP / Q and R_a = (1 - f_v) dP / Q.
          ///          Nothing absorbs a residual.
          Real venularPressureFraction = 0.13;

          /// @brief Fraction of LV pressure transmitted to the intramyocardial
          ///        compartment, p_im = intramyocardialFraction * p_LV.
          /// @details Bed-averaged transmission. Tissue pressure runs from
          ///          about p_LV in the subendocardium to a small fraction of
          ///          it in the subepicardium, so a volume-weighted average
          ///          over the wall is of order 0.5-0.7. It is the single
          ///          constant of the tissue-pressure block: the
          ///          shortening-induced term (alpha_a tau_c) and the
          ///          first-order lag (tau_f) have both been removed. Neither
          ///          had an independent measurement, and neither is needed:
          ///          p_im is not a source added to the balance but the
          ///          reference of the whole compartment, so a steep p_LV
          ///          translates the operating point instead of injecting an
          ///          unbounded flux. The model *is* sensitive to this value
          ///          and the systolic/diastolic inflow ratio should be checked
          ///          against it.
          Real intramyocardialFraction = 0.7;
          /// @brief Right atrial (coronary sinus) drainage pressure.
          Real rightAtrialPressure = 1800.0;

          /// @brief Operating right atrial pressure of the running 0D outlets
          ///        (Pa). <= 0 means "same as rightAtrialPressure".
          /// @details Separates the runtime drainage pressure from the value
          ///          used by the resting calibration. A venous-hypertension
          ///          scenario (right-heart failure, tricuspid regurgitation,
          ///          constrictive pericarditis; ~1800-2500 Pa, 13-19 mmHg)
          ///          should displace the operating point of a bed whose
          ///          (R_a, R_v, C) were calibrated at the healthy baseline:
          ///          set THIS value and leave rightAtrialPressure alone.
          ///          Setting rightAtrialPressure itself instead models a
          ///          chronically adapted bed, because the calibration budget
          ///          dP = p_ar(0) - P_RA then absorbs part of the change.
          ///          The first cycles carry a small transient because the
          ///          calibrated initial p_tm uses the baseline value; discard
          ///          them as usual.
          Real operatingRightAtrialPressure = 0.0;

          /// @brief Freeze the reduced outlet closure at its high-shear
          ///        plateau, giving a CONSTANT outlet resistance.
          /// @details The comparison the paper reports at fixed haematological
          ///          state: R_mu is replaced by
          ///          R_inf = 8 mu_inf L / (N pi R^4), i.e. the WRMS table is
          ///          built with mu(gammadot) == mu_inf so Phi = mu_inf/mu_N is
          ///          constant and the pressure-flow relation of every outlet
          ///          becomes linear. The plateau is taken from the ACTIVE law
          ///          (mu_inf for Carreau-Yasuda; mu_F (1 - k_inf phi/2)^-2 for
          ///          Quemada), so setting this is exactly "same mu_inf, zero
          ///          shear-thinning".
          ///
          ///          Only the 0D closure changes: Config::viscosity still
          ///          drives the 3D Carreau-Yasuda field, so a run with this
          ///          flag differs from its R_mu twin ONLY in the reduced
          ///          outlet model. That is what isolates the closure, and it
          ///          is why the flag exists instead of setting
          ///          viscosity.mu0 = viscosity.muInf, which would make the
          ///          resolved domain Newtonian as well.
          bool constantOutletResistance = false;

          Real inletTangentialDamping = 1e3;
          Real inletVelocityDamping = 0.0;
      };

      explicit CoupledLV0DCoronary3D(const Rodin::Context::MPI& context);
      CoupledLV0DCoronary3D(const Rodin::Context::MPI& context, const Config& cfg);
      ~CoupledLV0DCoronary3D();

      CoupledLV0DCoronary3D(const CoupledLV0DCoronary3D&) = delete;
      CoupledLV0DCoronary3D& operator=(const CoupledLV0DCoronary3D&) = delete;

      int run();
      CoupledLV0DCoronary3D& initialize();

      Config& getConfig() noexcept
      {
        return m_cfg;
      }
      const Config& getConfig() const noexcept
      {
        return m_cfg;
      }

      Model& getModel() noexcept
      {
        return m_model;
      }
      const Model& getModel() const noexcept
      {
        return m_model;
      }

    private:
      struct StepData
      {
          Real t = 0.0;
          Real pat = 0.0;
          Real psv = 0.0;

          Real y = 0.0;
          Real v = 0.0;
          Real radius = 0.0;
          Real lvVolume = 0.0;
          Real lvFlow = 0.0;
          Real pv = 0.0;
          Real par = 0.0;
          Real pd = 0.0;

          Real ec = 0.0;
          Real gamma = 0.0;
          Real beta = 0.0;
          Real w = 1.0;
          Real kc = 0.0;
          Real tauc = 0.0;

          Real qIn = 0.0;
          Real qOutSum = 0.0;
          Real qDistalSum = 0.0;
          Real qCapChargingSum = 0.0;
          Real flowBalance = 0.0;

          std::map<Attribute, Real> qOut;
          std::map<Attribute, Real> qDistal;
          std::map<Attribute, Real> pc;
          std::map<Attribute, Real> pOut;
          Real pim = 0.0;

          /// @brief Mechanism diagnostics, summed or averaged over the outlets.
          /// @details p_tm is the state, the stored volume C p_tm is the pump,
          ///          and the two viscosity ratios are the rheological
          ///          modulation of each limb. There is no radius any more:
          ///          the conductance modulation is now the Starling throat,
          ///          which is a switch on the drainage pressure and not a
          ///          calibre.
          Real ptm = 0.0;
          Real muARatio = 1.0;
          Real muVRatio = 1.0;
          Real storedVolume = 0.0;
      };

      static Model::Input makeInput(const Config& cfg);
      static MeshType makeMesh(const Rodin::Context::MPI& context, const Config& cfg);

      /// @brief Advance one coronary outlet by one time step.
      /// @details Implicit Euler and a scalar Newton on the transmural
      ///          pressure p_tm. The residual derivative is
      ///          C/dt + dq_v/dp_tm > 0 everywhere, so the iteration converges
      ///          globally from any iterate: no line search, no nested solve,
      ///          no bracketing.
      static void updateOutlet0D(const Config& cfg, const WRMSTable& wrms,
        const Model& model, RCR& bc, Real Q, Real dt);

      /// @brief Build the universal WRMS apparent-viscosity table.
      /// @details mu_ap(tau_w) = tau_w^4 / (4 I(tau_w)) with
      ///          I(tau_w) = int_0^{tau_w} tau^2 gammadot(tau) dtau. Independent
      ///          of radius, length and branch, so it is built once.
      static WRMSTable buildWRMSTable(const Config& cfg, const CarreauYasuda& visc);

      static Real periodic_activation(const Activation& cfg, Real t);
      static Real atrial_pressure(const AtrialPressure& cfg, Real t);

      void setupModel();
      void setupMeshAndSpaces();
      void setupProjectionSolvers();
      void setupDiagnostics();
      void printInitialState() const;
      bool isRoot() const;

      bool advance0D();
      void solveStatic();
      bool solve3D();
      void computeFluxes();
      void computeWallShear();
      void updateHistory();
      void writeOutputs();
      void writeCSVHeader();
      void writeCSVRow();
      StepData collectStepData() const;
      void printStepTiming(int step) const;

      Config m_cfg;
      Model::Input m_input;
      Model m_model;

      MeshType m_mesh;
      Rodin::IO::XDMF m_xdmf;

      VelocityFESType m_uh;
      PressureFESType m_ph;
      VMSFESType m_uph;
      TauFESType m_tauh;

      VelocityTrialFunctionType m_u;
      PressureTrialFunctionType m_p;
      PressureTrialFunctionType m_mu;

      VelocityTestFunctionType m_v;
      PressureTestFunctionType m_q;
      PressureTestFunctionType m_r;

      VelocityGridFunctionType m_uOld;
      PressureGridFunctionType m_pOld;
      PressureGridFunctionType m_one;
      PressureTestFunctionType m_qFlux;

      FluxLinearFormType m_flux;

      /// @brief Recovered wall shear stress field (P2 velocity space), written to
      ///        XDMF as "shearStress".
      VelocityGridFunctionType m_shearWall;

      VMSTrialFunctionType m_sub;
      VMSGridFunctionType m_subOld;
      VMSTrialFunctionType m_up;
      VMSTestFunctionType m_vp;

      TauGridFunctionType m_tauOld;
      TauTrialFunctionType m_tau;
      TauTestFunctionType m_t;

      FlowProblemType m_flow;
      VMSProblemType m_l2ConvU;
      VMSProblemType m_subProjection;
      TauProblemType m_tauProjection;
      Rodin::Solver::KSP m_flowKSP;
      Rodin::Solver::SNES m_flowSolver;
      ScalarProjectionProblemType m_viscosityProjection;
      Rodin::Solver::KSP m_viscosityProjectionKSP;
      Rodin::Solver::KSP m_l2ConvUSolver;
      Rodin::Solver::KSP m_subProjectionSolver;
      Rodin::Solver::KSP m_tauProjectionSolver;

      // Wall-shear recovery projections. Hoisted to members so their mass
      // matrix and preconditioner are built once and reused across every
      // output step instead of being reconstructed (and, under a direct
      // solver, refactorized) each call to computeWallShear().
      VelocityTrialFunctionType m_gradRecTrial;
      VelocityTestFunctionType m_gradRecTest;
      VelocityProjectionProblemType m_gradRecProj;
      Rodin::Solver::KSP m_gradRecKSP;
      VelocityTrialFunctionType m_wssTrial;
      VelocityTestFunctionType m_wssTest;
      VelocityProjectionProblemType m_wssProj;
      Rodin::Solver::KSP m_wssKSP;

      std::map<Attribute, RCR> m_wk;

      /// @brief Universal WRMS apparent-viscosity table. Built once at
      ///        calibration and shared by every outlet and both limbs.
      WRMSTable m_wrms;

      mutable StepData m_stepData;
      std::ofstream m_csv;
      StepTiming m_stepTiming;
      bool m_flowFieldSplitsSet = false;
      bool m_initialized = false;
  };
} // namespace Rodin::Examples::Heart

#endif
