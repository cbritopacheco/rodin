/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * Explicit (staggered) sequential PETSc coronary ALE FSI -- P1/P1 VARIANT.
 *
 * Equal-order P1 velocity / P1 pressure fluid (NOT inf-sup stable) and P1
 * solid.  Because Taylor-Hood is lost, the pressure mode is stabilized with an
 * additional pressure-gradient (PSPG / grad-p) VMS term on top of the
 * convective and grad-div stabilizations.  Everything else matches the P2/P1
 * solver (CoronaryArtery_FSI_Explicit_PETSc_Seq.cpp).  k = 1 in the tau
 * formulas.  Outputs go to results_p1p1/.
 *
 * Physics: Carreau-Yasuda blood flow (ALE, conservative BDF1) in a coronary
 * tree, coupled to a hyperelastic arterial wall (total Lagrangian, Newmark)
 * by the Robin-Robin loose coupling of Burman, Durst, Fernandez, Guzman &
 * Ruz (2025), incl. the interface convective stabilization and the
 * alpha = gamma sqrt(rho_s E) scaling.  0D heart model (CCMLC2014) at the
 * inlet.  Outlets: the R-mu intramyocardial closure of CoupledLV0DCoronary3D
 * -- Starling resistor on the transmural pressure p_tm, universal WRMS
 * apparent-viscosity table (Carreau-Yasuda or Quemada), and the arteriolar
 * resistance R_a Phi_a assembled IMPLICITLY as R_a Phi_a A (u.n)(v.n) so the
 * 0D-3D coupling is unconditionally stable in dt.
 *
 * Stabilization parameters (tau_K, tau_C, tau_p) are evaluated POINTWISE from
 * the local lagged Carreau-Yasuda viscosity (Atrium pattern); only the
 * orthogonal-subscale projections Pi[(grad u)u], u' and pi~ are L2-projected.
 *
 * Startup: static follower-pressure prestress to par(0) (exact load
 * stiffness, quadratic Newton).  The dynamic wall load is split into a
 * follower-pressure LEVEL + projected transfer REMAINDER, so the first
 * dynamic residual vanishes identically: no ramps, no startup transient.
 *
 * Step n -> n+1 (per coupling iterate; 1 = loose):
 *   solid (SNES) -> harmonic ALE lift + mesh move -> fluid (Oseen, KSP);
 *   commit, interface fluxes, RCR update.
 *
 * Attributes: 1 = FSI wall; caps: fluid 27 (inlet) / 24,25,26,28,30,31
 * (outlets), solid 153 (inlet) / 150,151,152,154,155,156 (outlets);
 * 99 = clamped FSI ring band; 110..120 = solid heart-contact patches.
 */

#include "Rodin/Solid/Integrators/InternalVirtualWorkResidual.h"
#include "Rodin/Solid/Integrators/InternalVirtualWorkTangent.h"
#include "Rodin/Variational/BoundaryIntegral.h"
#include "Rodin/Variational/ForwardDecls.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <unordered_set>

#include <petscsys.h>
#include <petscvec.h>

#include <Rodin/Alert.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/PETSc.h>
#include <Rodin/Solid.h>
#include <Rodin/Solver.h>
#include <Rodin/Variational.h>

#include "Rodin/Heart/CCMLC2014.h"

// Projected-VMS convective stabilization (lagged Oseen), shared with the
// working CoupledLV0DCoronary3D fluid solver.
#include "CoronaryArtery/VMSConvectionIntegrator.h"
// Position-graded Yeoh law (transmural intima/media/adventitia profile).
#include "CoronaryArtery/GradedYeoh.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Math;
using namespace Rodin::Solver;
using namespace Rodin::Variational;
using namespace Rodin::Heart;

namespace
{
  using Model = Rodin::Heart::CCMLC2014T<>;
  using MeshType = Geometry::Mesh<Context::Local>;

  struct BoundaryFluid
  {
      static constexpr Attribute FSI = 2;
      static constexpr Attribute Inlet = 4;
      static constexpr std::array<Attribute, 6> Outlets{{7, 8, 9, 10, 14, 15}};
      static constexpr Attribute FSIRing = 99;
  };

  struct BoundarySolid
  {
      static constexpr Attribute FSI = 1;
      static constexpr Attribute Inlet = 150;
      static constexpr std::array<Attribute, 6> Outlets{{151, 152, 153, 154, 155, 156}};
      static constexpr Attribute FSIRing = 99;
      static constexpr std::array<Attribute, 1> Contact{{101}};
      /// Outer (adventitial) surface pieces; xi = 1 there for the transmural
      /// coordinate (xi = 0 on the FSI/luminal surface).
      static constexpr std::array<Attribute, 2> Outer{{2, 101}};
  };

  // --------------------------------------------------------------------------
  // R-mu coronary outlet (ported from CoupledLV0DCoronary3D): one state per
  // outlet, the microvascular transmural pressure p_tm = p_c - p_im.  The
  // compartment is intramyocardial, so p_im is the REFERENCE of the whole
  // compartment (not a source added to the balance); the venous limb is a
  // Starling resistor whose throat closes at p_im (Permutt-Riley; Downey &
  // Kirk 1975).  The arteriolar resistance R_a Phi_a is assembled IMPLICITLY
  // on the 3D outlet boundary as R_a Phi_a A (u.n)(v.n), so p_out carries only
  // p_c and the coupling is unconditionally stable in dt.  Phi = mu_ap/mu_N is
  // the rheological modulation read off the universal WRMS table.
  // See RCR_formulacion_minima.tex and CoupledLV0DCoronary3D.{h,cpp}.
  // --------------------------------------------------------------------------
  struct RCR
  {
      /// Microvascular transmural pressure p_tm = p_c - p_im.  THE state.
      Real ptm = 0.0;
      /// Diagnostic: microvascular pressure, p_c = p_tm + p_im.
      Real pc = 0.0;
      /// Outlet pressure applied to the 3D model (the p_c part; R_a Phi_a Q is
      /// assembled implicitly, see Ra below).
      Real pout = 0.0;
      /// Flow leaving the compartment towards the right atrium (Starling).
      Real qd = 0.0;
      /// Intramyocardial (tissue) pressure, p_im = alpha p_LV.
      Real pim = 0.0;
      /// Arteriolar lumped resistance at the reference viscosity (implicit).
      Real Ra = 0.0;
      /// Venular lumped resistance at the reference viscosity.
      Real Rv = 0.0;
      /// Microvascular compliance, C_tot weighted by the Murray split.
      Real C = 0.0;
      /// Measured outlet area; scales the implicit term R_a Phi_a A (u.n)(v.n).
      Real area = 0.0;
      /// Calibrated branch flow, reference point of the modulation Phi(q).
      Real q0 = 0.0;
      /// Derived resting wall shear rate of the arteriolar limb (1/s).
      Real gammaA = 0.0;
      /// Derived resting wall shear rate of the venular limb (1/s).
      Real gammaV = 0.0;
      /// Diagnostic: apparent-viscosity multiplier of the arteriolar limb.
      Real muA = 1.0;
      /// Diagnostic: apparent-viscosity multiplier of the venular limb.
      Real muV = 1.0;
      /// Diagnostic: stored microvascular volume, V = C p_tm.
      Real vol = 0.0;
  };

  // Carreau-Yasuda blood viscosity, UNIFIED with CoupledLV0DCoronary3D's
  // defaults so the two solvers integrate the same fluid.
  struct CarreauYasuda
  {
      Real mu0 = 0.301;
      Real muInf = 0.0055;
      Real lambda = 16.15;
      Real n = 0.21;
      Real yasuda = 0.77;
      Real gammaRegularization = 1.0e-3;
  };

  /// Rheological model used by the 0D outlet closure.
  enum class RheologyModel
  {
    CarreauYasuda,
    Quemada
  };

  /// Quemada blood viscosity parameters (haematocrit axis); see
  /// CoupledLV0DCoronary3D.h for the rationale and the Cokelet correlations.
  struct Quemada
  {
      Real plasmaViscosity = 0.0017963;
      Real hematocrit = 0.45;
      Real k0 = 3.7;      // derived from phi if <= 0
      Real kInf = 1.66;   // derived from phi if <= 0
      Real gammaC = 2.29; // derived from phi if <= 0
  };

  /// Tabulated WRMS closure and scalar-solve tolerances for the 0D outlet.
  /// mu_ap(tau_w) = tau_w^4 / (4 I(tau_w)) is universal (independent of R, L
  /// and branch), so it is tabulated ONCE and the outlet only interpolates.
  struct OutletFlowLaw
  {
      Real tableTauMin = 1.0e-6;
      Real tableTauMax = 1.0e4;
      int tableNodes = 241;
      int integralSteps = 2000;
      Real shearStepTolerance = 1.0e-12;
      int shearMaxIterations = 200;
      Real outletStepTolerance = 1.0e-9;
      int outletMaxIterations = 50;
      Real zeroFlowTolerance = 1.0e-16;
  };

  /// Universal WRMS apparent-viscosity table, log-log interpolated in the
  /// nominal shear rate gammadot = 4Q/(pi R^3).  Clamped, not extrapolated.
  struct WRMSTable
  {
      std::vector<Real> logGamma;
      std::vector<Real> logMu;

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

  struct Config
  {
  // All outputs go under resultsDir so the P2/P1 and P1/P1 runs can execute
  // concurrently (different cores) without clobbering each other's files.
      std::string resultsDir = "results_p1p1_laplacian";
      std::string xdmfBasename =
        "results_p1p1_laplacian/CoronaryArtery_FSI_Explicit_P1P1";
      std::string csvPath = "results_p1p1_laplacian/CoronaryArtery_FSI_Explicit_P1P1.csv";
      Real meshScale = 1.0e-3;

      Real dt = 1.0e-3;
      size_t nsteps = 4 * static_cast<int>(0.85 / 1.0e-3);

      Real pressureDropScale = 1.0;

  // Fixed epicardial forward pressure drop [Pa] subtracted from EVERY outlet:
  //   p_out_3D = par - pressureDropScale*(par - p_out_RCR) - epicardialDrop.
      Real epicardialDrop = 0.0;

  // Automatic anatomical outlet calibration (CoupledLV0DCoronary3D scheme):
  // outlet areas measured on the mesh, Murray split Q_i ~ r_i^3, and the three
  // lumped quantities that appear in the balance -- R_a, R_v, C -- built from
  // the resting pressure budget dP = par(0) - P_RA:
  //   R_v,i = f_v dP / Q_i,  R_a,i = (1 - f_v) dP / Q_i,  C_i = C_tot w_i.
      bool autoCalibrateOutlets = true;
      Real lcaTargetFlow = 1.5e-6;     // m^3/s  (~90 mL/min; LCA rest ~150-250)

  // Newtonian calibration viscosity mu_N: R_a and R_v are the NEWTONIAN
  // resistances of the budget, so a change of blood properties moves Phi and
  // therefore the flow (normalizing by the running rheology would absorb it).
      Real newtonianCalibrationViscosity = 0.0035;
  // Total microvascular compliance (m^3/Pa), split by the Murray weight.
      Real coronaryComplianceTotal = 4e-10;
  // Venular share of the microvascular pressure budget (head-loss split).
      Real venularPressureFraction = 0.13;
  // Fraction of LV pressure transmitted to the intramyocardial compartment,
  // p_im = intramyocardialFraction * p_LV (bed-averaged, 0.5-0.7).
      Real intramyocardialFraction = 0.7;
  // Right atrial (coronary sinus) drainage pressure.
      Real rightAtrialPressure = 1800.0;

  // Morphometric operating point (r, v) of each limb: intravital-microscopy
  // measurables from which gamma_0 = 4v/r, N, L and T are derived.
      Real arteriolarRadius = 25.0e-6;
      Real arteriolarVelocity = 5.0e-3;
      Real venularRadius = 30.0e-6;
      Real venularVelocity = 3.0e-3;
  // Reference mean transit time for the calibration consistency check (s).
      Real referenceTransitTime = 1.5;

  // Rheological model of the 0D outlet closure.
      RheologyModel rheologyModel = RheologyModel::CarreauYasuda;
      Quemada quemada;

      Real fluidDensity = 1060.0;
      Real pressurePenalty = 0.0;
      Real inletBackflowStabilization = 1.0;
      Real outletBackflowStabilization = 1.0;
      CarreauYasuda viscosity;

      OutletFlowLaw outletFlowLaw;

      Real inletImpedance = 1e3;
      Real inletTangentialDamping = 1e3;

      Real outletResistanceScale = 1.0;

  // Static follower-pressure prestress: ramp 0 -> par(0) in this many
  // increments (0 = disabled)
      size_t prestressSteps = 50;

  // Prestress STATICALLY to prestressFraction*par, then ramp the remaining
  // (1 - prestressFraction) of the loads smoothly from prestressFraction -> 1

      Real prestressFraction = 0.975;
      size_t prestressRampSteps = 10;

  // Projected-VMS convective stabilization scale (lagged Oseen, ported from
  // the working CoupledLV0DCoronary3D fluid).  0 disables the VMS terms and
  // recovers the plain Temam-stabilized fluid; 1 is the nominal scaling.
      Real vmsScale = 1.0;

  // Grad-div (continuity / "div-u") stabilization scale.
      Real gradDivScale = 1.0;

  // ALE mesh-motion (harmonic extension) stiffening.  The lift is
  //   int  w(x) grad d : grad v ,   w(x) = (aleRefSize / h_K)^aleStiffPower ,
  // i.e. Jacobian/element-size stiffening (Stein-Tezduyar-Benney)
      Real aleStiffPower = 0.75;
      Real aleRefSize = 5.0e-4;

  // Kinematic interface velocity for the fluid no-slip / Robin data.  By the
  // geometric conservation law the ALE mesh boundary moves at the BDF1 rate
  //   w|FSI = (d^{n+1} - d^n)/dt = meshVelocity ,
  // so the FLUID at the wall must move with exactly that.
      bool meshConsistentInterfaceVelocity = true;

      bool subtractHeartHandoffOffset = true;

      Real pgpScale = 1.0;

      Real solidDensity = 1060.0;

      Real solidViscosity = 4.e3;

      Real newmarkBeta = 0.25;
      Real newmarkGamma = 0.5;

      size_t couplingIterations = 1; // Picard passes; snes.solve() re-evaluates each pass
      Real couplingTolerance = 1.0e-6; // relative interface-displacement tolerance

      Real robinAlpha = 0.0;
      Real robinGamma = 1.0; // scales the Robin alpha = gamma*sqrt(rho_s E); raise
        // to enforce no-slip harder (smaller interface slip)

      Real heartDisplacementPenalty = 1.e7;
      Real heartDisplacementScale = 1.0;

      // Viscoelastic outlet tethering (spring/dashpot per unit area on the cap
      // annuli).  With the ring band now confined to the INLET, these springs are
      // the ONLY thing holding the outlet ends: 5e1 Pa/m left them essentially
      // free-flying, so the defaults are raised to a perivascular-tethering level
      // (1e4-1e6 Pa/m is the physiological ballpark).
      Real aViscCondition = 1.0e5;
      Real bViscCondition = 1.0e3;

      // Transmural grading of the wall constitutive law, applied as a spatial
      // multiplier m(xi) on the Yeoh energy (see GradedYeoh.h), where xi in [0,1]
      // is the harmonic through-thickness coordinate (0 = lumen/intima side,
      // 1 = adventitia).  Bands: xi < 1/3 intima, middle third media (reference),
      // xi > 2/3 adventitia, blended with smoothsteps of half-width
      // gradeTransitionWidth.
      //
      // Defaults follow the classical healthy-artery picture (Holzapfel-Gasser-
      // Ogden 2000): a thin, mechanically insignificant intima (soft), the
      // load-bearing media as reference, and a stiffer collagenous adventitial
      // sleeve, normalized so the thickness-weighted mean is ~1 and the global
      // compliance of the tuned homogenized wall is preserved.  NOTE: for aged /
      // atherosclerotic coronaries Holzapfel et al. (2005) measured the OPPOSITE
      // low-strain ordering (ground-matrix mu: intima 27.9, media 1.27,
      // adventitia 7.56 kPa); to run that profile set roughly
      // (intima, media, adventitia) = (2.3, 0.10, 0.76).
      Real gradeIntima = 0.5;
      Real gradeMedia = 1.0;
      Real gradeAdventitia = 1.5;
      Real gradeTransitionWidth = 0.08;
  };

  // C1-continuous (smoothstep) LV activation, UNIFIED with the
  // CoupledLV0DCoronary3D waveform: same plateau values, smooth ramps, and a
  // smooth return from the negative plateau (the old piecewise-linear version
  // jumped instantaneously from -20 to 0 at tau = 0.6).
  static Real periodic_activation(Real t)
  {
    const Real T = 0.85;
    const Real tau = t - T * std::floor(t / T);

    auto ss = [](Real s) {
      s = s < 0.0 ? 0.0 : (s > 1.0 ? 1.0 : s);
      return s * s * (3.0 - 2.0 * s);
    };

    const Real tRampStart = 0.15, tRampEnd = 0.21, tPlateauEnd = 0.36;
    const Real tRelaxEnd = 0.45, tNegativeEnd = 0.6;
    const Real positiveValue = 35.0, negativeValue = -20.0;

    if (tau < tRampStart)
      return 0.0;
    if (tau < tRampEnd)
      return positiveValue * ss((tau - tRampStart) / (tRampEnd - tRampStart));
    if (tau < tPlateauEnd)
      return positiveValue;
    if (tau < tRelaxEnd)
      return positiveValue +
        (negativeValue - positiveValue) *
        ss((tau - tPlateauEnd) / (tRelaxEnd - tPlateauEnd));
    if (tau < tNegativeEnd)
      return negativeValue;
    return negativeValue * (1.0 - ss((tau - tNegativeEnd) / (T - tNegativeEnd)));
  }

  static Real load_dependent_relaxation_m0(Real ec)
  {
    const Real lowEc = 0.0;
    const Real highEc = 2.0;
    const Real lowValue = 1.6;
    const Real highValue = 1.0;
    if (ec <= lowEc)
      return lowValue;
    if (ec >= highEc)
      return highValue;
    const Real s = (ec - lowEc) / (highEc - lowEc);
    return (1.0 - s) * lowValue + s * highValue;
  }

  static Real load_dependent_relaxation_dm0(Real ec)
  {
    const Real lowEc = 0.0;
    const Real highEc = 2.0;
    const Real lowValue = 1.6;
    const Real highValue = 1.0;
    if (ec <= lowEc || ec >= highEc)
      return 0.0;
    return (highValue - lowValue) / (highEc - lowEc);
  }

  // C1-continuous (smoothstep) atrial pressure between the same plateau values
  // as before, so the prescribed atrial/venous pressure has no derivative
  // kinks (UNIFIED with CoupledLV0DCoronary3D).
  static Real atrial_pressure(Real t)
  {
    const Real T = 0.85;
    const Real tau = t - T * std::floor(t / T);
    const Real minValue = 500.0;
    const Real maxValue = 1000.0;
    const Real secondThreshold = 1250.0;
    const Real t1 = 0.02;
    const Real t2 = 0.15;
    const Real t3 = 0.17;
    const Real t4 = 0.56;
    const Real t5 = 0.62;
    const Real t6 = 0.85;

    auto ss = [](Real s) {
      s = s < 0.0 ? 0.0 : (s > 1.0 ? 1.0 : s);
      return s * s * (3.0 - 2.0 * s);
    };
    auto ramp = [&](Real a, Real b, Real s) { return a + (b - a) * ss(s); };

    if (tau < t1)
      return ramp(minValue, maxValue, tau / t1);
    if (tau < t2)
      return maxValue;
    if (tau < t3)
      return ramp(maxValue, minValue, (tau - t2) / (t3 - t2));
    if (tau < t4)
      return ramp(minValue, secondThreshold, (tau - t3) / (t4 - t3));
    if (tau < t5)
      return secondThreshold;
    if (tau < t6)
      return ramp(secondThreshold, minValue, (tau - t5) / (t6 - t5));
    return minValue;
  }

  static Model::Input makeModelInput()
  {
    Model::Input in;

    in.rho = 1.0e3;
    in.R0 = 2.36e-2;
    in.d0 = 1.42e-2;

    in.Es = 3.0e7;
    in.mu = 70.0;
    in.eta = 70.0;
    in.alpha = 1.5;
    in.alphaR = 0.12;
    in.k0 = 1.0e5;
    in.sigma0 = 1.5e5;

    in.Rp = 8.0e6;
    in.Cp = 8.e-9;
    in.Rd = 1.0e8;
    in.Cd = 5.0e-10;

    // Carreau-Yasuda set UNIFIED with CoupledLV0DCoronary3D (same fluid in the
    // 0D inlet windkessel, the 3D solver and the 0D outlet closure).
    in.mu_0 = 0.301;
    in.mu_Inf = 0.0055;
    in.lambda = 16.152;
    in.n = 0.21;
    in.m = 0.0035;
    in.yasuda = 0.77;
    in.mu_plasma = 0.0032704;
    in.k_0 = 3.5678;
    in.gamma_c = 10.2754;
    in.k_Inf = 1.5352;
    in.proximalRadius = 0.0125;
    in.proximalLength = 0.4;
    in.distalRadius = 0.00175;
    in.distalLength = 0.2;
    in.windkesselRheology =
      Rodin::Heart::CCMLC2014::Model::WindkesselRheology::CarreauYasuda;

    in.Kat = 2.0e-6;
    in.Kp = 5.0e-10;
    in.Kar = 2.0e-7;
    in.cavityCapacity = 5.0e-12;

    in.localTolerance = 1.0e-12;
    in.localMaxIterations = 50;
    in.localDamping = 1.0;
    in.absRegularization = 1.0e-14;

    in.initFibDef = 0.0;
    in.initActiveStiffness = 0.0;
    in.initActiveStress = 0.0;

    in.pSv = [](Real) { return 1.0e3; };
    in.pAt = atrial_pressure;
    in.u = periodic_activation;
    in.m0 = load_dependent_relaxation_m0;
    in.dm0 = load_dependent_relaxation_dm0;

    using PassiveEnergy = std::decay_t<decltype(in.passiveEnergy)>;
    typename PassiveEnergy::Parameters hp;
    hp.mu1 = 0.0;
    hp.mu2 = 0.0;
    hp.C0 = 1.9e3;
    hp.C1 = 1.1e-1;
    hp.C2 = 1.9e3;
    hp.C3 = 1.1e-1;
    in.passiveEnergy = PassiveEnergy(hp);

    return in;
  }

  static void initializeModel(Model& model, const Model::Input& in)
  {
    model.setMaxIterations(200)
      .setAbsoluteTolerance(1.0e-8)
      .setRelativeTolerance(1.0e-8)
      .setStepTolerance(1.0e-10)
      .setDampingFactor(1.0);

    Model::State s0;
    s0.t = 0.0;
    s0.y = 0.0;
    s0.v = 0.0;
    s0.pv = in.pAt(0.0) - 100.0;
    s0.par = 11000.0;
    s0.pd = 10000.0;
    s0.ec = in.initFibDef;
    s0.gamma = std::sqrt(std::max<Real>(in.initActiveStiffness, 0.0));
    s0.beta = (s0.gamma > 0.0) ? (in.initActiveStress / s0.gamma) : 0.0;
    s0.kc = s0.gamma * s0.gamma;
    s0.tauc = s0.gamma * s0.beta;
    s0.w = in.m0(s0.ec);

    model.initialize(s0);
  }

  // Load the coronary FSI mesh directly into a local (single-process) mesh.
  static MeshType makeMesh(const Config& cfg, const std::string& meshPath)
  {
    MeshType mesh;
    mesh.load(meshPath, IO::FileFormat::MEDIT);

    if (mesh.getSpaceDimension() != 3)
      throw std::runtime_error("Expected a 3D coronary FSI mesh.");

    mesh.scale(cfg.meshScale);

    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);
    mesh.getConnectivity().compute(D, 0);
    mesh.getConnectivity().compute(D, D - 1);
    mesh.getConnectivity().compute(D - 1, D);
    mesh.getConnectivity().compute(D - 1, 0);
    mesh.getConnectivity().compute(D - 1, 1);
    mesh.getConnectivity().compute(1, 0);

    return mesh;
  }

  template <class MeshT>
  static void saveReferenceVertices(
    const MeshT& mesh, std::vector<Math::SpatialPoint>& vertices)
  {
    vertices.resize(mesh.getVertexCount());
    for (auto it = mesh.getVertex(); it; ++it)
      vertices[it->getIndex()] = mesh.getVertexCoordinates(it->getIndex());
  }

  template <class MeshT>
  static void restoreMeshToReference(
    MeshT& mesh, const std::vector<Math::SpatialPoint>& referenceVertices)
  {
    assert(mesh.getVertexCount() == referenceVertices.size());
    for (auto it = mesh.getVertex(); it; ++it)
    {
      const Index vertex = it->getIndex();
      mesh.setVertexCoordinates(vertex, referenceVertices[vertex]);
    }
    mesh.flush();
  }

  template <class MeshT, class FESType, class GridFunctionType>
  static void moveMeshWithVertexDisplacement(MeshT& mesh,
    const std::vector<Math::SpatialPoint>& referenceVertices,
    const FESType& displacementFES, const GridFunctionType& displacement)
  {
    assert(mesh.getVertexCount() == referenceVertices.size());
    const size_t dim = mesh.getSpaceDimension();

    for (auto it = mesh.getVertex(); it; ++it)
    {
      const Index vertex = it->getIndex();
      auto x = referenceVertices[vertex];

      for (Index c = 0; c < static_cast<Index>(dim); ++c)
        x(c) += displacement[displacementFES.getGlobalIndex({0, vertex}, c)];

      mesh.setVertexCoordinates(vertex, x);
    }
    mesh.flush();
  }

  // --------------------------------------------------------------------------
  // FSI interface map between the fluid and solid meshes.
  // --------------------------------------------------------------------------
  struct InterfaceSegment
  {
      Index fluid;
      Index solid;
  };

  struct InterfaceMap
  {
      std::unordered_map<Index, Index> fluidToSolid;
      std::unordered_map<Index, Index> solidToFluid;
      std::vector<InterfaceSegment> segments;
  };

  static Math::SpatialPoint centroid(const MeshType& mesh, const Polytope& polytope)
  {
    Math::SpatialPoint c(mesh.getSpaceDimension());
    c.setZero();

    const auto& vertices = polytope.getVertices();
    for (const auto& v : vertices)
      c += mesh.getVertexCoordinates(v);

    c /= static_cast<Real>(vertices.size());
    return c;
  }

  static std::size_t tagFSIRingBand(MeshType& mesh, Attribute fsi, Attribute ring,
    Attribute inlet, const std::array<Attribute, 6>& outlets)
  {
    const std::size_t faceDim = mesh.getDimension() - 1;

    const auto isCap = [&](const Optional<Attribute>& a) {
      if (!a)
        return false;
      if (*a == inlet)
        return true;
      for (const Attribute o : outlets)
        if (*a == o)
          return true;
      return false;
    };

    // Vertices lying on any cap face; the ring curves are a subset of these.
    std::unordered_set<Index> capVertices;
    for (auto it = mesh.getBoundary(); it; ++it)
    {
      if (!isCap(it->getAttribute()))
        continue;
      for (const auto& v : it->getVertices())
        capVertices.insert(v);
    }

    // Collect (don't mutate during iteration) every FSI face touching a cap vertex.
    std::vector<Index> toRelabel;
    for (auto it = mesh.getBoundary(); it; ++it)
    {
      if (it->getAttribute() != fsi)
        continue;
      for (const auto& v : it->getVertices())
        if (capVertices.count(v))
        {
          toRelabel.push_back(it->getIndex());
          break;
        }
    }

    for (const Index f : toRelabel)
      mesh.setAttribute({faceDim, f}, ring);

    return toRelabel.size();
  }

  // Match each fluid FSI face to its geometric twin on the solid FSI surface by
  // centroid proximity.  The two meshes are assumed CONFORMING on the interface
  // (coincident faces
  static InterfaceMap buildInterfaceMap(
    const MeshType& fluidReferenceMesh, const MeshType& solidReferenceMesh)
  {
    InterfaceMap map;

    // Collect local solid FSI face centroids and the interface bounding box.
    std::vector<std::pair<Math::SpatialPoint, Index>> solidFaces;
    Math::SpatialPoint lo, hi;
    bool haveBox = false;
    for (auto it = solidReferenceMesh.getBoundary(); it; ++it)
    {
      if (it->getAttribute() != BoundarySolid::FSI)
        continue;
      Math::SpatialPoint c = centroid(solidReferenceMesh, *it);
      solidFaces.emplace_back(c, it->getIndex());
      if (!haveBox)
      {
        lo = c;
        hi = c;
        haveBox = true;
      }
      else
      {
        for (Index i = 0; i < static_cast<Index>(c.size()); ++i)
        {
          lo(i) = std::min(lo(i), c(i));
          hi(i) = std::max(hi(i), c(i));
        }
      }
    }

    // Matching tolerance: a small fraction of the interface bounding-box diagonal
    // (large enough to absorb round-off, far smaller than a face).
    const Real diag = haveBox ? Real((hi - lo).norm()) : Real(0);
    const Real tol = std::max(Real(1e-12), Real(1e-6) * diag);

    // Spatial hash of solid centroids on a grid of cell size 'tol'.
    auto cellOf = [&](const Math::SpatialPoint& c) {
      std::array<long long, 3> k{0, 0, 0};
      const Index n = std::min<Index>(3, static_cast<Index>(c.size()));
      for (Index i = 0; i < n; ++i)
        k[i] = static_cast<long long>(std::floor(c(i) / tol));
      return k;
    };
    std::map<std::array<long long, 3>, std::vector<std::size_t>> grid;
    for (std::size_t s = 0; s < solidFaces.size(); ++s)
      grid[cellOf(solidFaces[s].first)].push_back(s);

    // Match each local fluid FSI face to the nearest solid centroid within 'tol'.
    std::size_t fluidCount = 0;
    std::size_t unmatched = 0;
    Real worst = 0.0;
    for (auto it = fluidReferenceMesh.getBoundary(); it; ++it)
    {
      if (it->getAttribute() != BoundaryFluid::FSI)
        continue;
      ++fluidCount;

      const Math::SpatialPoint c = centroid(fluidReferenceMesh, *it);
      const auto base = cellOf(c);

      Real best = std::numeric_limits<Real>::max();
      Index bestSolid = 0;
      bool found = false;
      for (long long dx = -1; dx <= 1; ++dx)
        for (long long dy = -1; dy <= 1; ++dy)
          for (long long dz = -1; dz <= 1; ++dz)
          {
            const std::array<long long, 3> key{base[0] + dx, base[1] + dy, base[2] + dz};
            const auto g = grid.find(key);
            if (g == grid.end())
              continue;
            for (const std::size_t s : g->second)
            {
              const Real d = Real((solidFaces[s].first - c).norm());
              if (d < best)
              {
                best = d;
                bestSolid = solidFaces[s].second;
                found = true;
              }
            }
          }

      if (!found || best > tol)
      {
        ++unmatched;
        if (found)
          worst = std::max(worst, best);
        continue;
      }

      const Index fluidFace = it->getIndex();
      map.fluidToSolid.emplace(fluidFace, bestSolid);
      map.solidToFluid.emplace(bestSolid, fluidFace);
      map.segments.push_back({fluidFace, bestSolid});
    }

    if (unmatched > 0)
    {
      std::ostringstream os;
      os << "buildInterfaceMap: " << unmatched << " of " << fluidCount
         << " fluid FSI face(s) had no solid twin within tol=" << tol
         << " (nearest distance up to " << worst
         << "; solid FSI faces=" << solidFaces.size() << "). "
         << "If solid FSI faces=0 the FSI attribute is wrong/missing in the solid "
            "mesh. If the distance is ~one face size, the meshes are "
            "non-conforming on the interface.";
      throw std::runtime_error(os.str());
    }

    return map;
  }

  [[maybe_unused]] static Point forwardFluidPointToSolid(
    const Point& p, const MeshType& solidReferenceMesh, const InterfaceMap& map)
  {
    const auto found = map.fluidToSolid.find(p.getPolytope().getIndex());

    if (found == map.fluidToSolid.end())
      throw std::runtime_error("Fluid point is not on a mapped FSI face.");

    auto solidFace = solidReferenceMesh.getFace(found->second);

    // Geometric re-localization: matched faces need not share vertex ordering,
    // so the source face's reference coordinate cannot be reused directly.
    const Math::SpatialPoint pc = p.getPhysicalCoordinates();
    Math::SpatialPoint rc;
    solidFace->getTransformation().inverse(rc, pc);

    return Point(*solidFace, rc, pc);
  }

  static Point forwardSolidPointToFluid(
    const Point& p, const MeshType& fluidMesh, const InterfaceMap& map)
  {
    const auto found = map.solidToFluid.find(p.getPolytope().getIndex());

    if (found == map.solidToFluid.end())
      throw std::runtime_error("Solid point is not on a mapped FSI face.");

    auto fluidFace = fluidMesh.getFace(found->second);

    // Geometric re-localization on the twin face -- see forwardFluidPointToSolid.
    const Math::SpatialPoint pc = p.getPhysicalCoordinates();
    Math::SpatialPoint rc;
    fluidFace->getTransformation().inverse(rc, pc);

    return Point(*fluidFace, rc, pc);
  }

  // True iff every component is finite.
  static bool isFiniteVec(const Math::SpatialVector<Real>& x)
  {
    for (Index i = 0; i < x.size(); ++i)
      if (!std::isfinite(x(i)))
        return false;
    return true;
  }

  static void setPETScDefault(const char* key, const char* value)
  {
    PetscBool set = PETSC_FALSE;
    PetscOptionsHasName(PETSC_NULLPTR, PETSC_NULLPTR, key, &set);
    if (!set)
      PetscOptionsSetValue(PETSC_NULLPTR, key, value);
  }

  static void readOptions(Config& cfg)
  {
    char meshPath[PETSC_MAX_PATH_LEN] = {};
    PetscBool meshPathSet = PETSC_FALSE;
    PetscOptionsGetString(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_fsi_mesh", meshPath,
      sizeof(meshPath), &meshPathSet);

    PetscReal dt = cfg.dt;
    PetscBool dtSet = PETSC_FALSE;
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_dt", &dt, &dtSet);
    if (dtSet)
      cfg.dt = dt;

    PetscInt nsteps = static_cast<PetscInt>(cfg.nsteps);
    PetscBool nstepsSet = PETSC_FALSE;
    PetscOptionsGetInt(
      PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_nsteps", &nsteps, &nstepsSet);
    if (nstepsSet)
      cfg.nsteps = static_cast<size_t>(std::max<PetscInt>(0, nsteps));

    PetscInt couplingIterations = static_cast<PetscInt>(cfg.couplingIterations);
    PetscBool couplingIterationsSet = PETSC_FALSE;
    PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_coupling_iterations",
      &couplingIterations, &couplingIterationsSet);
    if (couplingIterationsSet)
      cfg.couplingIterations =
        static_cast<size_t>(std::max<PetscInt>(1, couplingIterations));

    PetscReal couplingTolerance = cfg.couplingTolerance;
    PetscBool couplingToleranceSet = PETSC_FALSE;
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_coupling_tolerance",
      &couplingTolerance, &couplingToleranceSet);
    if (couplingToleranceSet)
      cfg.couplingTolerance = couplingTolerance;

    PetscInt prestressSteps = static_cast<PetscInt>(cfg.prestressSteps);
    PetscBool prestressStepsSet = PETSC_FALSE;
    PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_prestress_steps",
      &prestressSteps, &prestressStepsSet);
    if (prestressStepsSet)
      cfg.prestressSteps = static_cast<size_t>(std::max<PetscInt>(0, prestressSteps));

    PetscReal prestressFraction = cfg.prestressFraction;
    PetscBool prestressFractionSet = PETSC_FALSE;
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_prestress_fraction",
      &prestressFraction, &prestressFractionSet);
    if (prestressFractionSet)
      cfg.prestressFraction = std::min<Real>(1.0, std::max<Real>(0.0, prestressFraction));

    PetscInt prestressRampSteps = static_cast<PetscInt>(cfg.prestressRampSteps);
    PetscBool prestressRampStepsSet = PETSC_FALSE;
    PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_prestress_ramp_steps",
      &prestressRampSteps, &prestressRampStepsSet);
    if (prestressRampStepsSet)
      cfg.prestressRampSteps =
        static_cast<size_t>(std::max<PetscInt>(0, prestressRampSteps));

    PetscReal vmsScale = cfg.vmsScale;
    PetscBool vmsScaleSet = PETSC_FALSE;
    PetscOptionsGetReal(
      PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_vms_scale", &vmsScale, &vmsScaleSet);
    if (vmsScaleSet)
      cfg.vmsScale = std::max<Real>(0.0, vmsScale);

    // Inlet normal impedance Z in  +Z (u.n)(v.n)  on the inlet.  This makes the
    // inlet a Robin BC: p_inlet ~ pin - Z|u.n|, so a large Z artificially LOWERS
    // the inlet pressure below the imposed pin.  Lower it (-> 0) to recover a
    // pure pressure inlet p_inlet = pin = par.
    PetscReal inletImpedance = cfg.inletImpedance;
    PetscBool inletImpedanceSet = PETSC_FALSE;
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_inlet_impedance",
      &inletImpedance, &inletImpedanceSet);
    if (inletImpedanceSet)
      cfg.inletImpedance = std::max<Real>(0.0, inletImpedance);

    PetscReal heartDisplacementPenalty = cfg.heartDisplacementPenalty;
    PetscBool heartDisplacementPenaltySet = PETSC_FALSE;
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_heart_disp_penalty",
      &heartDisplacementPenalty, &heartDisplacementPenaltySet);
    if (heartDisplacementPenaltySet)
      cfg.heartDisplacementPenalty = std::max<Real>(0.0, heartDisplacementPenalty);

    // Scale on the imposed 0D-heart radial displacement (s.y ~ 6 mm is the LV
    // scale, far too large to impose as a normal wall displacement on a coronary).
    PetscReal heartDisplacementScale = cfg.heartDisplacementScale;
    PetscBool heartDisplacementScaleSet = PETSC_FALSE;
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_heart_disp_scale",
      &heartDisplacementScale, &heartDisplacementScaleSet);
    if (heartDisplacementScaleSet)
      cfg.heartDisplacementScale = heartDisplacementScale;

    // Mesh-consistent (BDF1) interface velocity for the fluid no-slip / Robin.
    PetscBool meshConsistentVel =
      cfg.meshConsistentInterfaceVelocity ? PETSC_TRUE : PETSC_FALSE;
    PetscBool meshConsistentVelSet = PETSC_FALSE;
    PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR,
      "-coronary_mesh_consistent_interface_vel", &meshConsistentVel,
      &meshConsistentVelSet);
    if (meshConsistentVelSet)
      cfg.meshConsistentInterfaceVelocity = (meshConsistentVel == PETSC_TRUE);

    // Subtract the 0D-heart displacement present at the dynamic handoff.
    PetscBool subHeartOffset = cfg.subtractHeartHandoffOffset ? PETSC_TRUE : PETSC_FALSE;
    PetscBool subHeartOffsetSet = PETSC_FALSE;
    PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_subtract_heart_offset",
      &subHeartOffset, &subHeartOffsetSet);
    if (subHeartOffsetSet)
      cfg.subtractHeartHandoffOffset = (subHeartOffset == PETSC_TRUE);

    // Newmark gamma (>= 0.5; > 0.5 adds numerical damping to kill structural
    // ringing).
    PetscReal newmarkGamma = cfg.newmarkGamma;
    PetscBool newmarkGammaSet = PETSC_FALSE;
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_newmark_gamma",
      &newmarkGamma, &newmarkGammaSet);
    if (newmarkGammaSet)
    {
      cfg.newmarkGamma = std::max<Real>(0.5, newmarkGamma);
      const Real g = cfg.newmarkGamma + 0.5;
      cfg.newmarkBeta = 0.25 * g * g; // unconditional stability default
    }

    PetscReal newmarkBeta = cfg.newmarkBeta;
    PetscBool newmarkBetaSet = PETSC_FALSE;
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_newmark_beta",
      &newmarkBeta, &newmarkBetaSet);
    if (newmarkBetaSet)
      cfg.newmarkBeta = std::max<Real>(1.0e-6, newmarkBeta);

    PetscReal solidViscosity = cfg.solidViscosity;
    PetscBool solidViscositySet = PETSC_FALSE;
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_solid_viscosity",
      &solidViscosity, &solidViscositySet);
    if (solidViscositySet)
      cfg.solidViscosity = std::max<Real>(0.0, solidViscosity);

    PetscReal gradDivScale = cfg.gradDivScale;
    PetscBool gradDivScaleSet = PETSC_FALSE;
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_graddiv_scale",
      &gradDivScale, &gradDivScaleSet);
    if (gradDivScaleSet)
      cfg.gradDivScale = std::max<Real>(0.0, gradDivScale);

    // ALE mesh-motion stiffening (Jacobian/element-size weighting).
    PetscReal aleStiffPower = cfg.aleStiffPower;
    PetscBool aleStiffPowerSet = PETSC_FALSE;
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_ale_stiff_power",
      &aleStiffPower, &aleStiffPowerSet);
    if (aleStiffPowerSet)
      cfg.aleStiffPower = std::max<Real>(0.0, aleStiffPower);

    PetscReal aleRefSize = cfg.aleRefSize;
    PetscBool aleRefSizeSet = PETSC_FALSE;
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_ale_ref_size",
      &aleRefSize, &aleRefSizeSet);
    if (aleRefSizeSet)
      cfg.aleRefSize = std::max<Real>(1.0e-30, aleRefSize);

    PetscReal outletResistanceScale = cfg.outletResistanceScale;
    PetscBool outletResistanceScaleSet = PETSC_FALSE;
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_outlet_resistance_scale",
      &outletResistanceScale, &outletResistanceScaleSet);
    if (outletResistanceScaleSet)
      cfg.outletResistanceScale = std::max<Real>(0.0, outletResistanceScale);

    PetscReal pressureDropScale = cfg.pressureDropScale;
    PetscBool pressureDropScaleSet = PETSC_FALSE;
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_pressure_drop_scale",
      &pressureDropScale, &pressureDropScaleSet);
    if (pressureDropScaleSet)
      cfg.pressureDropScale = std::max<Real>(0.0, pressureDropScale);

    PetscReal epicardialDrop = cfg.epicardialDrop;
    PetscBool epicardialDropSet = PETSC_FALSE;
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_epicardial_drop",
      &epicardialDrop, &epicardialDropSet);
    if (epicardialDropSet)
      cfg.epicardialDrop = std::max<Real>(0.0, epicardialDrop);

    PetscBool autoCalib = cfg.autoCalibrateOutlets ? PETSC_TRUE : PETSC_FALSE;
    PetscBool autoCalibSet = PETSC_FALSE;
    PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_auto_calibrate",
      &autoCalib, &autoCalibSet);
    if (autoCalibSet)
      cfg.autoCalibrateOutlets = (autoCalib == PETSC_TRUE);

    PetscReal lcaTargetFlow = cfg.lcaTargetFlow;
    PetscBool lcaTargetFlowSet = PETSC_FALSE;
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_lca_flow",
      &lcaTargetFlow, &lcaTargetFlowSet);
    if (lcaTargetFlowSet)
      cfg.lcaTargetFlow = std::max<Real>(1e-12, lcaTargetFlow);

    PetscReal pgpScale = cfg.pgpScale;
    PetscBool pgpScaleSet = PETSC_FALSE;
    PetscOptionsGetReal(
      PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_pgp_scale", &pgpScale, &pgpScaleSet);
    if (pgpScaleSet)
      cfg.pgpScale = std::max<Real>(0.0, pgpScale);

    PetscReal robinAlpha = cfg.robinAlpha;
    PetscBool robinAlphaSet = PETSC_FALSE;
    PetscOptionsGetReal(
      PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_robin_alpha", &robinAlpha, &robinAlphaSet);
    if (robinAlphaSet)
      cfg.robinAlpha = robinAlpha;

    PetscReal robinGamma = cfg.robinGamma;
    PetscBool robinGammaSet = PETSC_FALSE;
    PetscOptionsGetReal(
      PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_robin_gamma", &robinGamma, &robinGammaSet);
    if (robinGammaSet)
      cfg.robinGamma = robinGamma;

    // Transmural grading multipliers (intima / media / adventitia).
    PetscReal gradeVal = 0.0;
    PetscBool gradeSet = PETSC_FALSE;
    PetscOptionsGetReal(
      PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_grade_intima", &gradeVal, &gradeSet);
    if (gradeSet)
      cfg.gradeIntima = std::max<Real>(0.0, gradeVal);
    gradeSet = PETSC_FALSE;
    PetscOptionsGetReal(
      PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_grade_media", &gradeVal, &gradeSet);
    if (gradeSet)
      cfg.gradeMedia = std::max<Real>(0.0, gradeVal);
    gradeSet = PETSC_FALSE;
    PetscOptionsGetReal(
      PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_grade_adventitia", &gradeVal, &gradeSet);
    if (gradeSet)
      cfg.gradeAdventitia = std::max<Real>(0.0, gradeVal);
  }

  // --------------------------------------------------------------------------
  // Universal WRMS apparent-viscosity table (ported from
  // CoupledLV0DCoronary3D::buildWRMSTable).  The
  // Weissenberg-Rabinowitsch-Mooney-Schofield closure gives, for a tube of
  // radius R and length L carrying a generalized-Newtonian fluid,
  //
  //   Q = pi R^3 I(tau_w) / tau_w^3,   I(tau_w) = int_0^{tau_w} tau^2 gd dtau,
  //
  // and writing the same flow as Hagen-Poiseuille with an apparent viscosity
  // the radius and the length cancel identically:
  //
  //   mu_ap(tau_w) = tau_w^4 / (4 I(tau_w)),  gd_nom = 4Q/(pi R^3) = tau_w/mu_ap.
  //
  // mu_ap is therefore a UNIVERSAL function of the wall shear stress for a
  // given rheology: build it once, share it across every outlet and both
  // limbs.  Newtonian check: I = tau_w^4/(4 mu) gives mu_ap = mu.
  // --------------------------------------------------------------------------
  static WRMSTable buildWRMSTable(const Config& cfg, const CarreauYasuda& visc)
  {
    const auto& law = cfg.outletFlowLaw;

    const Real mu0 = visc.mu0;
    const Real muInf = visc.muInf;
    const Real lambda = visc.lambda;
    const Real n = visc.n;
    const Real yasuda = visc.yasuda;
    const Real delta = mu0 - muInf;

    // Constitutive law.  Quemada separates what Carreau-Yasuda entangles: the
    // haematocrit sets the high-shear level and k_0 the low-shear aggregation
    // rise.  k_0, k_inf, gamma_c follow the Cokelet correlations when left
    // non-positive (the law has a packing limit phi_max = 2/k).
    const bool quemada = (cfg.rheologyModel == RheologyModel::Quemada);
    const auto& qp = cfg.quemada;

    const Real phi = std::clamp<Real>(qp.hematocrit, 0.0, 0.75);
    const Real p2 = phi * phi;
    const Real p3 = p2 * phi;

    const Real qk0 =
      (qp.k0 > 0.0) ? qp.k0 : std::exp(3.874 - 10.41 * phi + 13.8 * p2 - 6.738 * p3);
    const Real qkInf = (qp.kInf > 0.0)
      ? qp.kInf
      : std::exp(1.3435 - 2.803 * phi + 2.711 * p2 - 0.6479 * p3);
    const Real qgc = (qp.gammaC > 0.0)
      ? qp.gammaC
      : std::exp(-6.1508 + 27.923 * phi - 25.6 * p2 + 3.697 * p3);

    auto muQ = [&](Real g) -> Real {
      const Real s = std::sqrt(std::max<Real>(g, 0.0) / qgc);
      const Real k = (qk0 + qkInf * s) / (1.0 + s);
      const Real b = std::max<Real>(1.0 - 0.5 * k * phi, 1e-3);
      return qp.plasmaViscosity / (b * b);
    };

    auto mu = [&](Real g) -> Real {
      if (quemada)
        return muQ(g);
      return muInf +
        delta * std::pow(1.0 + std::pow(lambda * g, yasuda), (n - 1.0) / yasuda);
    };

    auto dmu = [&](Real g) -> Real {
      if (quemada)
      {
        const Real h = 1e-6 * std::max<Real>(g, 1e-12);
        return (muQ(g + h) - muQ(std::max<Real>(g - h, 0.0))) / (2.0 * h);
      }

      const Real base = 1.0 + std::pow(lambda * g, yasuda);

      return delta * (n - 1.0) * std::pow(base, (n - 1.0 - yasuda) / yasuda) *
        std::pow(lambda, yasuda) * std::pow(g, yasuda - 1.0);
    };

    // gd(tau_w): tau = mu(gd) gd is strictly increasing, so a bisection-
    // safeguarded Newton on log gd converges globally.  Run once per node.
    auto shearAt = [&](Real tauW) -> Real {
      Real lo = std::log(1e-14);
      Real hi = std::log(std::max<Real>(tauW / muInf, 1e-12)) + 1.0;
      Real s = 0.5 * (lo + hi);

      for (int it = 0; it < law.shearMaxIterations; ++it)
      {
        const Real g = std::exp(s);
        const Real f = mu(g) * g - tauW;

        if (f < 0.0)
          lo = s;
        else
          hi = s;

        const Real df = (mu(g) + g * dmu(g)) * g; // d f / d log g
        Real sNext = (std::abs(df) > 0.0) ? s - f / df : 0.5 * (lo + hi);

        if (!(sNext > lo && sNext < hi))
          sNext = 0.5 * (lo + hi);

        const Real step = std::abs(sNext - s);
        s = sNext;

        if (step < law.shearStepTolerance)
          break;
      }

      return std::exp(s);
    };

    // I = int_0^{gd_w} gd^3 mu^2 (mu + gd mu') dgd (same integral after the
    // change of variable tau = mu gd).  Composite Simpson.
    auto rheologicalIntegral = [&](Real gammaW) -> Real {
      const int m = 2 * (law.integralSteps / 2); // even
      const Real h = gammaW / static_cast<Real>(m);

      auto f = [&](Real g) -> Real {
        if (g <= 0.0)
          return 0.0;
        const Real mg = mu(g);
        return g * g * g * mg * mg * (mg + g * dmu(g));
      };

      Real sum = f(0.0) + f(gammaW);

      for (int i = 1; i < m; ++i)
        sum += ((i % 2 == 1) ? 4.0 : 2.0) * f(static_cast<Real>(i) * h);

      return sum * h / 3.0;
    };

    WRMSTable table;

    const int nodes = std::max(2, law.tableNodes);
    const Real logTauMin = std::log(law.tableTauMin);
    const Real logTauMax = std::log(law.tableTauMax);

    table.logGamma.reserve(static_cast<std::size_t>(nodes));
    table.logMu.reserve(static_cast<std::size_t>(nodes));

    for (int i = 0; i < nodes; ++i)
    {
      const Real tauW = std::exp(logTauMin +
        (logTauMax - logTauMin) * static_cast<Real>(i) / static_cast<Real>(nodes - 1));

      const Real gammaW = shearAt(tauW);
      const Real I = rheologicalIntegral(gammaW);

      if (!(I > 0.0) || !std::isfinite(I))
        continue;

      const Real muAp = std::pow(tauW, 4.0) / (4.0 * I);

      if (!(muAp > 0.0) || !std::isfinite(muAp))
        continue;

      const Real gammaNom = tauW / muAp;

      // Nodes must stay strictly increasing for the lookup to be well posed.
      if (!table.logGamma.empty() && std::log(gammaNom) <= table.logGamma.back())
        continue;

      table.logGamma.push_back(std::log(gammaNom));
      table.logMu.push_back(std::log(muAp));
    }

    if (table.logGamma.size() < 2)
    {
      // Degenerate rheology: fall back to a constant apparent viscosity.
      table.logGamma = {std::log(1e-6), std::log(1e6)};
      table.logMu = {std::log(mu0), std::log(mu0)};
    }

    return table;
  }

  // --------------------------------------------------------------------------
  // R-mu coronary outlet update (ported from
  // CoupledLV0DCoronary3D::updateOutlet0D):
  //
  //   3D --[R_a Phi_a, assembled implicitly]--> (p_c, C) --[R_v Phi_v]--> P_RA
  //                                                ^
  //                                             p_im(t)
  //
  // One state, p_tm = p_c - p_im:
  //   C d(p_tm)/dt = q_a - q_v,           q_a = Q  (the measured 3D flux)
  //   q_v = [p_c - max(p_im, P_RA)] / (R_v Phi_v)  if the lumen is open
  //   q_v = 0                                       if the throat is shut
  //   p_out = p_im + p_tm       (+ R_a Phi_a Q, assembled in the 3D form).
  //
  // The gate sits on the STATE of the lumen, not on the sign of the driving
  // pressure: with p_im <= P_RA nothing can collapse and q_v takes either
  // sign; with p_im > P_RA the throat is a check valve, shut once p_tm <= 0.
  // Implicit Euler + scalar Newton on p_tm; dq_v/dp_tm >= 0 everywhere, so
  // R' = C/dt + dq_v/dp_tm > 0 and the iteration converges globally.
  // --------------------------------------------------------------------------
  static void updateOutlet0D(const Config& cfg, const WRMSTable& wrms, const Model& model,
    RCR& bc, Real Q, Real dt)
  {
    const auto& s = model.getState();
    const auto& law = cfg.outletFlowLaw;

    const Real pim = cfg.intramyocardialFraction * s.pv;
    const Real ptmOld = bc.ptm;

    // Drainage pressure of the Starling throat: right atrium outside the
    // myocardium, the collapse pressure p_im inside.
    const Real pDrain = std::max<Real>(pim, cfg.rightAtrialPressure);

    const bool waterfall = pim > cfg.rightAtrialPressure;

    auto venousLumenOpen = [&](Real p) { return !waterfall || p > 0.0; };

    const Real Rv = std::max<Real>(bc.Rv, 1e-300);
    const Real C = std::max<Real>(bc.C, 1e-300);
    const Real q0 = std::max<Real>(std::abs(bc.q0), 1e-300);

    // Rheological modulation: the universal WRMS curve at the nominal shear
    // rate of the limb, normalized by the FIXED Newtonian calibration
    // viscosity mu_N (never the running rheology).
    const Real muN = std::max<Real>(cfg.newtonianCalibrationViscosity, 1e-300);

    auto viscosityFactorV = [&](Real q) -> Real {
      const Real aq = std::abs(q);
      const Real g =
        (aq < law.zeroFlowTolerance) ? law.zeroFlowTolerance : bc.gammaV * aq / q0;
      return wrms(g) / muN;
    };

    // Phi_v depends on |q_v|, so it is held fixed inside each Newton step and
    // refreshed between steps, seeded from the previous drainage (at zero flow
    // the CY plateau gives mu_0, ~ 80x mu_inf).
    Real ptm = ptmOld;
    Real qv = bc.qd;
    Real phiV = viscosityFactorV(qv);
    bool converged = false;

    for (int it = 0; it < law.outletMaxIterations; ++it)
    {
      phiV = viscosityFactorV(qv);

      const Real drive = (ptm + pim) - pDrain;
      const Real Gv = 1.0 / (Rv * phiV);
      const bool open = venousLumenOpen(ptm);

      qv = open ? drive * Gv : 0.0;

      const Real R = C * (ptm - ptmOld) / dt - Q + qv;
      const Real J = C / dt + (open ? Gv : 0.0);

      const Real d = -R / J;
      ptm += d;

      if (std::abs(d) < law.outletStepTolerance * (1.0 + std::abs(ptm)))
      {
        converged = true;
        break;
      }
    }

    if (!converged || !std::isfinite(ptm))
    {
      std::cerr << "Warning: coronary outlet solve did not converge. "
                << "Keeping previous state.\n";
      ptm = ptmOld;
    }

    // Final consistent evaluation at the converged state.
    const Real pc = ptm + pim;
    const Real drive = pc - pDrain;

    phiV = viscosityFactorV(qv);
    qv = venousLumenOpen(ptm) ? drive / (Rv * phiV) : 0.0;

    bc.ptm = ptm;
    bc.pim = pim;
    bc.pc = pc;
    bc.qd = qv;
    bc.vol = C * std::max<Real>(ptm, 0.0);
    bc.muV = phiV;

    // Arteriolar modulation, exported here and consumed by the next 3D solve
    // through the implicit boundary term R_a Phi_a A (u.n)(v.n).
    const Real aQ = std::abs(Q);
    const Real gammaA =
      (aQ < law.zeroFlowTolerance) ? law.zeroFlowTolerance : bc.gammaA * aQ / q0;
    bc.muA = wrms(gammaA) / muN;

    // Pressure applied to the 3D outlet as a Neumann traction.  The resistive
    // part R_a Phi_a Q is NOT included here: assembling it implicitly is what
    // removes the one-step lag on the dominant resistance.
    bc.pout = pc;
  }

} // namespace

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  // Robust direct sub-solvers by default (overridable on the command line).
  setPETScDefault("-ksp_type", "preonly");
  setPETScDefault("-pc_type", "lu");
  setPETScDefault("-pc_factor_mat_solver_type", "mumps");

  // Mass-matrix projections (VMS subscales, gradient recovery for the WSS)
  // are SPD: CG + Jacobi, never the global direct solver (Atrium pattern).
  setPETScDefault("-coronary_mass_ksp_type", "cg");
  setPETScDefault("-coronary_mass_pc_type", "jacobi");

  // Sequential build: no boost::mpi environment / communicator / context.
  // PETSc runs in serial (PETSC_COMM_SELF) and isRoot is trivially true.
  const bool isRoot = true;

  try
  {
    Config cfg;
    readOptions(cfg);

    // Create the per-case results directory (idempotent) so XDMF/CSV writes
    // land in their own folder.
    if (!cfg.resultsDir.empty())
    {
      std::error_code ec;
      std::filesystem::create_directories(cfg.resultsDir, ec);
      if (ec && isRoot)
        std::cerr << "Warning: could not create results directory '" << cfg.resultsDir
                  << "': " << ec.message() << '\n';
    }

    Model::Input modelInput = makeModelInput();
    Model model(modelInput);
    initializeModel(model, modelInput);

    const std::string fluidMesh =
      "../resources/examples/Heart/coronaria_estenosis_80.mesh";
    MeshType meshFluid = makeMesh(cfg, fluidMesh);
    const size_t dimFluid = meshFluid.getSpaceDimension();

    const std::string solidMesh =
      "../resources/examples/Heart/coronaria_80_prismatica.mesh";
    MeshType meshSolid = makeMesh(cfg, solidMesh);
    const size_t dimSolid = meshSolid.getSpaceDimension();

    // Spatial dimension shared by both blocks (used by the coupling fields).
    const size_t dim = dimFluid;

    // Ring band ONLY around the INLET.  The outlets are held by the
    // viscoelastic boundary condition (a, b below): ringing them too would
    // clamp the one-element wall band next to every outlet cap and pin the
    // outlet ends rigidly no matter how soft the springs are (with a single
    // wedge layer every cap node is one edge away from the ring).  The
    // "outlets" argument is therefore filled with the inlet attribute so only
    // inlet-adjacent FSI faces are relabeled to 99.
    const std::array<Attribute, 6> fluidInletOnly{
      {BoundaryFluid::Inlet, BoundaryFluid::Inlet, BoundaryFluid::Inlet,
        BoundaryFluid::Inlet, BoundaryFluid::Inlet, BoundaryFluid::Inlet}};
    const std::array<Attribute, 6> solidInletOnly{
      {BoundarySolid::Inlet, BoundarySolid::Inlet, BoundarySolid::Inlet,
        BoundarySolid::Inlet, BoundarySolid::Inlet, BoundarySolid::Inlet}};
    const std::size_t fluidRingFaces = tagFSIRingBand(meshFluid, BoundaryFluid::FSI,
      BoundaryFluid::FSIRing, BoundaryFluid::Inlet, fluidInletOnly);
    const std::size_t solidRingFaces = tagFSIRingBand(meshSolid, BoundarySolid::FSI,
      BoundarySolid::FSIRing, BoundarySolid::Inlet, solidInletOnly);
    if (isRoot)
      std::cout << "FSI ring band: fluid=" << fluidRingFaces
                << " face(s), solid=" << solidRingFaces << " face(s)\n";

    std::vector<Math::SpatialPoint> referenceVertices;
    saveReferenceVertices(meshFluid, referenceVertices);

    const InterfaceMap interfaceMap = buildInterfaceMap(meshFluid, meshSolid);

    // Fluid problem
    using VelocityFES = H1<1, Math::SpatialVector<Real>, MeshType>;
    using PressureFES = H1<1, Real, MeshType>;
    using DisplacementFluidFES = H1<1, Math::SpatialVector<Real>, MeshType>;

    // Solid problem
    using DisplacementFES = H1<1, Math::SpatialVector<Real>, MeshType>;

    // Laplacian
    using LaplacianFES = H1<1, Real, MeshType>;

    VelocityFES uh(std::integral_constant<size_t, 1>{}, meshFluid, dimFluid);
    PressureFES ph(std::integral_constant<size_t, 1>{}, meshFluid);
    DisplacementFluidFES dfh(std::integral_constant<size_t, 1>{}, meshFluid, dimFluid);

    DisplacementFES dh(std::integral_constant<size_t, 1>{}, meshSolid, dimSolid);

    LaplacianFES lh(std::integral_constant<size_t, 1>{}, meshSolid);

    // Fluid trial/test (velocity-pressure).
    PETSc::Variational::TrialFunction u(uh);
    PETSc::Variational::TrialFunction p(ph);
    PETSc::Variational::TestFunction v(uh);
    PETSc::Variational::TestFunction q(ph);

    // Solid (displacement increment) trial/test.
    PETSc::Variational::TrialFunction d(dh);
    PETSc::Variational::TestFunction w(dh);

    u.setName("u");
    p.setName("p");
    d.setName("disp");

    // Fluid state (uOld/pOld hold the most recent fluid solution).
    PETSc::Variational::GridFunction uOld(uh);
    PETSc::Variational::GridFunction pOld(ph);
    PETSc::Variational::GridFunction one(ph);
    PETSc::Variational::GridFunction shearWall(uh);
    PETSc::Variational::GridFunction gradRec0(uh);
    PETSc::Variational::GridFunction gradRec1(uh);
    PETSc::Variational::GridFunction gradRec2(uh);

    // ALE problem
    PETSc::Variational::GridFunction aleDisp(uh); // full mesh displacement (this step)
    PETSc::Variational::GridFunction aleDispOld(
      uh); // full mesh displacement (previous step)
    PETSc::Variational::GridFunction meshVelocity(uh); // (aleDisp - aleDispOld) / dt

    // Solid / mesh-motion state.
    PETSc::Variational::GridFunction dState(dh); // total solid displacement (current)
    PETSc::Variational::GridFunction dOld(dh); // total solid displacement (previous step)
    PETSc::Variational::GridFunction dIter(dh); // coupling iterate (= dState)
    PETSc::Variational::GridFunction etaState(
      dh); // solid displacement increment (SNES state)
    PETSc::Variational::GridFunction dPred(dh); // Newmark displacement predictor
    PETSc::Variational::GridFunction vPred(dh); // Newmark velocity predictor
    PETSc::Variational::GridFunction solidVelocity(dh);
    PETSc::Variational::GridFunction solidAcceleration(dh);
    PETSc::Variational::GridFunction solidVelocityOld(dh);
    PETSc::Variational::GridFunction solidAccelerationOld(dh);

    PETSc::Variational::GridFunction fluidTraction(dfh);
    PETSc::Variational::GridFunction tractionTransfer(dfh);
    PETSc::Variational::GridFunction uWall(dfh);

    // Laplacian  trial/test.
    PETSc::Variational::TrialFunction l(lh);
    PETSc::Variational::TestFunction t(lh);

    // Through-thickness (transmural) coordinate xi: harmonic between the
    // luminal (FSI) surface, xi = 0, and the adventitial surface (Outer),
    // xi = 1.  Solved once at startup (below, next to the heart-weight
    // laplacian); works for ANY number of wedge layers.  Initialized to 0.5
    // (media) so a premature evaluation before the solve grades to 1.
    PETSc::Variational::TrialFunction xiT(lh);
    PETSc::Variational::TestFunction xiTest(lh);
    PETSc::Variational::GridFunction xi(lh);
    xi = 0.5;

    auto zero = VectorFunction(dim, [&](const Point&) {
      Math::SpatialVector<Real> value(dim);
      value.setZero();
      return value;
    });

    uOld = zero;
    shearWall = zero;
    gradRec0 = zero;
    gradRec1 = zero;
    gradRec2 = zero;
    pOld = 0.0;
    one = 1.0;
    dState = zero;
    dOld = zero;
    dIter = zero;
    etaState = zero;
    dPred = zero;
    vPred = zero;
    solidVelocity = zero;
    solidAcceleration = zero;
    solidVelocityOld = zero;
    solidAccelerationOld = zero;
    aleDisp = zero;
    aleDispOld = zero;
    meshVelocity = zero;
    fluidTraction = zero;
    tractionTransfer = zero;
    uWall = zero;

    uOld.setName("FluidVelocity");
    shearWall.setName("shearStress");
    pOld.setName("FluidPressure");
    dState.setName("Displacement");
    solidVelocity.setName("SolidVelocity");
    meshVelocity.setName("ALEMeshVelocity");
    fluidTraction.setName("FluidTraction");
    aleDisp.setName("ALEDisp");

    IO::XDMF xdmf_fluid(cfg.xdmfBasename + "fluid");
    xdmf_fluid.setMesh(meshFluid);
    xdmf_fluid.add("FluidVelocity", uOld);
    xdmf_fluid.add("FluidPressure", pOld);
    xdmf_fluid.add("ALEMeshVelocity", meshVelocity);
    xdmf_fluid.add("FluidTraction", fluidTraction);
    xdmf_fluid.add("ALEDisp", aleDisp);
    xdmf_fluid.add("shearStress", shearWall);
    xdmf_fluid.write(0.0).flush();

    IO::XDMF xdmf_solid(cfg.xdmfBasename + "solid");
    xdmf_solid.setMesh(meshSolid);
    xdmf_solid.add("Displacement", dState);
    xdmf_solid.add("SolidVelocity", solidVelocity);
    xdmf_solid.write(0.0).flush();

    IO::XDMF xdmf_laplacian(cfg.xdmfBasename + "laplacian");
    xdmf_laplacian.setMesh(meshSolid);
    xdmf_laplacian.add("laplacian", l.getSolution());

    std::ofstream csv;
    if (isRoot)
    {
      csv.open(cfg.csvPath);
      csv << "t,lv_y,lv_v,lv_pv,lv_par,lv_pd,q_in,q_out_total";
      for (const Attribute outlet : BoundaryFluid::Outlets)
        csv << ",q_out_" << outlet << ",p_out_" << outlet;
      // R-mu mechanism diagnostics (averaged/summed over the outlets): p_im,
      // mean p_tm, the two viscosity ratios and the stored volume V = C p_tm.
      csv << ",pim,ptm_mean,phi_a_mean,phi_v_mean,stored_volume";
      csv << ",e_interface,coupling_rel\n";
    }

    std::map<Attribute, RCR> wk;
    for (const Attribute outlet : BoundaryFluid::Outlets)
      wk.emplace(outlet, RCR{});

    // Universal WRMS apparent-viscosity table: built once, shared by every
    // outlet and both limbs.
    const WRMSTable wrms = buildWRMSTable(cfg, cfg.viscosity);

    // ---- Anatomical outlet areas (needed by the implicit R_a Phi_a A term
    // regardless of the calibration path) --------------------------------
    const Real PI = std::numbers::pi_v<Real>;
    std::map<Attribute, Real> rEq;
    Real sumR3 = 0.0;
    {
      PETSc::Variational::TestFunction qCal(ph);
      LinearForm<PressureFES, ::Vec> calForm(qCal);
      for (const Attribute tag : BoundaryFluid::Outlets)
      {
        calForm = BoundaryIntegral(one, qCal).over(tag);
        calForm.assemble();
        const Real area = std::max<Real>(calForm(one), 1e-12);
        wk.at(tag).area = area;
        rEq[tag] = std::sqrt(area / PI);
        sumR3 += rEq[tag] * rEq[tag] * rEq[tag];
      }
    }

    // ---- Anatomical outlet calibration (CoupledLV0DCoronary3D scheme) ---
    //
    // Three divisions per outlet: exactly the three lumped quantities that
    // appear in the balance (R_a, R_v, C).  With dP = par(0) - P_RA the
    // resting budget and the Murray split w_i = r_i^3 / sum r_j^3:
    //   Q_i = Q_tot w_i,  C_i = C_tot w_i,
    //   R_v,i = f_v dP / Q_i,  R_a,i = (1 - f_v) dP / Q_i,
    //   p_tm,i(0) = P_RA + f_v dP Phi_v0 - alpha p_LV(0).
    if (cfg.autoCalibrateOutlets)
    {
      const Real par0 = model.getState().par;
      const Real dP = std::max<Real>(par0 - cfg.rightAtrialPressure, 1.0);

      const Real fv = cfg.venularPressureFraction;
      const Real dPv = std::max<Real>(fv * dP, 1.0);
      const Real dPa = std::max<Real>((1.0 - fv) * dP, 1.0);

      const Real pimRest = cfg.intramyocardialFraction * model.getState().pv;
      const Real muN = std::max<Real>(cfg.newtonianCalibrationViscosity, 1e-300);

      // Morphometric operating point (r, v) per limb; L and T are OUTPUTS:
      //   g_0 = 4 v / r,  N = Q_i/(pi r^2 v),  L = r dP_share/(2 mu_N g_0),
      //   T = L / v.
      const Real ra = std::max<Real>(cfg.arteriolarRadius, 1e-12);
      const Real va = std::max<Real>(cfg.arteriolarVelocity, 1e-12);
      const Real rv = std::max<Real>(cfg.venularRadius, 1e-12);
      const Real vv = std::max<Real>(cfg.venularVelocity, 1e-12);

      const Real gammaA0 = 4.0 * va / ra;
      const Real gammaV0 = 4.0 * vv / rv;

      const Real La = ra * dPa / (2.0 * muN * gammaA0);
      const Real Lv = rv * dPv / (2.0 * muN * gammaV0);

      const Real Ta = La / va;
      const Real Tv = Lv / vv;

      const Real phiA0 = wrms(gammaA0) / muN;
      const Real phiV0 = wrms(gammaV0) / muN;

      // Initial condition: steady state of the ACTUAL (non-Newtonian) network,
      // so the run does not open with a spurious transient.
      const Real ptmRest = cfg.rightAtrialPressure + dPv * phiV0 - pimRest;

      for (const Attribute tag : BoundaryFluid::Outlets)
      {
        const Real w = (rEq[tag] * rEq[tag] * rEq[tag]) / std::max(sumR3, 1e-30);
        const Real Qi = std::max<Real>(cfg.lcaTargetFlow * w, 1e-12);

        auto& bc = wk.at(tag);

        bc.q0 = Qi;
        bc.Ra = dPa / Qi;
        bc.Rv = dPv / Qi;
        bc.C = std::max<Real>(cfg.coronaryComplianceTotal * w, 1e-300);

        bc.gammaA = gammaA0;
        bc.gammaV = gammaV0;

        // Steady state of the calibrated network as initial condition.
        bc.ptm = ptmRest;
        bc.pim = pimRest;
        bc.pc = ptmRest + pimRest;
        bc.pout = bc.pc;
        bc.qd = Qi;
        bc.vol = bc.C * std::max<Real>(ptmRest, 0.0);
        bc.muA = phiA0;
        bc.muV = phiV0;

        if (isRoot)
          Alert::Info() << "  [calib] outlet " << tag << "  A=" << bc.area << " m^2"
                        << "  Q=" << (Qi * 6.0e7) << " mL/min"
                        << "  Ra=" << bc.Ra << "  Rv=" << bc.Rv << " Pa s/m^3"
                        << "  C=" << bc.C << " m^3/Pa"
                        << "  tau=C*Rv=" << (bc.C * bc.Rv) << " s"
                        << "  ptm0=" << ptmRest << "  pc0/pout0=" << bc.pc << "/"
                        << bc.pout << " Pa (par=" << par0 << ", pim=" << pimRest << ")"
                        << Alert::Raise;
      }

      if (isRoot)
      {
        Alert::Info() << "  [calib] WRMS nodes=" << wrms.logGamma.size() << "  reologia="
                      << (cfg.rheologyModel == RheologyModel::Quemada ? "Quemada"
                                                                      : "Carreau-Yasuda")
                      << "  mu_ap(" << gammaA0 << ")=" << wrms(gammaA0) << "  mu_ap("
                      << gammaV0 << ")=" << wrms(gammaV0) << "  (mu_N=" << muN << ")"
                      << "  Phi_a0=" << phiA0 << "  Phi_v0=" << phiV0 << "  |  T_a=" << Ta
                      << " s  T_v=" << Tv << " s  (referencia "
                      << cfg.referenceTransitTime << " s)" << Alert::Raise;

        const Real Tref = cfg.referenceTransitTime;
        if (Ta < 0.5 * Tref || Ta > 2.0 * Tref || Tv < 0.5 * Tref || Tv > 2.0 * Tref)
          Alert::Warning()
            << "  [calib] el tiempo de transito derivado se aparta mas de 2x "
            << "de la referencia: calibre, velocidad y reparto de presion no "
            << "son consistentes con un solo tramo efectivo por rama." << Alert::Raise;
      }
    }
    else
    {
      // Diagnostic path: neutral, non-degenerate constants so the outlet
      // still integrates without the calibration.
      for (const Attribute tag : BoundaryFluid::Outlets)
      {
        auto& bc = wk.at(tag);
        const Real Qi =
          cfg.lcaTargetFlow / static_cast<Real>(BoundaryFluid::Outlets.size());
        bc.q0 = Qi;
        bc.Ra = 4.5e9;
        bc.Rv = 6.8e8;
        bc.C =
          cfg.coronaryComplianceTotal / static_cast<Real>(BoundaryFluid::Outlets.size());
        bc.ptm = 1400.0;
        bc.pc = bc.ptm;
        bc.pout = bc.pc;
        bc.qd = Qi;
        bc.gammaA = 4.0 * cfg.arteriolarVelocity / cfg.arteriolarRadius;
        bc.gammaV = 4.0 * cfg.venularVelocity / cfg.venularRadius;
        bc.vol = bc.C * std::max<Real>(bc.ptm, 0.0);
      }
    }

    Real pinValue = model.getState().par;
    std::map<Attribute, Real> outletPressureValue;
    for (const Attribute outlet : BoundaryFluid::Outlets)
      outletPressureValue[outlet] = wk.at(outlet).pout;

    auto pin = RealFunction([&](const Point&) { return pinValue; });
    auto pout0 = RealFunction(
      [&](const Point&) { return outletPressureValue[BoundaryFluid::Outlets[0]]; });
    auto pout1 = RealFunction(
      [&](const Point&) { return outletPressureValue[BoundaryFluid::Outlets[1]]; });
    auto pout2 = RealFunction(
      [&](const Point&) { return outletPressureValue[BoundaryFluid::Outlets[2]]; });
    auto pout3 = RealFunction(
      [&](const Point&) { return outletPressureValue[BoundaryFluid::Outlets[3]]; });
    auto pout4 = RealFunction(
      [&](const Point&) { return outletPressureValue[BoundaryFluid::Outlets[4]]; });
    auto pout5 = RealFunction(
      [&](const Point&) { return outletPressureValue[BoundaryFluid::Outlets[5]]; });

    // ----------------------------------------------------------------------
    // Constitutive and time-integration constants.
    // ----------------------------------------------------------------------
    const Real dt = cfg.dt;
    const Real betaN = cfg.newmarkBeta;
    const Real gammaN = cfg.newmarkGamma;
    const Real solidMass = cfg.solidDensity / (betaN * dt * dt);
    const Real solidVelocityCoeff = gammaN / (betaN * dt);

    const Real yeohC1 = 80000.0;
    const Real yeohC2 = 400000.0;
    const Real yeohC3 = 5000000.0;
    const Real yeohKappa = 12000000.0;

    const Real solidShearEquiv = 2.0 * yeohC1;
    const Real solidYoungEquiv =
      9.0 * yeohKappa * solidShearEquiv / (3.0 * yeohKappa + solidShearEquiv);

    // alpha = gamma * sqrt(rho_s * E_eq) unless explicitly overridden (>0).
    const Real robinAlpha = (cfg.robinAlpha > 0.0)
      ? cfg.robinAlpha
      : cfg.robinGamma * std::sqrt(cfg.solidDensity * solidYoungEquiv);
    if (isRoot)
      Alert::Info() << "Robin parameter alpha = " << robinAlpha
                    << "  (gamma * sqrt(rho_s E_eq), gamma = " << cfg.robinGamma
                    << ", E_eq = " << solidYoungEquiv << " Pa)" << Alert::Raise;
    const Real robinVelocityCoeff = robinAlpha * solidVelocityCoeff;

    const auto& cy = cfg.viscosity;
    const Real gammaReg = cy.gammaRegularization;
    const Real deltaMu = cy.mu0 - cy.muInf;

    const auto normalFluid = BoundaryNormal(meshFluid);

    // ------------------------------------------------------------------
    // Transmural (intima/media/adventitia) grading of the wall law.  The
    // multiplier m(xi) blends the three band values with smoothsteps at
    // xi = 1/3 and 2/3 (half-width gradeTransitionWidth); the Yeoh energy is
    // jointly linear in (c1, c2, c3, kappa), so GradedYeoh scales stress and
    // tangent by exactly the same factor -- the Newton tangent stays
    // consistent.  xi is evaluated through the P1 GridFunction solved below,
    // so the profile is sampled AT EACH QUADRATURE POINT.
    //
    // NOTE: with a single wedge layer through the thickness the profile is
    // under-resolved (each element spans all three bands); for a faithful
    // three-band wall rebuild the solid with >= 3 wedge layers -- this code
    // needs no change, xi is harmonic and mesh-agnostic.
    // ------------------------------------------------------------------
    auto wallGrade = [&](const Point& p) -> Real {
      const Real s = std::clamp<Real>(xi.getValue(p), 0.0, 1.0);
      auto ss = [](Real u) {
        u = u < 0.0 ? 0.0 : (u > 1.0 ? 1.0 : u);
        return u * u * (3.0 - 2.0 * u);
      };
      const Real hw = std::max<Real>(cfg.gradeTransitionWidth, 1.0e-6);
      const Real f1 = ss((s - (1.0 / 3.0 - hw)) / (2.0 * hw)); // intima -> media
      const Real f2 = ss((s - (2.0 / 3.0 - hw)) / (2.0 * hw)); // media -> adventitia
      return cfg.gradeIntima + (cfg.gradeMedia - cfg.gradeIntima) * f1 +
        (cfg.gradeAdventitia - cfg.gradeMedia) * f2;
    };

    Rodin::Examples::Heart::GradedYeoh law(
      Solid::Yeoh(yeohC1, yeohC2, yeohC3, yeohKappa), wallGrade);
    Solid::InternalVirtualWorkTangent solidTangent(law, d, w, dState);
    Solid::InternalVirtualWorkResidual solidInternal(law, w, dState);

    // Current fluid solution (read back after each Oseen solve).
    PETSc::Variational::GridFunction uCur(uh);
    PETSc::Variational::GridFunction pCur(ph);
    uCur = zero;
    pCur = 0.0;

    // Fluid Cauchy traction on Gamma_FSI (current config):
    //   tractionFSI = p n - mu(grad u + grad u^T) n = -sigma_f n_f,
    const auto gradUfsi = Jacobian(uCur);
    const auto strainRateFsi = gradUfsi + Transpose(gradUfsi);
    const auto symUfsi = 0.5 * strainRateFsi;
    const auto shearFsi = Sqrt(gammaReg * gammaReg + 2.0 * Dot(symUfsi, symUfsi));
    const auto muFsi = cy.muInf +
      deltaMu * Pow(1.0 + Pow(cy.lambda * shearFsi, cy.yasuda), (cy.n - 1.0) / cy.yasuda);

    // PHYSICAL (unscaled) fluid traction t_f = -sigma_f n_f: the fluid's OWN
    // lagged Robin datum sigma_f^{lag} n and the traction output field.
    const auto tractionFSI =
      (1.0 * pCur) * normalFluid - (1.0 * muFsi) * Mult(strainRateFsi, normalFluid);

    // ------------------------------------------------------------------
    // Volume-recovered wall shear stress (fixes the surface-gradient bug).
    // ------------------------------------------------------------------

    PETSc::Variational::TrialFunction gradRecTrial(uh);
    PETSc::Variational::TestFunction gradRecTest(uh);
    const auto jacRow0 = VectorFunction(Component(Jacobian(uCur), 0, 0),
      Component(Jacobian(uCur), 0, 1), Component(Jacobian(uCur), 0, 2));
    const auto jacRow1 = VectorFunction(Component(Jacobian(uCur), 1, 0),
      Component(Jacobian(uCur), 1, 1), Component(Jacobian(uCur), 1, 2));
    const auto jacRow2 = VectorFunction(Component(Jacobian(uCur), 2, 0),
      Component(Jacobian(uCur), 2, 1), Component(Jacobian(uCur), 2, 2));
    Problem gradRecProj0(gradRecTrial, gradRecTest);
    gradRecProj0 = Integral(gradRecTrial, gradRecTest) - Integral(jacRow0, gradRecTest);
    Problem gradRecProj1(gradRecTrial, gradRecTest);
    gradRecProj1 = Integral(gradRecTrial, gradRecTest) - Integral(jacRow1, gradRecTest);
    Problem gradRecProj2(gradRecTrial, gradRecTest);
    gradRecProj2 = Integral(gradRecTrial, gradRecTest) - Integral(jacRow2, gradRecTest);

    // Build the WSS from the recovered gradient rows using only vector ops
    // (Component scalars don't compose under +).  (grad u . n)_i = gradRec_i . n
    // is the wall-NORMAL directional derivative of u; at a no-slip wall the
    // transpose part (grad u^T . n) vanishes (u_n ~ 0 -> d u_n/dn ~ 0), so this
    // equals the full strain-rate traction in the TANGENTIAL direction.  The
    // LOCAL Carreau-Yasuda viscosity muFsi is used (same expression as the
    // momentum equation): a constant muInf decouples the WSS from the
    // rheology and under-predicts tau_w exactly in the low-shear
    // recirculation zones where mu departs from muInf (see
    // CoupledLV0DCoronary3D::computeWallShear).
    const auto gradUn0 = Dot(gradRec0, normalFluid);
    const auto gradUn1 = Dot(gradRec1, normalFluid);
    const auto gradUn2 = Dot(gradRec2, normalFluid);
    const auto tracRec =
      VectorFunction(muFsi * gradUn0, muFsi * gradUn1, muFsi * gradUn2);
    const auto wallStressRec = tracRec - Dot(tracRec, normalFluid) * normalFluid;

    // Areal stretch J_a = A_t/A_0 at the CURRENT iterate dIter: pulls
    // per-current-area interface data back to the reference solid surface.
    auto arealStretchAt = [&](const Point& xs) -> Real {
      Real stretch = 1.0;
      const auto& verts = xs.getPolytope().getVertices();
      if (dim == 3 && verts.size() == 3)
      {
        Math::SpatialPoint X0 = meshSolid.getVertexCoordinates(verts[0]);
        Math::SpatialPoint X1 = meshSolid.getVertexCoordinates(verts[1]);
        Math::SpatialPoint X2 = meshSolid.getVertexCoordinates(verts[2]);

        Math::SpatialPoint x0 = X0, x1 = X1, x2 = X2;
        for (Index c = 0; c < 3; ++c)
        {
          x0(c) += dIter[dh.getGlobalIndex({0, verts[0]}, c)];
          x1(c) += dIter[dh.getGlobalIndex({0, verts[1]}, c)];
          x2(c) += dIter[dh.getGlobalIndex({0, verts[2]}, c)];
        }

        const auto triArea = [](const Math::SpatialPoint& a, const Math::SpatialPoint& b,
                               const Math::SpatialPoint& c) -> Real {
          const Math::SpatialPoint e1 = b - a;
          const Math::SpatialPoint e2 = c - a;
          const Real nx = e1(1) * e2(2) - e1(2) * e2(1);
          const Real ny = e1(2) * e2(0) - e1(0) * e2(2);
          const Real nz = e1(0) * e2(1) - e1(1) * e2(0);
          return 0.5 * std::sqrt(nx * nx + ny * ny + nz * nz);
        };

        const Real A0 = triArea(X0, X1, X2);
        const Real At = triArea(x0, x1, x2);
        if (A0 > 0.0)
          stretch = At / A0;
      }
      return stretch;
    };

    auto fluidStress = VectorFunction(dim, [&](const Point& xs) {
      const Point xf = forwardSolidPointToFluid(xs, meshFluid, interfaceMap);

      Math::SpatialVector<Real> value(dim);
      value.setZero();
      const auto force = tractionTransfer(xf);
      for (Index i = 0; i < static_cast<Index>(dim); ++i)
        value(i) = force(i);

      if (!isFiniteVec(value))
      {
        static bool reported = false;
        if (!reported)
        {
          reported = true;
          const auto& pc = xs.getPhysicalCoordinates();
          Alert::Warning() << "fluidStress: non-finite cross-mesh traction sample at xs=("
                           << pc.transpose() << "); using zero there (warned once)."
                           << Alert::Raise;
        }
        value.setZero();
        return value;
      }

      // Pull the current-config traction back to the reference configuration.
      value *= arealStretchAt(xs);
      return value;
    });

    // ----------------------------------------------------------------------
    // Fluid Oseen / BDF1-ALE problem
    // ----------------------------------------------------------------------
    // Interface velocity entering the fluid Robin data alpha*u_s.  Two choices:
    //   meshConsistentInterfaceVelocity == true  (recommended):
    //     u_s = meshVelocity(xf) -- the ALE mesh velocity at the wall, i.e. the
    //     BDF1 rate (d^{n+1}-d^n)/dt that the fluid domain boundary ACTUALLY
    //     moves with (GCL).
    auto interfaceSolidVelocity = VectorFunction(dim, [&](const Point& xf) {
      Math::SpatialVector<Real> value(dim);
      value.setZero();

      if (cfg.meshConsistentInterfaceVelocity)
      {
        const auto w = meshVelocity(xf);
        for (Index i = 0; i < static_cast<Index>(dim); ++i)
          value(i) = w(i);
      }
      else
      {
        const Point xs = forwardFluidPointToSolid(xf, meshSolid, interfaceMap);
        const auto us = solidVelocity(xs);
        for (Index i = 0; i < static_cast<Index>(dim); ++i)
          value(i) = us(i);
      }

      if (!isFiniteVec(value))
      {
        static bool reported = false;
        if (!reported)
        {
          reported = true;
          Alert::Warning()
            << "interfaceSolidVelocity: non-finite interface velocity sample; "
            << "using zero there (warned once)." << Alert::Raise;
        }
        value.setZero();
      }

      return value;
    });

    const auto transportLag = uOld - meshVelocity;
    const auto convU = Mult(Jacobian(u), transportLag);
    // Conservative two-mesh BDF1 mass split requires the FULL geometric
    // term -rho(div w)(u.v); with the 1/2 factor below the coefficient is
    // (1/2)div(u^n) - div(w)  (Temam part + geometric companion).
    const auto divGeomTemam = Div(uOld) - 2.0 * Div(meshVelocity);

    const auto duNormal = Dot(u, normalFluid) * normalFluid;

    const auto duTangential = u - duNormal;

    const auto symU = 0.5 * (Jacobian(u) + Transpose(Jacobian(u)));
    const auto symV = 0.5 * (Jacobian(v) + Transpose(Jacobian(v)));
    const auto symLag = 0.5 * (Jacobian(uOld) + Transpose(Jacobian(uOld)));
    const auto shearLag = Sqrt(gammaReg * gammaReg + 2.0 * Dot(symLag, symLag));
    const auto muLag = cy.muInf +
      deltaMu * Pow(1.0 + Pow(cy.lambda * shearLag, cy.yasuda), (cy.n - 1.0) / cy.yasuda);

    const auto outletBeta = Max(-Dot(transportLag, normalFluid), 0.0);
    const auto inletBeta = Max(Dot(transportLag, normalFluid), 0.0);
    const auto outletBackflow =
      0.5 * cfg.outletBackflowStabilization * cfg.fluidDensity * outletBeta;
    const auto inletBackflow =
      0.5 * cfg.inletBackflowStabilization * cfg.fluidDensity * inletBeta;

    using namespace Rodin::Examples::Heart;

    // ALE convecting velocity u^n - w as a grid function (refreshed each
    // coupling iterate, since the mesh velocity w changes with the iterate).
    PETSc::Variational::GridFunction uConv(uh);
    uConv = uOld;

    PressureFES tauFes(std::integral_constant<size_t, 1>{}, meshFluid);
    PETSc::Variational::TestFunction vmsScalarTest(tauFes);
    PETSc::Variational::TrialFunction vmsPiTilde(tauFes); // pi~ (projected)

    PETSc::Variational::TrialFunction vmsUp(uh); // Pi[(grad uConv) uConv]
    PETSc::Variational::TestFunction vmsVp(uh);
    PETSc::Variational::TrialFunction vmsSub(uh); // dynamic subscale u'^{n+1}
    PETSc::Variational::GridFunction vmsSubOld(uh); // subscale history u'^n
    vmsSubOld = zero;

    // Frozen convective acceleration (grad uConv) uConv.
    const auto vmsConvectionTarget = Mult(Jacobian(uConv), uConv);

    auto viscosityAt = [&](const Point& pp) -> Real {
      const auto sym = 0.5 * (Jacobian(uOld) + Transpose(Jacobian(uOld)));
      const Real shear =
        std::sqrt(gammaReg * gammaReg + 2.0 * Dot(sym, sym).getValue(pp));
      return cy.muInf +
        deltaMu *
        std::pow(1.0 + std::pow(cy.lambda * shear, cy.yasuda), (cy.n - 1.0) / cy.yasuda);
    };

    // Shared stabilization parameter tau1 (Codina), c1 = 4, c2 = 2, k = 1
    // (P1), h_K = |K|^{1/d}, nu = mu(x)/rho local; the transport speed is the
    // ALE relative velocity |u^n - w|.
    auto tau1At = [&](const Point& pp) -> Real {
      const auto uc = uConv.getValue(pp);
      const Real nu = viscosityAt(pp) / cfg.fluidDensity;
      const Real hK =
        std::pow(pp.getPolytope().getMeasure(), 1.0 / pp.getPolytope().getDimension());
      const Real speed = std::sqrt(Math::dot(uc, uc));
      return 1.0 / (4.0 * nu / (hK * hK) + 2.0 * speed / hK);
    };

    // Convective subscale parameter tau_K = vmsScale/(rho/dt + rho/tau1),
    // evaluated pointwise (fed straight into the VMS integrators).
    auto vmsTauAt = [&](const Point& pp) -> Real {
      return cfg.vmsScale *
        (1.0 / (cfg.fluidDensity / dt + cfg.fluidDensity / tau1At(pp)));
    };
    RealFunction vmsTauFn = [&](const Point& pp) -> Real { return vmsTauAt(pp); };

    // Grad-div parameter tau_C = gradDivScale * rho * h^2/(4 tau1) and its
    // square root, exact by construction (tau_C = (sqrt tau_C)^2 pointwise).
    auto sqrtTauCAt = [&](const Point& pp) -> Real {
      const Real hK =
        std::pow(pp.getPolytope().getMeasure(), 1.0 / pp.getPolytope().getDimension());
      return std::sqrt(std::max<Real>(
        0.0, cfg.gradDivScale * cfg.fluidDensity * hK * hK / (4.0 * tau1At(pp))));
    };
    RealFunction sqrtTauCFn = [&](const Point& pp) -> Real { return sqrtTauCAt(pp); };
    RealFunction tauCFn = [&](const Point& pp) -> Real {
      const Real s = sqrtTauCAt(pp);
      return s * s;
    };

    // PSPG / grad-p coefficient tau_p = pgpScale * tau1/rho, pointwise.
    RealFunction vmsTauPFn = [&](const Point& pp) -> Real {
      return cfg.pgpScale * tau1At(pp) / cfg.fluidDensity;
    };

    // Dynamic subscale update (L2-projected into vmsSub: it must be STORED,
    // it carries the history u'^n).
    auto vmsSubUpdate =
      VectorFunction(dim, [&](const Point& pp) -> Math::SpatialVector<Real> {
        const auto conv = vmsConvectionTarget.getValue(pp);
        const auto proj = vmsUp.getSolution().getValue(pp);
        const auto old = vmsSubOld.getValue(pp);
        const Real tau = vmsTauAt(pp);

        Math::SpatialVector<Real> out(dim);
        for (Index c = 0; c < static_cast<Index>(dim); ++c)
          out(c) = tau * cfg.fluidDensity * (1.0 / dt * old(c) - (conv(c) - proj(c)));
        return out;
      });

    // The two orthogonal-subscale L2 projections + the subscale history
    // (mass matrices reassembled each iterate as the mesh moves).
    Problem vmsL2Conv(vmsUp, vmsVp);
    vmsL2Conv = Integral(vmsUp, vmsVp) - Integral(vmsConvectionTarget, vmsVp);

    Problem vmsSubProj(vmsSub, vmsVp);
    vmsSubProj = Integral(vmsSub, vmsVp) - Integral(vmsSubUpdate, vmsVp);

    // pi~ = Pi( sqrt(tau_C) div u^n ).  The SAME sqrt(tau_C) multiplies
    // div(v) in the linear term and squares into the implicit coefficient;
    // lagging one of the two by a step leaves the explicit half larger than
    // the implicit one on the step where the viscosity drops (see Atrium).
    Problem vmsPiTildeProj(vmsPiTilde, vmsScalarTest);
    vmsPiTildeProj = Integral(vmsPiTilde, vmsScalarTest) -
      Integral(sqrtTauCFn * Div(uOld), vmsScalarTest);

    // Mass-matrix solves are SPD and cheap with CG + Jacobi; never the global
    // direct (MUMPS LU) solver used for the coupled flow system.  The prefixed
    // defaults are set in main() and remain overridable on the command line.
    auto solveMass = [](auto& problem) {
      problem.assemble();
      Solver::KSP ksp(problem);
      ksp.setPrefix("coronary_mass_");
      ksp.solve();
    };

    // The Neumann traction at an outlet is
    // p_out = p_c + R_a Phi_a Q, and for a flat profile the resistive part is
    //   int_G R_a Phi_a Q (v.n) = R_a Phi_a A int_G (u.n)(v.n),
    // symmetric and positive semidefinite for any profile.  Assembling it here
    // instead of lagging it is what makes the 0D-3D coupling unconditionally
    // stable: the amplification factor goes from |1 - R dt / L_3D| to
    // 1/(1 + R dt / L_3D) < 1 for any dt and any R.  Phi_a is refreshed by
    // updateOutlet0D after each step and read here at reassembly (the flow
    // form is reassembled every coupling iterate because the mesh moves).
    auto outletZAt = [&](size_t i) -> Real {
      const auto& bc = wk.at(BoundaryFluid::Outlets[i]);
      return cfg.outletResistanceScale * bc.Ra * bc.muA * bc.area;
    };
    auto zFn0 = RealFunction([&](const Point&) { return outletZAt(0); });
    auto zFn1 = RealFunction([&](const Point&) { return outletZAt(1); });
    auto zFn2 = RealFunction([&](const Point&) { return outletZAt(2); });
    auto zFn3 = RealFunction([&](const Point&) { return outletZAt(3); });
    auto zFn4 = RealFunction([&](const Point&) { return outletZAt(4); });
    auto zFn5 = RealFunction([&](const Point&) { return outletZAt(5); });

    // BDF1 mass split: (rho/dt)[(u,v)_{n+1} - (u^n,v)_n]; the implicit part
    // lives in 'flow', the explicit u^n part is 'massOld', assembled on the
    // PREVIOUS configuration and injected into the RHS at solve time.
    Problem flow(u, p, v, q);
    flow = (cfg.fluidDensity / dt) * Integral(u, v) +
      cfg.fluidDensity * Integral(Dot(convU, v)) +
      0.5 * cfg.fluidDensity * Integral(divGeomTemam * Dot(u, v)) +
      // Projected-VMS convective stabilization (bilinear + subtracted
      // linear); tau is evaluated POINTWISE (vmsTauFn folds cfg.vmsScale so
      // both terms vanish when vmsScale == 0).
      VMSConvectionBilinearIntegrator(u, v, uConv, vmsTauFn, cfg.fluidDensity) -
      VMSConvectionLinearIntegrator(v, vmsSub.getSolution(), uConv, vmsUp.getSolution(),
        vmsTauFn, cfg.fluidDensity, dt) +
      // Orthogonal grad-div / pressure subscale p~, IMEX split:
      //   + int tau_C (div u^{n+1})(div v)        [implicit, LHS]
      //   - int sqrt(tau_C) pi~ (div v)           [lagged, RHS]
      // with pi~ = Pi( sqrt(tau_C) div u^n ); tau_C = (sqrt tau_C)^2 exactly,
      // both evaluated pointwise at the SAME time level.
      VMSGradDivBilinearIntegrator(u, v, tauCFn) -
      VMSGradDivLinearIntegrator(v, vmsPiTilde.getSolution(), sqrtTauCFn) +
      2.0 * Integral(muLag * symU, symV) - Integral(p, Div(v)) + Integral(Div(u), q) +
      cfg.pressurePenalty * Integral(p, q) +
      // Pressure-gradient (PSPG / grad-p) stabilization, REQUIRED for the
      // equal-order P1/P1 pair; tau_p pointwise.
      Integral(vmsTauPFn * Grad(p), Grad(q)) +
      BoundaryIntegral(inletBackflow * Dot(u, v)).over(BoundaryFluid::Inlet) +
      BoundaryIntegral(outletBackflow * Dot(u, v))
        .over(BoundaryFluid::Outlets[0], BoundaryFluid::Outlets[1],
          BoundaryFluid::Outlets[2], BoundaryFluid::Outlets[3], BoundaryFluid::Outlets[4],
          BoundaryFluid::Outlets[5]) +
      BoundaryIntegral(pin * Dot(v, normalFluid)).over(BoundaryFluid::Inlet) +
      BoundaryIntegral(pout0 * Dot(v, normalFluid)).over(BoundaryFluid::Outlets[0]) +
      BoundaryIntegral(pout1 * Dot(v, normalFluid)).over(BoundaryFluid::Outlets[1]) +
      BoundaryIntegral(pout2 * Dot(v, normalFluid)).over(BoundaryFluid::Outlets[2]) +
      BoundaryIntegral(pout3 * Dot(v, normalFluid)).over(BoundaryFluid::Outlets[3]) +
      BoundaryIntegral(pout4 * Dot(v, normalFluid)).over(BoundaryFluid::Outlets[4]) +
      BoundaryIntegral(pout5 * Dot(v, normalFluid)).over(BoundaryFluid::Outlets[5])
      // Implicit outlet resistance R_a Phi_a A (u.n)(v.n), per outlet.  The
      // resistive part of p_out is NOT in the Neumann datum (pout = p_c), so
      // this term IS the arteriolar resistance, fully implicit -- no lagged
      // cancellation needed, unconditionally stable (CoupledLV0DCoronary3D).
      + BoundaryIntegral(zFn0 * Dot(Dot(u, normalFluid) * normalFluid, v))
          .over(BoundaryFluid::Outlets[0]) +
      BoundaryIntegral(zFn1 * Dot(Dot(u, normalFluid) * normalFluid, v))
        .over(BoundaryFluid::Outlets[1]) +
      BoundaryIntegral(zFn2 * Dot(Dot(u, normalFluid) * normalFluid, v))
        .over(BoundaryFluid::Outlets[2]) +
      BoundaryIntegral(zFn3 * Dot(Dot(u, normalFluid) * normalFluid, v))
        .over(BoundaryFluid::Outlets[3]) +
      BoundaryIntegral(zFn4 * Dot(Dot(u, normalFluid) * normalFluid, v))
        .over(BoundaryFluid::Outlets[4]) +
      BoundaryIntegral(zFn5 * Dot(Dot(u, normalFluid) * normalFluid, v))
        .over(BoundaryFluid::Outlets[5]) +
      cfg.inletImpedance *
        BoundaryIntegral(Dot(Dot(u, normalFluid) * normalFluid, v))
          .over(BoundaryFluid::Inlet) +
      cfg.inletTangentialDamping *
        BoundaryIntegral(Dot(duTangential, v)).over(BoundaryFluid::Inlet)
      // Robin-Robin transmission (fluid side), Burman et al. (2025):
      //   sigma_f n + alpha u = alpha d_dot_s + lambda^{k-1},
      // with lambda^{k-1} = tractionFSI at the PREVIOUS correction (= the
      // previous time step at the first correction, i.e. lambda^{n-1} in
      // the loose kappa = 0 scheme).
      + robinAlpha * BoundaryIntegral(u, v).over(BoundaryFluid::FSI) -
      robinAlpha * BoundaryIntegral(interfaceSolidVelocity, v).over(BoundaryFluid::FSI) +
      BoundaryIntegral(tractionFSI, v).over(BoundaryFluid::FSI)
      // Interface convective stabilization (Burman et al. 2025, eq. 13):
      // -(rho/2) (transportLag.n)(u.v) on Sigma
      - 0.5 * cfg.fluidDensity *
        BoundaryIntegral(Dot(transportLag, normalFluid) * Dot(u, v))
          .over(BoundaryFluid::FSI)
      // Strong no-slip on the cap rings: the one-element FSI band touching
      // the inlet/outlet caps is pinned to zero, consistent with the solid
      // ring clamp.
      + DirichletBC(u, zero).on(BoundaryFluid::FSIRing);

    PETSc::Variational::TestFunction vMass(uh);
    LinearForm<VelocityFES, ::Vec> massOld(vMass);
    massOld = (cfg.fluidDensity / dt) * Integral(uOld, vMass);

    // ALE Dirichlet datum: solid displacement iterate sampled cross-mesh.
    auto interfaceSolidDisplacement = VectorFunction(dim, [&](const Point& xf) {
      const Point xs = forwardFluidPointToSolid(xf, meshSolid, interfaceMap);

      Math::SpatialVector<Real> value(dim);
      value.setZero();

      const auto ds = dIter(xs);
      for (Index i = 0; i < static_cast<Index>(dim); ++i)
        value(i) = ds(i);

      if (!isFiniteVec(value))
      {
        static bool reported = false;
        if (!reported)
        {
          reported = true;
          Alert::Warning() << "interfaceSolidDisplacement: non-finite cross-mesh sample; "
                           << "using zero there (warned once)." << Alert::Raise;
        }
        value.setZero();
      }

      return value;
    });

    // Harmonic ALE lift (reference config).
    PETSc::Variational::TrialFunction dMove(uh);
    PETSc::Variational::TestFunction vMove(uh);

    // Jacobian/element-size stiffening weight w(x) = (aleRefSize/h_K)^power:
    RealFunction aleStiffFn = [&](const Point& pp) -> Real {
      if (cfg.aleStiffPower <= 0.0)
        return 1.0;
      const Real hK =
        std::pow(pp.getPolytope().getMeasure(), 1.0 / pp.getPolytope().getDimension());
      return std::pow(cfg.aleRefSize / std::max(hK, 1.0e-30), cfg.aleStiffPower);
    };

    Problem ale(dMove, vMove);
    ale = Integral(aleStiffFn * Jacobian(dMove), Jacobian(vMove)) +
      DirichletBC(dMove, interfaceSolidDisplacement).on(BoundaryFluid::FSI) +
      DirichletBC(dMove, zero).on(BoundaryFluid::Inlet, BoundaryFluid::FSIRing);

    // Solid Robin data, per reference area:
    //   J_a [ rVC (dState - dPred) + alpha vPred - alpha u_f^{k-1} ],
    // with u_f sampled from the projected uWall (latest fluid iterate).
    auto robinInterfaceData =
      VectorFunction(
        dim, [&](const Point& xs) {
          const Point xf = forwardSolidPointToFluid(xs, meshFluid, interfaceMap);
          //const Real Ja = arealStretchAt(xs);

          Math::SpatialVector<Real> value(dim);
          value.setZero();

          // Guarded cross-mesh fluid velocity sample (see fluidStress note).
          Math::SpatialVector<Real> uf(dim);
          uf.setZero();
          {
            const auto ufRaw = uWall(xf);
            for (Index i = 0; i < static_cast<Index>(dim); ++i)
              uf(i) = ufRaw(i);
            if (!isFiniteVec(uf))
            {
              static bool reported = false;
              if (!reported)
              {
                reported = true;
                const auto& pc = xs.getPhysicalCoordinates();
                Alert::Warning()
                  << "robinInterfaceData: non-finite cross-mesh velocity sample "
                  << "at xs=(" << pc.transpose() << "); using zero there (warned once)."
                  << Alert::Raise;
              }
              uf.setZero();
            }
          }
          const auto vp = vPred(xs);
          const auto dS = dState(xs);
          const auto dP = dPred(xs);
          for (Index i = 0; i < static_cast<Index>(dim); ++i)
            value(i) = (robinVelocityCoeff * (dS(i) - dP(i)) + robinAlpha * vp(i) -
              robinAlpha * uf(i));

          return value;
        });

    // ----------------------------------------------------------------------
    // Solid Newmark / NeoHookean problem (nonlinear: SNES on the increment
    // etaState, with dState = dOld + etaState).
    // ----------------------------------------------------------------------
    // 0D heart-motion penalty (CoronarySolid style)
    const auto normalSolid = BoundaryNormal(meshSolid);
    Real disp0D = 0.0;
    Real disp0DOffset = 0.0; // 0D displacement at the dynamic handoff (step 1)
    auto disp0DFn = RealFunction([&](const Point&) { return -disp0D; });
    const Real heartK = cfg.heartDisplacementPenalty;

    // Kelvin-Voigt viscous damping force  int eta_s grad(d_dot):grad(w), with
    // the Newmark velocity  d_dot = vPred + solidVelocityCoeff (dState - dPred).
    const Real solidViscImpl = cfg.solidViscosity * solidVelocityCoeff;

    const Real a = cfg.aViscCondition;
    const Real b = cfg.bViscCondition;
    const Real aVel = b * gammaN / (betaN * dt);

    Problem solid(d, w);
    solid = solidMass * Integral(d, w) + solidTangent + solidMass * Integral(dState, w) -
      solidMass * Integral(dPred, w) +
      solidInternal
      // viscous (strain-rate) damping, Newmark-consistent:
      + solidViscImpl * Integral(Jacobian(d), Jacobian(w)) +
      solidViscImpl * Integral(Jacobian(dState), Jacobian(w)) -
      solidViscImpl * Integral(Jacobian(dPred), Jacobian(w)) +
      cfg.solidViscosity * Integral(Jacobian(vPred), Jacobian(w)) +
      DirichletBC(d, zero).on(BoundarySolid::Inlet) +
      DirichletBC(d, zero).on(BoundarySolid::FSIRing) +
      a *
        BoundaryIntegral(d, w).over(BoundarySolid::Outlets[0], BoundarySolid::Outlets[1],
          BoundarySolid::Outlets[2], BoundarySolid::Outlets[3], BoundarySolid::Outlets[4],
          BoundarySolid::Outlets[5]) +
      aVel *
        BoundaryIntegral(d, w).over(BoundarySolid::Outlets[0], BoundarySolid::Outlets[1],
          BoundarySolid::Outlets[2], BoundarySolid::Outlets[3], BoundarySolid::Outlets[4],
          BoundarySolid::Outlets[5]) +
      a *
        BoundaryIntegral(dState, w).over(BoundarySolid::Outlets[0],
          BoundarySolid::Outlets[1], BoundarySolid::Outlets[2], BoundarySolid::Outlets[3],
          BoundarySolid::Outlets[4], BoundarySolid::Outlets[5]) +
      aVel *
        BoundaryIntegral(dState, w).over(BoundarySolid::Outlets[0],
          BoundarySolid::Outlets[1], BoundarySolid::Outlets[2], BoundarySolid::Outlets[3],
          BoundarySolid::Outlets[4], BoundarySolid::Outlets[5]) -
      aVel *
        BoundaryIntegral(dPred, w).over(BoundarySolid::Outlets[0],
          BoundarySolid::Outlets[1], BoundarySolid::Outlets[2], BoundarySolid::Outlets[3],
          BoundarySolid::Outlets[4], BoundarySolid::Outlets[5]) +
      b *
        BoundaryIntegral(vPred, w).over(BoundarySolid::Outlets[0],
          BoundarySolid::Outlets[1], BoundarySolid::Outlets[2], BoundarySolid::Outlets[3],
          BoundarySolid::Outlets[4], BoundarySolid::Outlets[5])
      // 0D heart motion: weakly impose d.n = disp_0D on attributes 110..120
      // (k (d.n)(w.n) [tangent] + k (dState.n)(w.n) - k disp_0D (w.n)
      //  [residual] = k (dState.n - disp_0D)(w.n)).
      + heartK *
        BoundaryIntegral(Dot(d, normalSolid), Dot(w, normalSolid))
          .over(BoundarySolid::Contact[0]) +
      heartK *
        BoundaryIntegral(Dot(dState, normalSolid), Dot(w, normalSolid))
          .over(BoundarySolid::Contact[0]) -
      heartK *
        BoundaryIntegral(disp0DFn * l.getSolution(), Dot(w, normalSolid))
          .over(BoundarySolid::Contact[0])
      // Robin-Robin transmission (solid side), per current area (J_a):
      //   sigma_s n_s + alpha d_dot = alpha u_f^{lag} + t_f^{lag},
      //   d_dot = vPred + solidVelocityCoeff (dState - dPred)  (Newmark).
      + robinVelocityCoeff * BoundaryIntegral(Dot(d, w)).over(BoundarySolid::FSI) +
      BoundaryIntegral(robinInterfaceData, w).over(BoundarySolid::FSI)
      // Neumann load: the traction exerted by the fluid on the solid is
      //   t_f = sigma_f n_s = -sigma_f n_f = fluidStress
      - BoundaryIntegral(fluidStress, w).over(BoundarySolid::FSI);

    solid.assemble();
    Solver::KSP kspSolid(solid);
    Solver::SNES snes(kspSolid);
    snes.setTolerances(1.0e-10, 1.0e-8, 1.0e-10, 50, 10000);
    snes.setStateUpdate([&](const PETSc::Math::Vector& state) {
      etaState.setData(state, 0);
      dState = dOld;
      dState += etaState;
    });

    Problem laplacian(l, t);
    laplacian = Integral(Grad(l), Grad(t)) +
      DirichletBC(l, RealFunction(0.0)).on(BoundarySolid::Inlet) +
      DirichletBC(l, RealFunction(0.95)).on(BoundarySolid::Outlets[3]) +
      DirichletBC(l, RealFunction(0.85)).on(BoundarySolid::Outlets[4]) +
      DirichletBC(l, RealFunction(0.75))
        .on(BoundarySolid::Outlets[0], BoundarySolid::Outlets[1],
          BoundarySolid::Outlets[2]) +
      DirichletBC(l, RealFunction(0.45)).on(BoundarySolid::Outlets[5]);

    if (isRoot)
    {
      PETSc::Variational::GridFunction oneSolid(lh);
      oneSolid = 1.0;
      LinearForm<LaplacianFES, ::Vec> solidArea(t);
      auto faceArea = [&](Attribute tag) -> Real {
        solidArea = BoundaryIntegral(oneSolid, t).over(tag);
        solidArea.assemble();
        return solidArea(oneSolid);
      };
      Alert::Info() << "  [solid-bdr] Inlet(" << BoundarySolid::Inlet
                    << ") area=" << faceArea(BoundarySolid::Inlet) << "  FSI("
                    << BoundarySolid::FSI << ") area=" << faceArea(BoundarySolid::FSI)
                    << "  FSIRing(" << BoundarySolid::FSIRing
                    << ") area=" << faceArea(BoundarySolid::FSIRing) << Alert::Raise;
      for (size_t i = 0; i < BoundarySolid::Outlets.size(); ++i)
        Alert::Info() << "  [solid-bdr] Outlet(" << BoundarySolid::Outlets[i]
                      << ") area=" << faceArea(BoundarySolid::Outlets[i]) << Alert::Raise;
    }

    laplacian.assemble();
    Solver::KSP(laplacian).solve();

    // Transmural coordinate: Delta xi = 0, xi = 0 on the luminal (FSI)
    // surface, xi = 1 on the adventitial surface (Outer).  Natural BC on the
    // caps gives the linear through-thickness profile.  Must be solved BEFORE
    // the prestress: the graded law reads xi at every quadrature point.
    Problem thickness(xiT, xiTest);
    thickness = Integral(Grad(xiT), Grad(xiTest)) +
      DirichletBC(xiT, RealFunction(0.0)).on(BoundarySolid::FSI) +
      DirichletBC(xiT, RealFunction(1.0))
        .on(BoundarySolid::Outer[0], BoundarySolid::Outer[1]);
    thickness.assemble();
    Solver::KSP(thickness).solve();
    xi.setData(xiT.getSolution().getData());
    if (isRoot)
      Alert::Info() << "  [grade] transmural xi solved; multipliers"
                    << "  intima=" << cfg.gradeIntima << "  media=" << cfg.gradeMedia
                    << "  adventitia=" << cfg.gradeAdventitia
                    << "  (transition half-width " << cfg.gradeTransitionWidth << ")"
                    << Alert::Raise;

    xdmf_laplacian.add("xi", xi);
    xdmf_laplacian.write().flush();
    xdmf_laplacian.close();

    if (cfg.prestressSteps > 0)
    {
      // Prestress statically to prestressFraction*par; the dynamic loop ramps
      // the remaining (1 - prestressFraction) over prestressRampSteps.
      const Real p0 = cfg.prestressFraction * model.getState().par;
      Real prestressPressure = 0.0;

      PETSc::Variational::TrialFunction dPre(dh);
      PETSc::Variational::TestFunction wPre(dh);
      Solid::InternalVirtualWorkTangent preTangent(law, dPre, wPre, dState);
      Solid::InternalVirtualWorkResidual preInternal(law, wPre, dState);

      // Follower pressure (exact deformed-surface load + consistent
      // tangent): full quadratic Newton up to total pressure.
      Solid::FollowerPressureForce preLoad(prestressPressure, wPre, dState);
      preLoad.over(BoundarySolid::FSI);
      Solid::FollowerPressureTangent preLoadK(prestressPressure, dPre, wPre, dState);
      preLoadK.over(BoundarySolid::FSI);

      // The outlets are NOT clamped: they carry the same elastic tethering
      // spring a used by the dynamic problem (tangent on dPre + residual on
      // dState), so the prestressed state hands off to the dynamics without a
      // jump when the springs take over.  Only the inlet and its ring band
      // are pinned (the anchor of the tree).
      Problem prestress(dPre, wPre);
      prestress = preTangent + preInternal + preLoadK + preLoad +
        a *
          BoundaryIntegral(dPre, wPre)
            .over(BoundarySolid::Outlets[0], BoundarySolid::Outlets[1],
              BoundarySolid::Outlets[2], BoundarySolid::Outlets[3],
              BoundarySolid::Outlets[4], BoundarySolid::Outlets[5]) +
        a *
          BoundaryIntegral(dState, wPre)
            .over(BoundarySolid::Outlets[0], BoundarySolid::Outlets[1],
              BoundarySolid::Outlets[2], BoundarySolid::Outlets[3],
              BoundarySolid::Outlets[4], BoundarySolid::Outlets[5]) +
        DirichletBC(dPre, zero).on(BoundarySolid::Inlet, BoundarySolid::FSIRing);

      prestress.assemble();
      Solver::KSP kspPre(prestress);
      Solver::SNES snesPre(kspPre);
      snesPre.setTolerances(1.0e-10, 1.0e-8, 1.0e-10, 50, 10000);
      snesPre.setStateUpdate([&](const PETSc::Math::Vector& state) {
        etaState.setData(state, 0);
        dState = dOld; // dOld stays zero: etaState accumulates the full lift
        dState += etaState;
      });

      for (size_t k = 1; k <= cfg.prestressSteps; ++k)
      {
        prestressPressure =
          (static_cast<Real>(k) / static_cast<Real>(cfg.prestressSteps)) * p0;
        snesPre.solve();
        if (!snesPre.converged())
        {
          if (isRoot)
            std::cerr << "Prestress SNES failed at increment " << k << " / "
                      << cfg.prestressSteps << "; continuing with the last "
                      << "converged (partial) prestress state.\n";
          break;
        }
      }
      if (isRoot)
        Alert::Info() << "Prestressed wall to " << prestressPressure << " Pa in "
                      << cfg.prestressSteps << " increment(s)" << Alert::Raise;

      // Commit the prestressed state as the dynamic initial condition (zero
      // velocity / acceleration are already set) and lift the fluid mesh.
      dOld.setData(dState.getData());
      dIter.setData(dState.getData());
      restoreMeshToReference(meshFluid, referenceVertices);
      ale.assemble();
      Solver::KSP(ale).solve();
      aleDisp.setData(dMove.getSolution().getData());
      aleDispOld.setData(aleDisp.getData());
      moveMeshWithVertexDisplacement(meshFluid, referenceVertices, uh, aleDisp);

      // Seed the fluid pressure at the prestress level (p0 = prestressFraction*par)
      // so the transferred wall traction balances the prestressed wall at the
      // handoff; the dynamic ramp then takes the loads up to full par.
      pOld = p0;
      pCur = p0;
      tractionTransfer.project(Region::Faces, tractionFSI, BoundaryFluid::FSI);
    }

    // Interface flux functional q_flux = \int_Gamma (u . n) for the RCR/0D
    // Windkessel coupling at the inlet and outlets.
    using FluxLinearForm = LinearForm<PressureFES, ::Vec>;
    PETSc::Variational::TestFunction qFlux(ph);
    FluxLinearForm flux(qFlux);

    // ======================================================================
    // Time loop: staggered Robin-Neumann FSI.  Per step:
    //   0D heart model advance -> inlet/outlet pressures; Newmark predictors.
    //   Coupling iteration k = 1 .. couplingIterations:
    //     (a) solid Newmark/NeoHookean SNES solve on the reference mesh;
    //     (b) under-relax + convergence test on the interface displacement;
    //     (c) harmonic ALE extension; mesh velocity; move the fluid mesh;
    //     (d) fluid Oseen solve on the moved mesh.
    //   Commit state, compute the interface flux, update the RCR, write output.
    //
    //   couplingIterations == 1 is loosely coupled; > 1 is strong coupling
    // ======================================================================

    for (size_t step = 1; step <= cfg.nsteps; ++step)
    {
      const auto rep = model.step(dt);
      if (!rep.converged)
      {
        if (isRoot)
          std::cerr << "0D model failed at step " << step << '\n';
        break;
      }

      const auto& s = model.getState();

      if (cfg.subtractHeartHandoffOffset)
      {
        if (step == 1)
          disp0DOffset = s.y;
        disp0D = cfg.heartDisplacementScale * (s.y - disp0DOffset);
      }
      else
      {
        disp0D = cfg.heartDisplacementScale * s.y;
      }

      // Prestress-handoff ramp
      Real loadRamp = 1.0;
      if (cfg.prestressFraction < 1.0 && cfg.prestressRampSteps > 0)
      {
        const Real sRamp = std::min(Real(1.0),
          static_cast<Real>(step - 1) / static_cast<Real>(cfg.prestressRampSteps));
        const Real ss = sRamp * sRamp * (3.0 - 2.0 * sRamp); // smoothstep 0 -> 1
        loadRamp = cfg.prestressFraction + (1.0 - cfg.prestressFraction) * ss;
      }

      pinValue = loadRamp * s.par;
      for (const auto& [tag, bc] : wk)
        outletPressureValue[tag] = loadRamp *
          (s.par - cfg.pressureDropScale * (s.par - bc.pout) - cfg.epicardialDrop);

      // Newmark predictors.
      dPred = dOld;
      auto tmp = solidVelocityOld;
      tmp *= dt;
      dPred += tmp;
      tmp = solidAccelerationOld;
      tmp *= dt * dt * (0.5 - betaN);
      dPred += tmp;

      vPred = solidVelocityOld;
      tmp = solidAccelerationOld;
      tmp *= dt * (1.0 - gammaN);
      vPred += tmp;

      // Initialize the coupling iterate with the predictor.
      dIter = dPred;

      if (isRoot)
      {
        Alert::Info() << "Coronary explicit ALE FSI step " << step << " / " << cfg.nsteps
                      << "  (coupling iterations: " << cfg.couplingIterations << ")"
                      << Alert::Raise;
        Alert::Info() << "  [0D-heart] s.y = " << s.y << " m"
                      << "   imposed disp0D = " << disp0D << " m"
                      << "   (radial vel s.v = " << s.v << " m/s,"
                      << "  penalty k = " << heartK << ")" << Alert::Raise;
      }

      std::map<Attribute, Real> qOut;
      for (const Attribute outlet : BoundaryFluid::Outlets)
        qOut[outlet] = 0.0;
      Real qOutSum = 0.0;
      Real qIn = 0.0;
      bool stepFailed = false;
      Real lastRel = 0.0; // final relative interface change (coupling)
      size_t couplesDone = 0; // coupling iterates actually performed

      // Omega^n and u^n are fixed within the step: assemble massOld once.
      moveMeshWithVertexDisplacement(meshFluid, referenceVertices, uh, aleDispOld);
      massOld.assemble();

      {
        for (size_t couple = 1; couple <= cfg.couplingIterations; ++couple)
        {
          // Solve the solid.  snes.solve() re-evaluates the residual/Jacobian
          // (incl. the lagged fluidStress / robinInterfaceData updated by the
          // previous pass's fluid solve), so pass > 1 is a proper Picard
          // sub-iteration with NO re-assembly needed.
          snes.solve();
          if (!snes.converged())
          {
            if (isRoot)
              std::cerr << "Solid SNES failed at step " << step << " (coupling iterate "
                        << couple << ") after " << snes.getIterationNumber()
                        << " iterations.\n";
            stepFailed = true;
            break;
          }

          // Full (un-relaxed) interface update: dIter <- dState.
          auto delta = dState;
          delta -= dIter;
          PetscReal deltaNorm = 0.0;
          PetscReal stateNorm = 0.0;
          VecNorm(delta.getData(), NORM_2, &deltaNorm);
          VecNorm(dState.getData(), NORM_2, &stateNorm);
          const Real rel = (stateNorm > 0.0)
            ? (static_cast<Real>(deltaNorm) / static_cast<Real>(stateNorm))
            : static_cast<Real>(deltaNorm);
          lastRel = rel;
          couplesDone = couple;
          dIter = dState;

          // Newmark kinematics from the updated iterate dIter (= dState).
          solidAcceleration = dIter;
          solidAcceleration -= dPred;
          solidAcceleration *= 1.0 / (betaN * dt * dt);
          solidVelocity = vPred;
          tmp = solidAcceleration;
          tmp *= gammaN * dt;
          solidVelocity += tmp;

          if (isRoot)
          {
            Alert::Info() << "  coupling iterate " << couple << " / "
                          << cfg.couplingIterations
                          << "  relative interface change = " << rel
                          << "  |d - dPrev| = " << deltaNorm << Alert::Raise;
          }

          // Converged: skip the redundant ALE+fluid solve (iterate 1 must
          // always run the fluid).
          if (couple > 1 && rel < cfg.couplingTolerance)
            break;

          restoreMeshToReference(meshFluid, referenceVertices);
          ale.assemble();
          Solver::KSP(ale).solve();
          aleDisp.setData(dMove.getSolution().getData());

          {
            PetscReal aleNorm = 0.0;
            VecNorm(aleDisp.getData(), NORM_INFINITY, &aleNorm);
            if (!std::isfinite(static_cast<Real>(aleNorm)))
            {
              if (isRoot)
                std::cerr << "ALE lift non-finite at step " << step
                          << " (coupling iterate " << couple
                          << "): the fluid mesh is tangling -- reduce the "
                          << "imposed interface displacement (e.g. "
                          << "-coronary_heart_disp_scale) or the time step.\n";
              stepFailed = true;
              break;
            }
          }

          // ALE mesh velocity.
          meshVelocity = aleDisp;
          meshVelocity -= aleDispOld;
          meshVelocity *= 1.0 / dt;

          moveMeshWithVertexDisplacement(meshFluid, referenceVertices, uh, aleDisp);

          uConv = uOld;
          uConv -= meshVelocity;
          // Only the orthogonal-subscale projections are solved: the taus are
          // pointwise and enter the form directly.  CG + Jacobi mass solves.
          solveMass(vmsL2Conv);
          solveMass(vmsSubProj);
          solveMass(vmsPiTildeProj);

          flow.assemble().setFieldSplits();

          // Inject (rho/dt)(u^n, v)|_{Omega^n} into the velocity block.
          {
            ::Vec b = flow.getLinearSystem().getVector();
            const PetscInt vOff =
              static_cast<PetscInt>(flow.getTestOffsets()[0]); // velocity block
            const ::Vec& mOld = massOld.getVector();
            PetscInt lo = 0, hi = 0;
            VecGetOwnershipRange(mOld, &lo, &hi);
            const PetscScalar* arr = nullptr;
            VecGetArrayRead(mOld, &arr);
            for (PetscInt i = lo; i < hi; ++i)
              if (arr[i - lo] != PetscScalar(0))
                VecSetValue(b, vOff + i, arr[i - lo], ADD_VALUES);
            VecRestoreArrayRead(mOld, &arr);
            VecAssemblyBegin(b);
            VecAssemblyEnd(b);
          }

          Solver::KSP(flow).solve();
          uCur.setData(u.getSolution().getData());
          pCur.setData(p.getSolution().getData());

          fluidTraction.project(Region::Faces, tractionFSI, BoundaryFluid::FSI);
          // Native projections consumed by the NEXT solid solve (lagged
          // Robin-Robin data): the FULL transferred fluid traction + wall
          // velocity.
          tractionTransfer.project(Region::Faces, tractionFSI, BoundaryFluid::FSI);
          uWall.project(Region::Faces, uCur, BoundaryFluid::FSI);
        }
      }

      if (stepFailed)
        break;

      // ------------------------------------------------------------------
      // Wall shear stress -- an OUTPUT diagnostic, computed ONCE per step
      // from the converged fluid state (not once per coupling iterate).
      // The gradient recovery stays an L2 projection: a P1 gradient is
      // elementwise-constant and multivalued at the nodes, and evaluating
      // Jacobian(u) on a boundary face keeps only the TANGENTIAL (surface)
      // derivatives -- exactly the bug the volume recovery fixes -- so this
      // is NOT a candidate for pointwise evaluation (same conclusion as
      // Atrium::computeWallShear and CoupledLV0DCoronary3D).
      // ------------------------------------------------------------------
      {
        solveMass(gradRecProj0);
        gradRec0.setData(gradRecTrial.getSolution().getData());
        solveMass(gradRecProj1);
        gradRec1.setData(gradRecTrial.getSolution().getData());
        solveMass(gradRecProj2);
        gradRec2.setData(gradRecTrial.getSolution().getData());

        PETSc::Variational::TestFunction wssTest(uh);
        const auto onesVec = VectorFunction(dim, [&](const Point&) {
          Math::SpatialVector<Real> o(dim);
          for (Index c = 0; c < static_cast<Index>(dim); ++c)
            o(c) = 1.0;
          return o;
        });
        LinearForm<VelocityFES, ::Vec> wssLoad(wssTest);
        wssLoad = BoundaryIntegral(wallStressRec, wssTest).over(BoundaryFluid::FSI);
        wssLoad.assemble();
        LinearForm<VelocityFES, ::Vec> wssArea(wssTest);
        wssArea = BoundaryIntegral(onesVec, wssTest).over(BoundaryFluid::FSI);
        wssArea.assemble();

        ::Vec bvec = wssLoad.getVector();
        ::Vec mvec = wssArea.getVector();
        ::Vec svec = shearWall.getData();
        // Bound the loop by the MINIMUM local size of the three vectors so
        // a mismatch can never overrun an array (silent heap corruption ->
        // a SEGV many steps later).  All three are on uh and should match;
        // the min is just a hard safety bound.
        PetscInt nb = 0, nm = 0, ns = 0;
        VecGetLocalSize(bvec, &nb);
        VecGetLocalSize(mvec, &nm);
        VecGetLocalSize(svec, &ns);
        const PetscInt n = std::min(ns, std::min(nb, nm));
        const PetscScalar *barr = nullptr, *marr = nullptr;
        PetscScalar* sarr = nullptr;
        VecGetArrayRead(bvec, &barr);
        VecGetArrayRead(mvec, &marr);
        VecGetArray(svec, &sarr);
        for (PetscInt i = 0; i < n; ++i)
          sarr[i] = (std::abs(marr[i]) > 1.0e-30) ? (barr[i] / marr[i]) : PetscScalar(0);
        VecRestoreArray(svec, &sarr);
        VecRestoreArrayRead(mvec, &marr);
        VecRestoreArrayRead(bvec, &barr);
      }

      // Converged interface displacement and consistent kinematics.
      dState = dIter;
      solidAcceleration = dState;
      solidAcceleration -= dPred;
      solidAcceleration *= 1.0 / (betaN * dt * dt);
      solidVelocity = vPred;
      tmp = solidAcceleration;
      tmp *= gammaN * dt;
      solidVelocity += tmp;

      flux = BoundaryIntegral(Dot(uCur, normalFluid), qFlux).over(BoundaryFluid::Inlet);
      flux.assemble();
      qIn = flux(one);
      qOutSum = 0.0;
      for (const Attribute outlet : BoundaryFluid::Outlets)
      {
        flux = BoundaryIntegral(Dot(uCur, normalFluid), qFlux).over(outlet);
        flux.assemble();
        const Real qo = flux(one);
        qOut[outlet] = qo;
        qOutSum += qo;
      }

      // ---- Interface power (Robin-Robin energy consistency) -------------
      // E_Gamma = int_FSI (sigma_f . n_f) . (u_f - d_dot_s) dGamma, the rate
      // of work the fluid does against the kinematic mismatch on Gamma_FSI.
      const auto sigmaFn = muFsi * Mult(strainRateFsi, normalFluid) - pCur * normalFluid;
      flux = BoundaryIntegral(Dot(sigmaFn, uCur), qFlux).over(BoundaryFluid::FSI);
      flux.assemble();
      const Real ePowerFluid = flux(one);
      flux = BoundaryIntegral(Dot(sigmaFn, interfaceSolidVelocity), qFlux)
               .over(BoundaryFluid::FSI);
      flux.assemble();
      const Real ePowerSolid = flux(one);
      const Real eInterface = ePowerFluid - ePowerSolid;

      // ---- Interface slip diagnostic ------------------------------------
      // RMS of |u_f - u_s| on the FSI wall.
      Real slipRms = 0.0;
      {
        const auto slipVec = uCur - interfaceSolidVelocity;
        flux = BoundaryIntegral(Dot(slipVec, slipVec), qFlux).over(BoundaryFluid::FSI);
        flux.assemble();
        const Real slipSq = flux(one);
        flux = BoundaryIntegral(RealFunction(1.0), qFlux).over(BoundaryFluid::FSI);
        flux.assemble();
        const Real fsiArea = flux(one);
        slipRms =
          (fsiArea > 0.0) ? std::sqrt(std::max(Real(0.0), slipSq / fsiArea)) : 0.0;
      }

      // ---- Startup/stability diagnostics --------------------------------
      {
        PetscReal pMin = 0.0, pMax = 0.0;
        VecMin(pCur.getData(), PETSC_NULLPTR, &pMin);
        VecMax(pCur.getData(), PETSC_NULLPTR, &pMax);
        if (isRoot)
          Alert::Info() << "  [diag] p in [" << pMin << ", " << pMax << "] Pa"
                        << "  mass(qIn+qOut) = " << (qIn + qOutSum)
                        << "  E_iface = " << eInterface << " W"
                        << "  slip(RMS|u_f-u_s|) = " << slipRms << " m/s"
                        << "  | coupling: iters = " << couplesDone << "/"
                        << cfg.couplingIterations << "  interface change = " << lastRel
                        << Alert::Raise;
      }

      // R-mu outlet update from the converged interface flux: advances p_tm
      // (implicit Euler + scalar Newton), exports pout = p_c and the
      // rheological factors Phi_a / Phi_v for the next 3D assembly.
      for (const Attribute outlet : BoundaryFluid::Outlets)
        updateOutlet0D(cfg, wrms, model, wk[outlet], qOut[outlet], dt);

      // Commit the step-(n+1) state.
      uOld.setData(uCur.getData());
      pOld.setData(pCur.getData());
      dOld.setData(dState.getData());
      solidVelocityOld.setData(solidVelocity.getData());
      solidAccelerationOld.setData(solidAcceleration.getData());
      // Carry the dynamic VMS subscale u'^{n+1} -> u'^n for the next step.
      vmsSubOld.setData(vmsSub.getSolution().getData());

      aleDispOld.setData(aleDisp.getData());

      xdmf_fluid.write(s.t).flush();
      xdmf_solid.write(s.t).flush();

      if (isRoot && csv)
      {
        csv << s.t << ',' << s.y << ',' << s.v << ',' << s.pv << ',' << s.par << ','
            << s.pd << ',' << qIn << ',' << qOutSum;
        for (const Attribute outlet : BoundaryFluid::Outlets)
          csv << ',' << qOut[outlet] << ',' << wk[outlet].pout;
        // R-mu mechanism diagnostics.  p_tm must stay positive over the
        // whole cycle: it is the model's own check that the Starling throat
        // is doing its job.
        {
          Real pimNow = 0.0, ptmSum = 0.0, muASum = 0.0, muVSum = 0.0;
          Real volSum = 0.0;
          for (const auto& [tag, bc] : wk)
          {
            pimNow = bc.pim;
            ptmSum += bc.ptm;
            muASum += bc.muA;
            muVSum += bc.muV;
            volSum += bc.vol;
          }
          const Real nOut = static_cast<Real>(wk.size());
          csv << ',' << pimNow << ',' << (ptmSum / nOut) << ',' << (muASum / nOut) << ','
              << (muVSum / nOut) << ',' << volSum;
        }
        csv << ',' << eInterface << ',' << lastRel << '\n';
        csv.flush();
      }
    }

    xdmf_fluid.close();
    xdmf_solid.close();
    if (isRoot && csv)
      csv.close();
  }
  catch (const std::exception& e)
  {
    std::cerr << "CoronaryArtery_FSI_Explicit_PETSc_Seq failed: " << e.what() << '\n';
    PetscFinalize();
    return 1;
  }

  PetscFinalize();
  return 0;
}
