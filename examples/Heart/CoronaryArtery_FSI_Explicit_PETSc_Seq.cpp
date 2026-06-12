/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * Explicit (staggered) sequential PETSc coronary ALE FSI.
 *
 * Physics: Carreau-Yasuda blood flow (ALE, conservative BDF1) in a coronary
 * tree, coupled to a hyperelastic arterial wall (total Lagrangian, Newmark)
 * by the Robin-Robin loose coupling of Burman, Durst, Fernandez, Guzman &
 * Ruz (2025), incl. the interface convective stabilization and the
 * alpha = gamma sqrt(rho_s E) scaling.  0D heart model (CCMLC2014) at the
 * inlet; RCR windkessels with implicit outlet impedance at the outlets.
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
 * Attributes: 2 = FSI wall, inlet/outlet caps, 99 = clamped FSI ring band.
 */

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
#include <Rodin/Math/RootFinding/NewtonRaphson.h>
#include <Rodin/Math/RungeKutta/RK4.h>
#include <Rodin/PETSc.h>
#include <Rodin/Solid.h>
#include <Rodin/Solver.h>
#include <Rodin/Variational.h>

#include "Rodin/Heart/CCMLC2014.h"

// Projected-VMS convective stabilization (lagged Oseen), shared with the
// working CoupledLV0DCoronary3D fluid solver.
#include "CoronaryArtery/VMSConvectionIntegrator.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Math;
using namespace Rodin::Solver;
using namespace Rodin::Variational;
using namespace Rodin::Heart;

namespace {
using Model = Rodin::Heart::CCMLC2014T<>;
using MeshType = Geometry::Mesh<Context::Local>;

struct BoundaryFluid {
  static constexpr Attribute FSI = 2;
  static constexpr Attribute Inlet = 4;
  static constexpr std::array<Attribute, 6> Outlets{{7, 8, 9, 10, 14, 15}};
  static constexpr Attribute FSIRing = 99;
};

struct BoundarySolid {
  static constexpr Attribute FSI = 2;
  static constexpr Attribute Inlet = 17;
  static constexpr std::array<Attribute, 6> Outlets{{18, 19, 20, 21, 22, 31}};
  static constexpr Attribute FSIRing = 99;
};

struct RCR {
  // RCR windkessel parameters matched to the WORKING CoupledLV0DCoronary3D
  // setup (Rp, C, Rd, pd0, pc0, pout0) = (5e8, 5e-11, 1e9, 400, 10500, 11000).
  // The previous (C = 1e-12, Rd = 5e10) gave a windkessel with ~no capacitive
  // damping (cap = C/dt ~ 0, so pc jumps instantly with every flux ripple)
  // and a 50x too-large distal resistance: that combination is what makes the
  // explicit 0D-3D outlet coupling oscillate and blow up.
  Real Rp = 5.0e8;
  Real C = 5.0e-11;
  Real Rd = 1.0e9;
  Real pd = 400.0;
  Real pc = 10500.0;
  Real pout = 11000.0;
  Real qd = 0.0;
};

struct CarreauYasuda {
  Real mu0 = 0.0186058;
  Real muInf = 0.0042963;
  Real lambda = 0.2435;
  Real n = 0.2079;
  Real yasuda = 1.5410;
  Real gammaRegularization = 1.0e-3;
};

struct GeoArtery {
  Real Rp;
  Real Lp;
  Real Rd;
  Real Ld;
};

struct OutletFlowLaw {
  /// @brief Proximal surrogate vessel radius.
  Real proximalRadius = 6.e-4;
  /// @brief Proximal surrogate vessel length.
  Real proximalLength = 0.00075;
  /// @brief Distal surrogate vessel radius.
  Real distalRadius = 1e-4;
  /// @brief Distal surrogate vessel length.
  Real distalLength = 0.0025;
  /// @brief Array with radius and large for each branch
  std::unordered_map<Attribute, GeoArtery> geometricParam{
      {7,  {6.e-4,  0.0125, 3e-4, 0.0025}},
      {8,  {4.e-4,  0.01,  2e-4, 0.0025 }},
      {9,  {4.e-4,  0.01,  2e-4, 0.0025 }},
      {10, {4.e-4,  0.01,  2e-4, 0.0025 }},
      {14, {6.e-4,  0.0125, 3e-4, 0.0025}},
      {15, {6.e-4,  0.0125, 3e-4, 0.0025}},
  };
  /// @brief Pressure-drop threshold for the Poiseuille fallback.
  Real pressureDropTolerance = 1.0e-12;
  /// @brief Minimum shear-rate bracket.
  Real minShearRate = 1.0e-8;
  /// @brief Number of RK4 substeps for the WRMS flow integral.
  int integralSteps = 100;
  /// @brief Maximum bracketing expansions for outlet scalar solves.
  int maxBracketIterations = 100;
  /// @brief Wall shear root solver absolute tolerance.
  Real shearAbsoluteTolerance = 1.0e-12;
  /// @brief Wall shear root solver relative tolerance.
  Real shearRelativeTolerance = 1.0e-10;
  /// @brief Wall shear root solver step tolerance.
  Real shearStepTolerance = 1.0e-12;
  /// @brief Wall shear root solver maximum iterations.
  int shearMaxIterations = 50;
  /// @brief Flow inversion root solver absolute tolerance.
  Real flowAbsoluteTolerance = 1.0e-10;
  /// @brief Flow inversion root solver relative tolerance.
  Real flowRelativeTolerance = 1.0e-9;
  /// @brief Flow inversion root solver step tolerance.
  Real flowStepTolerance = 1.0e-12;
  /// @brief Flow inversion root solver maximum iterations.
  int flowMaxIterations = 50;
  /// @brief Flow magnitude treated as zero in pressure-drop inversion.
  Real zeroFlowTolerance = 1.0e-16;
  /// @brief Minimum pressure-drop bracket.
  Real pressureDropBracketMin = 1.0;
  /// @brief Distal capacitor bracket pressure pad.
  Real distalPressureBracketPad = 1000.0;
};

struct Config {
  // All outputs go under resultsDir so the P2/P1 and P1/P1 runs can execute
  // concurrently (different cores) without clobbering each other's files.
  std::string resultsDir = "results_p2p1";
  std::string xdmfBasename = "results_p2p1/CoronaryArtery_FSI_Explicit_P2P1";
  std::string csvPath = "results_p2p1/CoronaryArtery_FSI_Explicit_P2P1.csv";
  Real meshScale = 1.0e-3;

  Real dt = 1.0e-3;
  size_t nsteps = 3 * static_cast<int>(0.85 / 1.0e-3);

  // Fraction of the RCR drop imposed at the 3D outlet faces:
  //   p_out_3D = par - pressureDropScale * (par - p_out_RCR).
  // The RCR capacitor pressure drains toward the venous/myocardial level
  // (pim ~ alpha*pv ~ a few hundred Pa), so imposing the FULL RCR pressure
  // (pressureDropScale = 1) pins the 3D outlets near venous, ~9 kPa BELOW the
  // inlet par -- the entire pressure drop then falls across the epicardial
  // artery, which is wrong (epicardial arteries are low-resistance conduits;
  // the big drop is DOWNSTREAM in the microcirculation, i.e. inside the 0D).
  // A small value references the 3D outlets to par minus a small epicardial
  // drop (~500-1000 Pa), so the artery stays near-isobaric at ~par and the
  // velocities (hence the divergence) drop too.  Tune with
  // -coronary_pressure_drop_scale; 1 = full RCR coupling at the face.
  Real pressureDropScale = 1.0;

  Real fluidDensity = 1060.0;
  Real pressurePenalty = 1.0e-9;
  Real inletBackflowStabilization = 1.0;
  Real outletBackflowStabilization = 1.0;
  CarreauYasuda viscosity;
  // Geometry + nonlinear-solve parameters for the non-Newtonian outlet flow
  // law used by updateRCRNonNew (per-outlet proximal/distal surrogate vessels).
  OutletFlowLaw outletFlowLaw;

  Real inletImpedance = 1.e3;
  Real inletTangentialDamping = 1.e3;

  // Implicit (semi-implicit) outlet resistance scale.  The explicit RCR
  // imposes p_out = p_c + R*Q with the LAGGED flux Q; for the high-resistance
  // small branches (R ~ mu L / r^4 is huge) p_out is hypersensitive to flux
  // ripple and oscillates.  We add a consistency-preserving IMPLICIT local
  // resistance in delta form,
  //    + Z_out int_outlet (u.n)(v.n)  -  Z_out int_outlet (u^n.n)(v.n),
  // with Z_out = outletResistanceScale * 8 mu0 L_p / R_p^2  (the local proximal
  // impedance, per outlet, units Pa s/m -- same as inletImpedance).  Because
  // it is the difference (u - u^n), it VANISHES at steady state, so it does NOT
  // change the converged RCR pressure (the full non-Newtonian p_out stays
  // explicit) -- it only damps the step-to-step flux oscillation implicitly.
  // 0 disables it.
  Real outletResistanceScale = 1.0;

  // Static follower-pressure prestress: ramp 0 -> par(0) in this many
  // increments (0 = disabled); the dynamic loop then starts from the
  // pressurized equilibrium with an exactly balanced first residual.
  size_t prestressSteps = 80;

  // No-prestress startup ramp ("rampage"): when prestressSteps == 0 there is
  // no pressurized equilibrium to hand off, so applying full par on step 1
  // slams the whole system from rest and diverges.  Instead the imposed loads
  // (inlet pin, outlet pout, wall-follower LEVEL) are scaled by a factor that
  // climbs linearly 0 -> 1 over the first rampSteps steps.  Ignored (factor
  // == 1 from step 1) whenever prestressSteps > 0, so the residual-free
  // prestress handoff is never perturbed.  Set to 0 to disable the ramp.
  // Inflating a compliant vessel from rest is a violent transient, so ramp
  // SLOWLY (a few hundred steps).  Prefer the prestress when possible.
  size_t rampSteps = 0;

  // Projected-VMS convective stabilization scale (lagged Oseen, ported from
  // the working CoupledLV0DCoronary3D fluid).  0 disables the VMS terms and
  // recovers the plain Temam-stabilized fluid; 1 is the nominal scaling.
  Real vmsScale = 1.0;

  // Grad-div (continuity / "div-u") stabilization scale.  Adds the consistent
  // term  int tau_C (div u)(div v)  with tau_C = gradDivScale * rho |u| h_K,
  // which improves mass conservation and damps the spurious pressure
  // component / open-boundary pressure layers on the inf-sup-stable P2/P1
  // pair.  Vanishes for divergence-free solutions, so it is consistency-safe.
  // This is NOT a pressure-gradient (PSPG/grad-p) term and is not needed for
  // inf-sup stability (P2/P1 already is stable).  0 disables it.
  // tau_C = gradDivScale * rho * tau_2,  tau_2 = (h/k^2)^2 / (c1 tau_1).
  Real gradDivScale = 1.0;

  Real solidDensity = 1060.0;
  Real solidYoungModulus = 2.0e6;
  Real solidPoissonRatio = 0.4;
  // Dissipative Newmark (gamma > 1/2, beta = (gamma + 1/2)^2 / 4):
  // unconditionally stable and damps the high-frequency content of the
  // startup transient (prestress handoff mismatch absorbed in the first
  // steps) so the acceleration spike cannot blow up the predictor.
  // Classic non-dissipative trapezoidal pair: beta = 0.25, gamma = 0.5.
  // Mildly dissipative Newmark (gamma = 0.6, beta = (gamma + 1/2)^2 / 4).
  // With the level/remainder load split the handoff residual is zero by
  // construction, so only a LIGHT touch of dissipation is kept as a safety
  // margin for the loose coupling; (0.25, 0.5) is worth trying once a clean
  // start is confirmed, or with >= 2 coupling iterations.
  Real newmarkBeta = 0.25;
  Real newmarkGamma = 0.5;

  size_t couplingIterations = 1;   // loose Robin-Robin (no Aitken)
  Real couplingRelaxation = 0.5;   // (unused: Aitken removed)
  Real couplingTolerance = 1.0e-6; // relative interface-displacement tolerance
  // Robin transmission parameter (both sides of Gamma_FSI).  Optimal
  // scaling alpha = gamma*sqrt(rho_s E) (Burman et al. 2025, Sec. 4.1).
  //   robinAlpha <= 0 -> auto from robinGamma; > 0 -> explicit override.
  Real robinAlpha = 0.0;
  // Robin scaling gamma in alpha = gamma sqrt(rho_s E_eq).  Kept at 1
  // (gamma = 2 was tried to lift the Aitken omega off its 0.05 floor but did
  // NOT help -- a stiffer Robin can destabilize the loose scheme).  Tune via
  // -coronary_robin_gamma / -coronary_robin_alpha instead.
  Real robinGamma = 1.0;
};

static Real periodic_activation(Real t) {
  const Real T = 0.85;
  const Real tau = t - T * std::floor(t / T);
  if (tau < 0.13)
    return 0.0;
  if (tau < 0.141)
    return 35.0 * ((tau - 0.13) / 0.011);
  if (tau < 0.281)
    return 35.0;
  if (tau < 0.361)
    return 35.0 - 55.0 * ((tau - 0.281) / 0.08);
  if (tau < 0.45)
    return -20.0;
  return 0.0;
}

static Real load_dependent_relaxation_m0(Real ec) {
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

static Real load_dependent_relaxation_dm0(Real ec) {
  const Real lowEc = 0.0;
  const Real highEc = 2.0;
  const Real lowValue = 1.6;
  const Real highValue = 1.0;
  if (ec <= lowEc || ec >= highEc)
    return 0.0;
  return (highValue - lowValue) / (highEc - lowEc);
}

static Real atrial_pressure(Real t) {
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

  Real alpha = 0.0;
  Real value = minValue;
  if (tau < t1) {
    alpha = -(tau - t1) / t1;
    value = alpha * minValue + (1.0 - alpha) * maxValue;
  } else if (tau < t2) {
    value = maxValue;
  } else if (tau < t3) {
    alpha = -(tau - t3) / (t3 - t2);
    value = alpha * maxValue + (1.0 - alpha) * minValue;
  } else if (tau < t4) {
    alpha = -(tau - t4) / (t4 - t3);
    value = alpha * minValue + (1.0 - alpha) * secondThreshold;
  } else if (tau < t5) {
    value = secondThreshold;
  } else if (tau < t6) {
    alpha = -(tau - t6) / (t6 - t5);
    value = alpha * secondThreshold + (1.0 - alpha) * minValue;
  }
  return value;
}

static Model::Input makeModelInput() {
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
  in.sigma0 = 1.24e5;

  in.Rp = 8.0e6;
  in.Cp = 2.5e-9;
  in.Rd = 1.0e8;
  in.Cd = 1.0e-8;

  in.mu_0 = 0.0526559;
  in.mu_Inf = 0.0052704;
  in.lambda = 0.2435;
  in.n = 0.2079;
  in.m = 0.0035;
  in.yasuda = 1.541;
  in.mu_plasma = 0.0032704;
  in.k_0 = 3.5678;
  in.gamma_c = 10.2754;
  in.k_Inf = 1.5352;
  in.proximalRadius = 0.015;
  in.proximalLength = 0.4;
  in.distalRadius = 0.0007;
  in.distalLength = 0.004;
  in.windkesselRheology =
      Rodin::Heart::CCMLC2014::Model::WindkesselRheology::CarreauYasuda;

  in.Kat = 9.0e-6;
  in.Kp = 5.0e-10;
  in.Kar = 1.3e-5;
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

static void initializeModel(Model &model, const Model::Input &in) {
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
// No partitioner / sharder: this is the sequential build.
static MeshType makeMesh(const Config &cfg, const std::string &meshPath) {
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
static void saveReferenceVertices(const MeshT &mesh,
                                  std::vector<Math::SpatialPoint> &vertices) {
  vertices.resize(mesh.getVertexCount());
  for (auto it = mesh.getVertex(); it; ++it)
    vertices[it->getIndex()] = mesh.getVertexCoordinates(it->getIndex());
}

template <class MeshT>
static void restoreMeshToReference(
    MeshT &mesh, const std::vector<Math::SpatialPoint> &referenceVertices) {
  assert(mesh.getVertexCount() == referenceVertices.size());
  for (auto it = mesh.getVertex(); it; ++it) {
    const Index vertex = it->getIndex();
    mesh.setVertexCoordinates(vertex, referenceVertices[vertex]);
  }
  mesh.flush();
}

template <class MeshT, class FESType, class GridFunctionType>
static void moveMeshWithVertexDisplacement(
    MeshT &mesh, const std::vector<Math::SpatialPoint> &referenceVertices,
    const FESType &displacementFES, const GridFunctionType &displacement) {
  assert(mesh.getVertexCount() == referenceVertices.size());
  const size_t dim = mesh.getSpaceDimension();

  for (auto it = mesh.getVertex(); it; ++it) {
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
struct InterfaceSegment {
  Index fluid;
  Index solid;
};

struct InterfaceMap {
  std::unordered_map<Index, Index> fluidToSolid;
  std::unordered_map<Index, Index> solidToFluid;
  std::vector<InterfaceSegment> segments;
};

static Math::SpatialPoint centroid(const MeshType &mesh,
                                   const Polytope &polytope) {
  Math::SpatialPoint c(mesh.getSpaceDimension());
  c.setZero();

  const auto &vertices = polytope.getVertices();
  for (const auto &v : vertices)
    c += mesh.getVertexCoordinates(v);

  c /= static_cast<Real>(vertices.size());
  return c;
}

// Relabel the one-element-wide band of FSI faces that touch an inlet/outlet cap
// to a distinct "ring" attribute.  Rodin's DirichletBC is codimension-one: it
// constrains whole faces, never a bare edge, so to strongly pin a field to zero
// exactly on the cap rings (the curves where the FSI surface meets the caps) we
// promote the row of FSI faces adjacent to each cap into their own attribute.
// The caller then clamps that attribute and keeps the Robin transmission over
// the remaining (bulk) FSI faces.
static std::size_t tagFSIRingBand(MeshType &mesh, Attribute fsi, Attribute ring,
                                  Attribute inlet,
                                  const std::array<Attribute, 6> &outlets) {
  const std::size_t faceDim = mesh.getDimension() - 1;

  const auto isCap = [&](const Optional<Attribute> &a) {
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
  for (auto it = mesh.getBoundary(); it; ++it) {
    if (!isCap(it->getAttribute()))
      continue;
    for (const auto &v : it->getVertices())
      capVertices.insert(v);
  }

  // Collect (don't mutate during iteration) every FSI face touching a cap vertex.
  std::vector<Index> toRelabel;
  for (auto it = mesh.getBoundary(); it; ++it) {
    if (it->getAttribute() != fsi)
      continue;
    for (const auto &v : it->getVertices())
      if (capVertices.count(v)) {
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
// (coincident faces); a tolerance absorbs floating-point / ASCII round-off.
static InterfaceMap buildInterfaceMap(const MeshType &fluidReferenceMesh,
                                      const MeshType &solidReferenceMesh) {
  InterfaceMap map;

  // Collect local solid FSI face centroids and the interface bounding box.
  std::vector<std::pair<Math::SpatialPoint, Index>> solidFaces;
  Math::SpatialPoint lo, hi;
  bool haveBox = false;
  for (auto it = solidReferenceMesh.getBoundary(); it; ++it) {
    if (it->getAttribute() != BoundarySolid::FSI)
      continue;
    Math::SpatialPoint c = centroid(solidReferenceMesh, *it);
    solidFaces.emplace_back(c, it->getIndex());
    if (!haveBox) {
      lo = c;
      hi = c;
      haveBox = true;
    } else {
      for (Index i = 0; i < static_cast<Index>(c.size()); ++i) {
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
  auto cellOf = [&](const Math::SpatialPoint &c) {
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
  for (auto it = fluidReferenceMesh.getBoundary(); it; ++it) {
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
        for (long long dz = -1; dz <= 1; ++dz) {
          const std::array<long long, 3> key{base[0] + dx, base[1] + dy,
                                             base[2] + dz};
          const auto g = grid.find(key);
          if (g == grid.end())
            continue;
          for (const std::size_t s : g->second) {
            const Real d = Real((solidFaces[s].first - c).norm());
            if (d < best) {
              best = d;
              bestSolid = solidFaces[s].second;
              found = true;
            }
          }
        }

    if (!found || best > tol) {
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

  if (unmatched > 0) {
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

[[maybe_unused]] static Point
forwardFluidPointToSolid(const Point &p, const MeshType &solidReferenceMesh,
                         const InterfaceMap &map) {
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

static Point forwardSolidPointToFluid(const Point &p, const MeshType &fluidMesh,
                                      const InterfaceMap &map) {
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

// True iff every component is finite.  Cross-mesh interface samples on the
// Local mesh can return non-finite values at edge-case quadrature points
// (cell inverse-mapping near face boundaries);
static bool isFiniteVec(const Math::SpatialVector<Real> &x) {
  for (Index i = 0; i < x.size(); ++i)
    if (!std::isfinite(x(i)))
      return false;
  return true;
}


static void setPETScDefault(const char *key, const char *value) {
  PetscBool set = PETSC_FALSE;
  PetscOptionsHasName(PETSC_NULLPTR, PETSC_NULLPTR, key, &set);
  if (!set)
    PetscOptionsSetValue(PETSC_NULLPTR, key, value);
}

static void readOptions(Config &cfg) {
  char meshPath[PETSC_MAX_PATH_LEN] = {};
  PetscBool meshPathSet = PETSC_FALSE;
  PetscOptionsGetString(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_fsi_mesh",
                        meshPath, sizeof(meshPath), &meshPathSet);

  PetscReal dt = cfg.dt;
  PetscBool dtSet = PETSC_FALSE;
  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_dt", &dt,
                      &dtSet);
  if (dtSet)
    cfg.dt = dt;

  PetscInt nsteps = static_cast<PetscInt>(cfg.nsteps);
  PetscBool nstepsSet = PETSC_FALSE;
  PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_nsteps", &nsteps,
                     &nstepsSet);
  if (nstepsSet)
    cfg.nsteps = static_cast<size_t>(std::max<PetscInt>(0, nsteps));


  PetscInt couplingIterations = static_cast<PetscInt>(cfg.couplingIterations);
  PetscBool couplingIterationsSet = PETSC_FALSE;
  PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR,
                     "-coronary_coupling_iterations", &couplingIterations,
                     &couplingIterationsSet);
  if (couplingIterationsSet)
    cfg.couplingIterations =
        static_cast<size_t>(std::max<PetscInt>(1, couplingIterations));

  PetscReal couplingRelaxation = cfg.couplingRelaxation;
  PetscBool couplingRelaxationSet = PETSC_FALSE;
  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
                      "-coronary_coupling_relaxation", &couplingRelaxation,
                      &couplingRelaxationSet);
  if (couplingRelaxationSet)
    cfg.couplingRelaxation = couplingRelaxation;

  PetscReal couplingTolerance = cfg.couplingTolerance;
  PetscBool couplingToleranceSet = PETSC_FALSE;
  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
                      "-coronary_coupling_tolerance", &couplingTolerance,
                      &couplingToleranceSet);
  if (couplingToleranceSet)
    cfg.couplingTolerance = couplingTolerance;

  PetscInt prestressSteps = static_cast<PetscInt>(cfg.prestressSteps);
  PetscBool prestressStepsSet = PETSC_FALSE;
  PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_prestress_steps",
                     &prestressSteps, &prestressStepsSet);
  if (prestressStepsSet)
    cfg.prestressSteps =
        static_cast<size_t>(std::max<PetscInt>(0, prestressSteps));

  PetscInt rampSteps = static_cast<PetscInt>(cfg.rampSteps);
  PetscBool rampStepsSet = PETSC_FALSE;
  PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_ramp_steps",
                     &rampSteps, &rampStepsSet);
  if (rampStepsSet)
    cfg.rampSteps = static_cast<size_t>(std::max<PetscInt>(0, rampSteps));

  PetscReal vmsScale = cfg.vmsScale;
  PetscBool vmsScaleSet = PETSC_FALSE;
  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_vms_scale",
                      &vmsScale, &vmsScaleSet);
  if (vmsScaleSet)
    cfg.vmsScale = std::max<Real>(0.0, vmsScale);

  PetscReal gradDivScale = cfg.gradDivScale;
  PetscBool gradDivScaleSet = PETSC_FALSE;
  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_graddiv_scale",
                      &gradDivScale, &gradDivScaleSet);
  if (gradDivScaleSet)
    cfg.gradDivScale = std::max<Real>(0.0, gradDivScale);

  PetscReal outletResistanceScale = cfg.outletResistanceScale;
  PetscBool outletResistanceScaleSet = PETSC_FALSE;
  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
                      "-coronary_outlet_resistance_scale",
                      &outletResistanceScale, &outletResistanceScaleSet);
  if (outletResistanceScaleSet)
    cfg.outletResistanceScale = std::max<Real>(0.0, outletResistanceScale);

  PetscReal pressureDropScale = cfg.pressureDropScale;
  PetscBool pressureDropScaleSet = PETSC_FALSE;
  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
                      "-coronary_pressure_drop_scale", &pressureDropScale,
                      &pressureDropScaleSet);
  if (pressureDropScaleSet)
    cfg.pressureDropScale = std::max<Real>(0.0, pressureDropScale);

  PetscReal robinAlpha = cfg.robinAlpha;
  PetscBool robinAlphaSet = PETSC_FALSE;
  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_robin_alpha",
                      &robinAlpha, &robinAlphaSet);
  if (robinAlphaSet)
    cfg.robinAlpha = robinAlpha;

  PetscReal robinGamma = cfg.robinGamma;
  PetscBool robinGammaSet = PETSC_FALSE;
  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_robin_gamma",
                      &robinGamma, &robinGammaSet);
  if (robinGammaSet)
    cfg.robinGamma = robinGamma;
}

// Non-Newtonian (Carreau-Yasuda) RCR outlet update.
static void updateRCRNonNew(const Config &cfg, const Attribute &tag,
                            const Model &model, RCR &bc, Real Q, Real dt) {
  const auto &s = model.getState();
  const auto &h = model.getHistory();

  const Real cap = bc.C / dt;
  const Real pcOld = bc.pc;

  const auto &law = cfg.outletFlowLaw;
  const auto &cy = cfg.viscosity;

  const Real radiusP = law.geometricParam.at(tag).Rp;
  const Real lengthP = law.geometricParam.at(tag).Lp;

  const Real radiusD = law.geometricParam.at(tag).Rd;
  const Real lengthD = law.geometricParam.at(tag).Ld;

  auto flowLaw = [&](Real dp, Real L, Real radius) -> std::pair<Real, Real> {
    const Real mu0 = cy.mu0;
    const Real muInf = cy.muInf;
    const Real lambda = cy.lambda;
    const Real n = cy.n;
    const Real yasuda = cy.yasuda;
    const Real delta = mu0 - muInf;

    const Real sgn = (dp >= 0.0) ? 1.0 : -1.0;
    const Real adp = std::abs(dp);

    const Real R0 =
        8.0 * mu0 * L / (std::numbers::pi_v<Real> * std::pow(radius, 4.0));

    if (adp < law.pressureDropTolerance)
      return {dp / R0, 1.0 / R0};

    const Real tauW = radius * adp / (2.0 * L);

    auto mu = [&](Real g) -> Real {
      return muInf + delta * std::pow(1.0 + std::pow(lambda * g, yasuda),
                                      (n - 1.0) / yasuda);
    };

    auto dmu = [&](Real g) -> Real {
      const Real base = 1.0 + std::pow(lambda * g, yasuda);

      return delta * (n - 1.0) * std::pow(base, (n - 1.0 - yasuda) / yasuda) *
             std::pow(lambda, yasuda) * std::pow(g, yasuda - 1.0);
    };

    auto tauMinusTauW = [&](Real g) -> std::pair<Real, Real> {
      const Real m = mu(g);
      const Real dm = dmu(g);
      return {g * m - tauW, m + g * dm};
    };

    Math::RootFinding::NewtonRaphson<Real> rootFinder(
        law.shearAbsoluteTolerance, law.shearRelativeTolerance,
        law.shearStepTolerance, law.shearMaxIterations);

    Real gHi = std::max<Real>(tauW / muInf, law.minShearRate);

    for (int k = 0;
         k < law.maxBracketIterations && tauMinusTauW(gHi).first < 0.0; ++k)
      gHi *= 2.0;

    if (tauMinusTauW(gHi).first < 0.0) {
      std::cerr << "Warning: failed to bracket wall shear rate. "
                << "Using Poiseuille fallback.\n";
      return {dp / R0, 1.0 / R0};
    }

    const auto gammaRoot =
        rootFinder.solve(tauMinusTauW, 0.5 * gHi, law.shearStepTolerance, gHi);

    if (!gammaRoot) {
      std::cerr << "Warning: failed to solve wall shear rate. "
                << "Using Poiseuille fallback.\n";
      return {dp / R0, 1.0 / R0};
    }

    const Real gammaW = *gammaRoot;

    auto integrand = [&](Real g) -> Real {
      if (g <= 0.0)
        return 0.0;

      const Real m = mu(g);
      const Real dm = dmu(g);
      const Real dtau = m + g * dm;

      return std::pow(g, 3.0) * m * m * dtau;
    };

    Math::RungeKutta::RK4 integrator;

    const int steps = law.integralSteps;
    const Real h = gammaW / static_cast<Real>(steps);

    Real I = 0.0;

    auto rhs = [&](Real g, Real y) -> Real {
      (void)y;
      return integrand(g);
    };

    for (int i = 0; i < steps; ++i) {
      const Real g = static_cast<Real>(i) * h;
      integrator.step(I, g, h, I, rhs);
    }

    if (I <= 0.0 || !std::isfinite(I)) {
      std::cerr << "Warning: invalid WRMS integral. "
                << "Using Poiseuille fallback.\n";
      return {dp / R0, 1.0 / R0};
    }

    const Real qAbs = std::numbers::pi_v<Real> * std::pow(radius, 3.0) * I /
                      std::pow(tauW, 3.0);

    const Real dqAbs =
        (std::numbers::pi_v<Real> * std::pow(radius, 3.0) * gammaW -
         3.0 * qAbs) /
        adp;

    if (!std::isfinite(qAbs) || !std::isfinite(dqAbs) || dqAbs <= 0.0) {
      std::cerr << "Warning: invalid WRMS flow derivative. "
                << "Using Poiseuille fallback.\n";
      return {dp / R0, 1.0 / R0};
    }

    return {sgn * qAbs, dqAbs};
  };

  auto solvePressureDropForFlow = [&](Real targetQ, Real L, Real radius,
                                      Real guess) -> Real {
    if (std::abs(targetQ) < law.zeroFlowTolerance)
      return 0.0;

    const Real sgn = (targetQ >= 0.0) ? 1.0 : -1.0;
    const Real qAbs = std::abs(targetQ);

    auto F = [&](Real x) -> std::pair<Real, Real> {
      const auto [q, dq] = flowLaw(sgn * x, L, radius);
      return {sgn * q - qAbs, dq};
    };

    Real hi = std::max<Real>(std::abs(guess), law.pressureDropBracketMin);

    for (int k = 0; k < law.maxBracketIterations && F(hi).first < 0.0; ++k)
      hi *= 2.0;

    if (F(hi).first < 0.0) {
      std::cerr << "Warning: failed to bracket pressure drop for targetQ = "
                << targetQ << ". Returning last upper bound.\n";
      return sgn * hi;
    }

    Math::RootFinding::NewtonRaphson<Real> solver(
        law.flowAbsoluteTolerance, law.flowRelativeTolerance,
        law.flowStepTolerance, law.flowMaxIterations);

    const auto root = solver.solve(F, std::min(std::abs(guess), hi), 0.0, hi);

    if (!root) {
      std::cerr << "Warning: failed to invert flow law for targetQ = "
                << targetQ << ". Returning bracket upper bound.\n";
      return sgn * hi;
    }

    return sgn * (*root);
  };

  const Real alpha = 0.7;
  const Real dPim = alpha * (s.pv - h.nm1.pv) / dt;
  const Real Qim = bc.C * dPim;

  auto distalResidual = [&](Real pc) -> std::pair<Real, Real> {
    const Real pim = alpha * s.pv;
    const Real x = std::max(pc - pim, Real(0.0));
    const auto [qd, dqd] = flowLaw(x, lengthD, radiusD);

    const Real f = cap * (pc - pcOld) - Qim + qd - Q;
    const Real df = cap + (pc > pim ? dqd : Real(0.0));
    return {f, df};
  };

  Math::RootFinding::NewtonRaphson<Real> solver(
      law.flowAbsoluteTolerance, law.flowRelativeTolerance,
      law.flowStepTolerance, law.flowMaxIterations);

  Real span = std::max<Real>(std::abs(Q) / cap + law.distalPressureBracketPad,
                             law.distalPressureBracketPad);

  Real lo = std::min(pcOld, s.pv) - span;
  Real hi = std::max(pcOld, s.pv) + span;

  for (int k = 0; k < law.maxBracketIterations &&
                  distalResidual(lo).first * distalResidual(hi).first > 0.0;
       ++k) {
    span *= 2.0;
    lo = std::min(pcOld, s.pv) - span;
    hi = std::max(pcOld, s.pv) + span;
  }

  if (distalResidual(lo).first * distalResidual(hi).first > 0.0) {
    std::cerr << "Warning: failed to bracket distal capacitor pressure. "
              << "Keeping previous pc.\n";
    bc.pc = pcOld;
  } else {
    const auto pcNew = solver.solve(distalResidual, pcOld, lo, hi);

    if (!pcNew) {
      std::cerr << "Warning: failed to solve distal capacitor equation. "
                << "Keeping previous pc.\n";
      bc.pc = pcOld;
    } else {
      bc.pc = *pcNew;
    }
  }

  const Real pim_f = alpha * s.pv;
  const auto [qd, dqd_f] = flowLaw(std::max(bc.pc - pim_f, Real(0.0)), lengthD, radiusD);
  (void)dqd_f;
  bc.qd = qd;

  const Real oldGuess = bc.pout - bc.pc;
  const Real dpP = solvePressureDropForFlow(Q, lengthP, radiusP, oldGuess);

  bc.pout = bc.pc + dpP;
}


} // namespace

int main(int argc, char **argv) {
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  // Robust direct sub-solvers by default (overridable on the command line).
  setPETScDefault("-ksp_type", "preonly");
  setPETScDefault("-pc_type", "lu");
  setPETScDefault("-pc_factor_mat_solver_type", "mumps");

  // Sequential build: no boost::mpi environment / communicator / context.
  // PETSc runs in serial (PETSC_COMM_SELF) and isRoot is trivially true.
  const bool isRoot = true;

  try {
    Config cfg;
    readOptions(cfg);

    // Create the per-case results directory (idempotent) so XDMF/CSV writes
    // land in their own folder.
    if (!cfg.resultsDir.empty()) {
      std::error_code ec;
      std::filesystem::create_directories(cfg.resultsDir, ec);
      if (ec && isRoot)
        std::cerr << "Warning: could not create results directory '"
                  << cfg.resultsDir << "': " << ec.message() << '\n';
    }

    Model::Input modelInput = makeModelInput();
    Model model(modelInput);
    initializeModel(model, modelInput);

    const std::string fluidMesh =
        "../resources/examples/Heart/CoronaryArtery_FSI_fluid.mesh";
    MeshType meshFluid = makeMesh(cfg, fluidMesh);
    const size_t dimFluid = meshFluid.getSpaceDimension();

    const std::string solidMesh =
        "../resources/examples/Heart/CoronaryArtery_FSI_solid.mesh";
    MeshType meshSolid = makeMesh(cfg, solidMesh);
    const size_t dimSolid = meshSolid.getSpaceDimension();

    // Spatial dimension shared by both blocks (used by the coupling fields).
    const size_t dim = dimFluid;

    const std::size_t fluidRingFaces =
        tagFSIRingBand(meshFluid, BoundaryFluid::FSI, BoundaryFluid::FSIRing,
                       BoundaryFluid::Inlet, BoundaryFluid::Outlets);
    const std::size_t solidRingFaces =
        tagFSIRingBand(meshSolid, BoundarySolid::FSI, BoundarySolid::FSIRing,
                       BoundarySolid::Inlet, BoundarySolid::Outlets);
    if (isRoot)
      std::cout << "FSI ring band: fluid=" << fluidRingFaces
                << " face(s), solid=" << solidRingFaces << " face(s)\n";

    std::vector<Math::SpatialPoint> referenceVertices;
    saveReferenceVertices(meshFluid, referenceVertices);

    const InterfaceMap interfaceMap = buildInterfaceMap(meshFluid, meshSolid);

    // Fluid problem
    using VelocityFES = H1<2, Math::SpatialVector<Real>, MeshType>;
    using PressureFES = H1<1, Real, MeshType>;
    using DisplacementFES = H1<1, Math::SpatialVector<Real>, MeshType>;

    // Solid problem
    using DisplacementFluidFES = H1<1, Math::SpatialVector<Real>, MeshType>;

    VelocityFES uh(std::integral_constant<size_t, 2>{}, meshFluid, dimFluid);
    PressureFES ph(std::integral_constant<size_t, 1>{}, meshFluid);
    DisplacementFluidFES dfh(std::integral_constant<size_t, 1>{}, meshFluid,
                             dimFluid);

    DisplacementFES dh(std::integral_constant<size_t, 1>{}, meshSolid,
                       dimSolid);

    // Fluid trial/test (velocity-pressure).
    PETSc::Variational::TrialFunction u(uh);
    PETSc::Variational::TrialFunction p(ph);
    PETSc::Variational::TestFunction v(uh);
    PETSc::Variational::TestFunction q(ph);


    // Harmonic ALE trial/test functions are created fresh inside
    // solveHarmonicALE() on each call (see below), mirroring the reference.

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

    // ALE problem
    PETSc::Variational::GridFunction aleDisp(
        uh); // full mesh displacement (this step)
    PETSc::Variational::GridFunction aleDispOld(
        uh); // full mesh displacement (previous step)
    PETSc::Variational::GridFunction meshVelocity(
        uh); // (aleDisp - aleDispOld) / dt

    // Solid / mesh-motion state.
    PETSc::Variational::GridFunction dState(
        dh); // total solid displacement (current)
    PETSc::Variational::GridFunction dOld(
        dh); // total solid displacement (previous step)
    PETSc::Variational::GridFunction dIter(dh); // coupling iterate (= dState)
    PETSc::Variational::GridFunction etaState(
        dh); // solid displacement increment (SNES state)
    PETSc::Variational::GridFunction dPred(
        dh); // Newmark displacement predictor
    PETSc::Variational::GridFunction vPred(dh); // Newmark velocity predictor
    PETSc::Variational::GridFunction solidVelocity(dh);
    PETSc::Variational::GridFunction solidAcceleration(dh);
    PETSc::Variational::GridFunction solidVelocityOld(dh);
    PETSc::Variational::GridFunction solidAccelerationOld(dh);

    PETSc::Variational::GridFunction fluidTraction(dfh);
    // Interface transfer fields: projected natively on the fluid FSI
    // faces after each fluid solve; the solid samples only P1 VALUES
    // cross-mesh (gradient evaluation at hand-built face Points is
    // unreliable on the Local mesh).
    PETSc::Variational::GridFunction tractionTransfer(dfh);
    PETSc::Variational::GridFunction uWall(dfh);

    auto zero = VectorFunction(dim, [&](const Point &) {
      Math::SpatialVector<Real> value(dim);
      value.setZero();
      return value;
    });

    uOld = zero;
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
    xdmf_fluid.write(0.0).flush();

    IO::XDMF xdmf_solid(cfg.xdmfBasename + "solid");
    xdmf_solid.setMesh(meshSolid);
    xdmf_solid.add("Displacement", dState);
    xdmf_solid.add("SolidVelocity", solidVelocity);
    xdmf_solid.write(0.0).flush();

    std::ofstream csv;
    if (isRoot) {
      csv.open(cfg.csvPath);
      csv << "t,lv_y,lv_v,lv_pv,lv_par,lv_pd,q_in,q_out_total";
      for (const Attribute outlet : BoundaryFluid::Outlets)
        csv << ",q_out_" << outlet << ",p_out_" << outlet;
      csv << ",e_interface,coupling_rel\n";
    }

    std::map<Attribute, RCR> wk;
    for (const Attribute outlet : BoundaryFluid::Outlets)
      wk.emplace(outlet, RCR{});

    Real pinValue = model.getState().par;
    std::map<Attribute, Real> outletPressureValue;
    for (const Attribute outlet : BoundaryFluid::Outlets)
      outletPressureValue[outlet] = wk.at(outlet).pout;

    auto pin = RealFunction([&](const Point &) { return pinValue; });
    auto pout0 = RealFunction([&](const Point &) {
      return outletPressureValue[BoundaryFluid::Outlets[0]];
    });
    auto pout1 = RealFunction([&](const Point &) {
      return outletPressureValue[BoundaryFluid::Outlets[1]];
    });
    auto pout2 = RealFunction([&](const Point &) {
      return outletPressureValue[BoundaryFluid::Outlets[2]];
    });
    auto pout3 = RealFunction([&](const Point &) {
      return outletPressureValue[BoundaryFluid::Outlets[3]];
    });
    auto pout4 = RealFunction([&](const Point &) {
      return outletPressureValue[BoundaryFluid::Outlets[4]];
    });
    auto pout5 = RealFunction([&](const Point &) {
      return outletPressureValue[BoundaryFluid::Outlets[5]];
    });

    // ----------------------------------------------------------------------
    // Constitutive and time-integration constants.
    // ----------------------------------------------------------------------
    const Real solidLambda =
        cfg.solidYoungModulus * cfg.solidPoissonRatio /
        ((1.0 + cfg.solidPoissonRatio) * (1.0 - 2.0 * cfg.solidPoissonRatio));
    const Real solidMu =
        cfg.solidYoungModulus / (2.0 * (1.0 + cfg.solidPoissonRatio));

    const Real dt = cfg.dt;
    const Real betaN = cfg.newmarkBeta;
    const Real gammaN = cfg.newmarkGamma;
    const Real solidMass = cfg.solidDensity / (betaN * dt * dt);
    const Real solidVelocityCoeff = gammaN / (betaN * dt);

    const Real yeohC1 = 70000.0;
    const Real yeohC2 = 200000.0;
    const Real yeohC3 = 19000000.0;
    const Real yeohKappa = 1400000.0;

    // Equivalent small-strain isotropic moduli of the Yeoh law:
    //   shear  mu_eq = 2 c1,
    //   Young  E_eq  = 9 kappa mu_eq / (3 kappa + mu_eq).
    // E_eq is the representative stiffness used by the optimal Robin scaling
    // alpha = gamma sqrt(rho_s E)  (Burman, Durst, Fernandez, Guzman & Ruz
    // 2025, Sec. 4.1).  For (c1, kappa) = (7e4, 1.4e6): mu_eq = 1.4e5,
    // E_eq ~ 4.06e5 Pa, so alpha ~ gamma sqrt(1060 * 4.06e5) ~ 2.08e4.
    const Real solidShearEquiv = 2.0 * yeohC1;
    const Real solidYoungEquiv =
        9.0 * yeohKappa * solidShearEquiv / (3.0 * yeohKappa + solidShearEquiv);

    // alpha = gamma * sqrt(rho_s * E_eq) unless explicitly overridden (>0).
    const Real robinAlpha =
        (cfg.robinAlpha > 0.0)
            ? cfg.robinAlpha
            : cfg.robinGamma *
                  std::sqrt(cfg.solidDensity * solidYoungEquiv);
    if (isRoot)
      Alert::Info() << "Robin parameter alpha = " << robinAlpha
                    << "  (gamma * sqrt(rho_s E_eq), gamma = " << cfg.robinGamma
                    << ", E_eq = " << solidYoungEquiv << " Pa)" << Alert::Raise;
    const Real robinVelocityCoeff = robinAlpha * solidVelocityCoeff;


    const auto &cy = cfg.viscosity;
    const Real gammaReg = cy.gammaRegularization;
    const Real deltaMu = cy.mu0 - cy.muInf;

    const auto normalFluid = BoundaryNormal(meshFluid);

    Solid::Yeoh law(yeohC1, yeohC2, yeohC3, yeohKappa);
    Solid::MaterialTangent solidTangent(law, d, w, dState);
    Solid::InternalForce solidInternal(law, w, dState);

    // Current fluid solution (read back after each Oseen solve).
    PETSc::Variational::GridFunction uCur(uh);
    PETSc::Variational::GridFunction pCur(ph);
    uCur = zero;
    pCur = 0.0;

    // Fluid Cauchy traction on Gamma_FSI (current config):
    //   tractionFSI = p n - mu(grad u + grad u^T) n = -sigma_f n_f,
    // i.e. the traction the fluid exerts on the wall (p>0 pushes outward).
    const auto gradUfsi = Jacobian(uCur);
    const auto strainRateFsi = gradUfsi + Transpose(gradUfsi);
    const auto symUfsi = 0.5 * strainRateFsi;
    const auto shearFsi =
        Sqrt(gammaReg * gammaReg + 2.0 * Dot(symUfsi, symUfsi));
    const auto muFsi =
        cy.muInf + deltaMu * Pow(1.0 + Pow(cy.lambda * shearFsi, cy.yasuda),
                                 (cy.n - 1.0) / cy.yasuda);

    // PHYSICAL (unscaled) fluid traction t_f = -sigma_f n_f: the fluid's OWN
    // lagged Robin datum sigma_f^{lag} n and the traction output field.  The
    // traction-scale knobs must NOT appear here: scaling the fluid-side datum
    // moves the fluid fixed point to (1 - scale) sigma_f n = alpha (u_s - u),
    // i.e. it violates kinematic continuity whenever scale != 1.
    const auto tractionFSI =
        (1.0 * pCur) * normalFluid -
        (1.0 * muFsi) * Mult(strainRateFsi, normalFluid);


    // Wall load = the FULL fluid traction transferred across the interface,
    // nothing else.  The solid feels only sigma_f n_f from the fluid (pressure
    // + viscous); there is NO separately-imposed follower-pressure level.  The
    // traction is projected onto the FSI faces (tractionTransfer) and pulled
    // back to the reference solid surface by J_a in fluidStress below.

    // Areal stretch J_a = A_t/A_0 at the CURRENT iterate dIter: pulls
    // per-current-area interface data back to the reference solid surface.
    auto arealStretchAt = [&](const Point &xs) -> Real {
      Real stretch = 1.0;
      const auto &verts = xs.getPolytope().getVertices();
      if (dim == 3 && verts.size() == 3) {
        Math::SpatialPoint X0 = meshSolid.getVertexCoordinates(verts[0]);
        Math::SpatialPoint X1 = meshSolid.getVertexCoordinates(verts[1]);
        Math::SpatialPoint X2 = meshSolid.getVertexCoordinates(verts[2]);

        Math::SpatialPoint x0 = X0, x1 = X1, x2 = X2;
        for (Index c = 0; c < 3; ++c) {
          x0(c) += dIter[dh.getGlobalIndex({0, verts[0]}, c)];
          x1(c) += dIter[dh.getGlobalIndex({0, verts[1]}, c)];
          x2(c) += dIter[dh.getGlobalIndex({0, verts[2]}, c)];
        }

        const auto triArea = [](const Math::SpatialPoint &a,
                                const Math::SpatialPoint &b,
                                const Math::SpatialPoint &c) -> Real {
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
    const auto arealStretch = RealFunction(arealStretchAt);

    auto fluidStress = VectorFunction(dim, [&](const Point &xs) {
      const Point xf = forwardSolidPointToFluid(xs, meshFluid, interfaceMap);

      Math::SpatialVector<Real> value(dim);
      value.setZero();
      const auto force = tractionTransfer(xf);
      for (Index i = 0; i < static_cast<Index>(dim); ++i)
        value(i) = force(i);

      if (!isFiniteVec(value)) {
        static bool reported = false;
        if (!reported) {
          reported = true;
          const auto &pc = xs.getPhysicalCoordinates();
          Alert::Warning()
              << "fluidStress: non-finite cross-mesh traction sample at xs=("
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
    // Solid interface velocity u_s evaluated at a FLUID interface point: this
    // is the fluid Robin data alpha*u_s.  Solid is solved first in the
    // staggered order, so we map the fluid point to its matching solid face and
    // sample the just-computed solid interface velocity 'solidVelocity' there
    // (cross-mesh).
    auto interfaceSolidVelocity = VectorFunction(dim, [&](const Point &xf) {
      const Point xs = forwardFluidPointToSolid(xf, meshSolid, interfaceMap);

      Math::SpatialVector<Real> value(dim);
      value.setZero();

      const auto us = solidVelocity(xs);
      for (Index i = 0; i < static_cast<Index>(dim); ++i)
        value(i) = us(i);

      if (!isFiniteVec(value)) {
        static bool reported = false;
        if (!reported) {
          reported = true;
          Alert::Warning()
              << "interfaceSolidVelocity: non-finite cross-mesh sample; "
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
    const auto muLag =
        cy.muInf + deltaMu * Pow(1.0 + Pow(cy.lambda * shearLag, cy.yasuda),
                                 (cy.n - 1.0) / cy.yasuda);

    const auto outletBeta = Max(-Dot(transportLag, normalFluid), 0.0);
    const auto inletBeta = Max(Dot(transportLag, normalFluid), 0.0);
    const auto outletBackflow =
        0.5 * cfg.outletBackflowStabilization * cfg.fluidDensity * outletBeta;
    const auto inletBackflow =
        0.5 * cfg.inletBackflowStabilization * cfg.fluidDensity * inletBeta;

    // ----------------------------------------------------------------------
    // Projected-VMS convective stabilization (lagged Oseen), ported from the
    // working CoupledLV0DCoronary3D fluid.  Adds streamline-oriented
    // dissipation along the ALE transport velocity uConv = u^n - w without
    // isotropic viscosity:
    //
    //   + int_K tau rho^2 ((grad u) uConv) . ((grad v) uConv)
    //   - int_K rho (tau rho Pi[(grad uConv) uConv] + u'/dt) . ((grad v) uConv)
    //
    // with the dynamic subscale
    //   u' = tau rho (1/dt u'^n - ((grad uConv)uConv - Pi[(grad uConv)uConv])).
    //
    // tau folds cfg.vmsScale, so vmsScale == 0 zeroes tau and the whole VMS
    // contribution (the cheap projection solves still run, harmlessly).  The
    // frozen velocity is the SAME ALE transport velocity used by the Oseen
    // convection term above, so the stabilization is consistent on the moving
    // mesh.
    // ----------------------------------------------------------------------
    using namespace Rodin::Examples::Heart;

    // ALE convecting velocity u^n - w as a grid function (refreshed each
    // coupling iterate, since the mesh velocity w changes with the iterate).
    PETSc::Variational::GridFunction uConv(uh);
    uConv = uOld;

    // P1 scalar tau and P2 vector projection/subscale fields (P2 = velocity
    // space uh, matching CoupledLV0DCoronary3D's VMSFES).
    PressureFES tauFes(std::integral_constant<size_t, 1>{}, meshFluid);
    PETSc::Variational::TrialFunction vmsTau(tauFes);
    PETSc::Variational::TestFunction vmsTauTest(tauFes);

    // Grad-div (continuity) parameter tau_C = rho*tau2, projected onto P1.
    PETSc::Variational::TrialFunction vmsTauC(tauFes);
    PETSc::Variational::TestFunction vmsTauCTest(tauFes);

    PETSc::Variational::TrialFunction vmsUp(uh);    // Pi[(grad uConv) uConv]
    PETSc::Variational::TestFunction vmsVp(uh);
    PETSc::Variational::TrialFunction vmsSub(uh);    // dynamic subscale u'^{n+1}
    PETSc::Variational::GridFunction vmsSubOld(uh);  // subscale history u'^n
    vmsSubOld = zero;

    // Orthogonal grad-div (p~) IMEX fields (Codina/Principe stable split):
    //   tauC      = tau_C^{n+1} = rho tau2          (P1; bilinear coeff, LHS)
    //   sqrtTauC  = sqrt(tau_C^{n+1})               (P1; linear coeff, LHS tau)
    //   sqrtTauCOld = sqrt(tau_C^n)                 (P1; lagged, for pi~)
    //   piTilde   = Pi( sqrt(tau_C^n) div u^n )     (P1; standard L2, RHS)
    // The stabilization assembled is
    //   + int tau_C^{n+1} (div u^{n+1})(div v)
    //   - int sqrt(tau_C^{n+1}) piTilde (div v),
    // i.e. tau^{n+1} implicit (LHS), the projection lagged to time n (RHS),
    // with the symmetric sqrt(tau) weighting that gives the telescopic energy
    // bound (temporal stability).  tau is a coefficient, never an unknown.
    PETSc::Variational::TrialFunction vmsSqrtTauC(tauFes);
    PETSc::Variational::GridFunction vmsSqrtTauCOld(tauFes);
    vmsSqrtTauCOld = 0.0;
    PETSc::Variational::TrialFunction vmsPiTilde(tauFes);

    // Frozen convective acceleration (grad uConv) uConv.
    const auto vmsConvectionTarget = Mult(Jacobian(uConv), uConv);

    // Shared stabilization parameter tau1 (Codina), with c1 = 4, c2 = 2,
    // k = 2 (P2), h_K = |K|^{1/d}, nu = mu0/rho:
    //   tau1 = ( c1 k^4 nu / h^2 + c2 k |u| / h )^-1.
    // A representative low-shear viscosity mu0 is used (conservative).
    auto tau1At = [&](const Point &pp) -> Real {
      const auto uc = uConv.getValue(pp);
      const Real nu = cy.mu0 / cfg.fluidDensity;
      const Real hK = std::pow(pp.getPolytope().getMeasure(),
                               1.0 / pp.getPolytope().getDimension());
      const Real k = 2.0;
      const Real speed = std::sqrt(Math::dot(uc, uc));
      return 1.0 / (4.0 * std::pow(k, 4.0) * nu / (hK * hK) +
                    2.0 * k * speed / hK);
    };

    // Convective subscale parameter tau_K = vmsScale/(rho/dt + rho/tau1)
    //   = vmsScale (1/rho)(1/dt + 1/tau1)^-1.
    RealFunction vmsTauFn = [&, tau1At](const Point &pp) -> Real {
      const Real tau1 = tau1At(pp);
      return cfg.vmsScale *
             (1.0 / (cfg.fluidDensity / dt + cfg.fluidDensity / tau1));
    };

    // Dynamic subscale update (L2-projected into vmsSub).
    auto vmsSubUpdate = VectorFunction(
        dim, [&](const Point &pp) -> Math::SpatialVector<Real> {
      const auto conv = vmsConvectionTarget.getValue(pp);
      const auto proj = vmsUp.getSolution().getValue(pp);
      const auto old = vmsSubOld.getValue(pp);
      const Real tau = vmsTau.getSolution().getValue(pp);

      Math::SpatialVector<Real> out(dim);
      for (Index c = 0; c < static_cast<Index>(dim); ++c)
        out(c) = tau * cfg.fluidDensity *
                 (1.0 / dt * old(c) - (conv(c) - proj(c)));
      return out;
    });

    // L2 projections (mass matrices reassembled each iterate as uConv moves).
    Problem vmsL2Conv(vmsUp, vmsVp);
    vmsL2Conv =
        Integral(vmsUp, vmsVp) - Integral(vmsConvectionTarget, vmsVp);

    Problem vmsTauProj(vmsTau, vmsTauTest);
    vmsTauProj = Integral(vmsTau, vmsTauTest) - Integral(vmsTauFn, vmsTauTest);

    Problem vmsSubProj(vmsSub, vmsVp);
    vmsSubProj = Integral(vmsSub, vmsVp) - Integral(vmsSubUpdate, vmsVp);

    // Grad-div / pressure-subscale parameter tau_C = gradDivScale * rho * tau2,
    //   tau2 = (h/k^2)^2 / (c1 tau1),  c1 = 4, k = 2.
    // The pressure subscale is p~ = -tau_C P^perp(div u).  gradDivScale folds
    // into tau_C, so gradDivScale == 0 zeroes the div-u stabilization.
    // Plain helper (callable lambda) so both the tau_C and sqrt(tau_C)
    // RealFunctions evaluate the SAME expression without invoking a
    // RealFunction as a callable.
    auto rhoTau2At = [&, tau1At](const Point &pp) -> Real {
      const Real tau1 = tau1At(pp);
      const Real hK = std::pow(pp.getPolytope().getMeasure(),
                               1.0 / pp.getPolytope().getDimension());
      const Real k = 2.0, c1 = 4.0;
      const Real lref = hK / (k * k); // h / k^2
      const Real tau2 = (lref * lref) / (c1 * tau1);
      return cfg.gradDivScale * cfg.fluidDensity * tau2;
    };
    // sqrt(tau_C) = sqrt(gradDivScale rho tau2) >= 0 is the ONLY tau field we
    // project.  The bilinear coefficient tau_C^{n+1} is then formed as the
    // POINTWISE SQUARE of this projected sqrt-field (vmsTauC = vmsSqrtTauC.^2,
    // computed in the loop via VecPointwiseMult).  This guarantees the
    // discrete identity tau_C = (sqrt(tau_C))^2, so the implicit bilinear
    //   int (sqrt tau)^2 (div u)(div v)
    // and the lagged linear  int (sqrt tau) pi~ (div v)  use a CONSISTENT
    // sqrt(tau) -- which is exactly what the IMEX telescopic energy bound
    // requires (projecting tau and sqrt(tau) separately would NOT satisfy
    // tau = (sqrt tau)^2 and would break the bound).
    RealFunction sqrtRhoTau2Fn = [rhoTau2At](const Point &pp) -> Real {
      return std::sqrt(std::max<Real>(0.0, rhoTau2At(pp)));
    };

    // sqrt(tau_C^{n+1}) (current).
    Problem vmsSqrtTauCProj(vmsSqrtTauC, vmsTauCTest);
    vmsSqrtTauCProj = Integral(vmsSqrtTauC, vmsTauCTest) -
                      Integral(sqrtRhoTau2Fn, vmsTauCTest);

    // pi~^n = Pi( sqrt(tau_C^n) div u^n ): STANDARD L2 projection of the
    // sqrt(tau)-weighted divergence of the PREVIOUS-step velocity (uOld) with
    // the PREVIOUS-step sqrt(tau) (vmsSqrtTauCOld).  Fully lagged (time n).
    Problem vmsPiTildeProj(vmsPiTilde, vmsTauCTest);
    vmsPiTildeProj =
        Integral(vmsPiTilde, vmsTauCTest) -
        Integral(vmsSqrtTauCOld * Div(uOld), vmsTauCTest);

    // Per-outlet implicit local resistance coefficient
    //   Z_out = outletResistanceScale * 8 mu0 L_p / R_p^2   [Pa s/m],
    // the local proximal impedance (R_p, L_p = proximal surrogate radius and
    // length per outlet; mu0 = low-shear viscosity -> the LARGEST, hence the
    // most conservatively-damping, Newtonian estimate).  Used in delta form so
    // the non-Newtonian RCR pressure is untouched at convergence.
    std::array<Real, 6> outletZ{};
    for (size_t i = 0; i < BoundaryFluid::Outlets.size(); ++i) {
      const auto &g =
          cfg.outletFlowLaw.geometricParam.at(BoundaryFluid::Outlets[i]);
      outletZ[i] = cfg.outletResistanceScale * 8.0 * cfg.viscosity.mu0 * g.Lp /
                   (g.Rp * g.Rp);
    }

    // BDF1 mass split: (rho/dt)[(u,v)_{n+1} - (u^n,v)_n]; the implicit part
    // lives in 'flow', the explicit u^n part is 'massOld', assembled on the
    // PREVIOUS configuration and injected into the RHS at solve time.
    Problem flow(u, p, v, q);
    flow =
        (cfg.fluidDensity / dt) * Integral(u, v) +
        cfg.fluidDensity * Integral(Dot(convU, v)) +
        0.5 * cfg.fluidDensity * Integral(divGeomTemam * Dot(u, v)) +
        // Projected-VMS convective stabilization (bilinear + subtracted
        // linear); tau folds cfg.vmsScale so both vanish when vmsScale == 0.
        VMSConvectionBilinearIntegrator(u, v, uConv, vmsTau.getSolution(),
                                        cfg.fluidDensity) -
        VMSConvectionLinearIntegrator(v, vmsSub.getSolution(), uConv,
                                      vmsUp.getSolution(), vmsTau.getSolution(),
                                      cfg.fluidDensity, dt) +
        // Orthogonal grad-div / pressure subscale p~, STABLE IMEX split:
        //   + int tau_C^{n+1} (div u^{n+1})(div v)          [implicit, LHS]
        //   - int sqrt(tau_C^{n+1}) pi~^n (div v)           [lagged, RHS]
        // with pi~^n = Pi( sqrt(tau_C^n) div u^n ).  tau_C folds
        // cfg.gradDivScale so both vanish when gradDivScale == 0.
        VMSGradDivBilinearIntegrator(u, v, vmsTauC.getSolution()) -
        VMSGradDivLinearIntegrator(v, vmsPiTilde.getSolution(),
                                   vmsSqrtTauC.getSolution()) +
        2.0 * Integral(muLag * symU, symV) - Integral(p, Div(v)) +
        Integral(Div(u), q) + cfg.pressurePenalty * Integral(p, q) +
        BoundaryIntegral(inletBackflow * Dot(u, v)).over(BoundaryFluid::Inlet) +
        BoundaryIntegral(outletBackflow * Dot(u, v))
            .over(BoundaryFluid::Outlets[0], BoundaryFluid::Outlets[1],
                  BoundaryFluid::Outlets[2], BoundaryFluid::Outlets[3],
                  BoundaryFluid::Outlets[4], BoundaryFluid::Outlets[5]) +
        BoundaryIntegral(pin * Dot(v, normalFluid)).over(BoundaryFluid::Inlet) +
        BoundaryIntegral(pout0 * Dot(v, normalFluid))
            .over(BoundaryFluid::Outlets[0]) +
        BoundaryIntegral(pout1 * Dot(v, normalFluid))
            .over(BoundaryFluid::Outlets[1]) +
        BoundaryIntegral(pout2 * Dot(v, normalFluid))
            .over(BoundaryFluid::Outlets[2]) +
        BoundaryIntegral(pout3 * Dot(v, normalFluid))
            .over(BoundaryFluid::Outlets[3]) +
        BoundaryIntegral(pout4 * Dot(v, normalFluid))
            .over(BoundaryFluid::Outlets[4]) +
        BoundaryIntegral(pout5 * Dot(v, normalFluid))
            .over(BoundaryFluid::Outlets[5])
        // Implicit local outlet resistance (delta form, per outlet):
        //   + Z_out (u.n)(v.n)  -  Z_out (u^n.n)(v.n).
        // The implicit (u) part damps a step-to-step flux change; the lagged
        // (u^n) part cancels it at steady state (u = u^n), so the converged
        // RCR outlet pressure is unchanged -- only the oscillation is killed.
        + outletZ[0] * BoundaryIntegral(Dot(Dot(u, normalFluid) * normalFluid, v)).over(BoundaryFluid::Outlets[0])
        - outletZ[0] * BoundaryIntegral(Dot(Dot(uOld, normalFluid) * normalFluid, v)).over(BoundaryFluid::Outlets[0])
        + outletZ[1] * BoundaryIntegral(Dot(Dot(u, normalFluid) * normalFluid, v)).over(BoundaryFluid::Outlets[1])
        - outletZ[1] * BoundaryIntegral(Dot(Dot(uOld, normalFluid) * normalFluid, v)).over(BoundaryFluid::Outlets[1])
        + outletZ[2] * BoundaryIntegral(Dot(Dot(u, normalFluid) * normalFluid, v)).over(BoundaryFluid::Outlets[2])
        - outletZ[2] * BoundaryIntegral(Dot(Dot(uOld, normalFluid) * normalFluid, v)).over(BoundaryFluid::Outlets[2])
        + outletZ[3] * BoundaryIntegral(Dot(Dot(u, normalFluid) * normalFluid, v)).over(BoundaryFluid::Outlets[3])
        - outletZ[3] * BoundaryIntegral(Dot(Dot(uOld, normalFluid) * normalFluid, v)).over(BoundaryFluid::Outlets[3])
        + outletZ[4] * BoundaryIntegral(Dot(Dot(u, normalFluid) * normalFluid, v)).over(BoundaryFluid::Outlets[4])
        - outletZ[4] * BoundaryIntegral(Dot(Dot(uOld, normalFluid) * normalFluid, v)).over(BoundaryFluid::Outlets[4])
        + outletZ[5] * BoundaryIntegral(Dot(Dot(u, normalFluid) * normalFluid, v)).over(BoundaryFluid::Outlets[5])
        - outletZ[5] * BoundaryIntegral(Dot(Dot(uOld, normalFluid) * normalFluid, v)).over(BoundaryFluid::Outlets[5])
        + cfg.inletImpedance *
              BoundaryIntegral(Dot(Dot(u, normalFluid) * normalFluid, v))
                  .over(BoundaryFluid::Inlet)
        + cfg.inletTangentialDamping *
              BoundaryIntegral(Dot(duTangential, v)).over(BoundaryFluid::Inlet)
        // Robin-Robin transmission (fluid side), Burman et al. (2025):
        //   sigma_f n + alpha u = alpha d_dot_s + lambda^{k-1},
        // with lambda^{k-1} = tractionFSI at the PREVIOUS correction (= the
        // previous time step at the first correction, i.e. lambda^{n-1} in
        // the loose kappa = 0 scheme).  Enters with PLUS (a minus sign would
        // impose a traction-free wall and diverge).
        + robinAlpha * BoundaryIntegral(u,v).over(BoundaryFluid::FSI)
        - robinAlpha * BoundaryIntegral(interfaceSolidVelocity,v).over(BoundaryFluid::FSI)
        + BoundaryIntegral(tractionFSI,v).over(BoundaryFluid::FSI)
        // Interface convective stabilization (Burman et al. 2025, eq. 13):
        // -(rho/2) (transportLag.n)(u.v) on Sigma; controls the convective
        // energy on the moving wall (omitting it -> added-mass growth).
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
    auto interfaceSolidDisplacement = VectorFunction(dim, [&](const Point& xf)
    {
      const Point xs =
        forwardFluidPointToSolid(xf, meshSolid, interfaceMap);

      Math::SpatialVector<Real> value(dim);
      value.setZero();

      const auto ds = dIter(xs);
      for (Index i = 0; i < static_cast<Index>(dim); ++i)
        value(i) = ds(i);

      if (!isFiniteVec(value)) {
        static bool reported = false;
        if (!reported) {
          reported = true;
          Alert::Warning()
              << "interfaceSolidDisplacement: non-finite cross-mesh sample; "
              << "using zero there (warned once)." << Alert::Raise;
        }
        value.setZero();
      }

      return value;
    });

    // Harmonic ALE lift (reference config).  Must be RE-ASSEMBLED on every
    // coupling iterate: assemble() bakes the FSI Dirichlet values (dIter)
    // into the system; KSP::solve() alone does not re-evaluate BCs.
      PETSc::Variational::TrialFunction dMove(uh);
      PETSc::Variational::TestFunction vMove(uh);

      Problem ale(dMove, vMove);
      ale = Integral(Jacobian(dMove), Jacobian(vMove)) +
            DirichletBC(dMove, interfaceSolidDisplacement).on(BoundaryFluid::FSI) +
            DirichletBC(dMove, zero)
                .on(BoundaryFluid::Inlet, BoundaryFluid::Outlets[0],
                    BoundaryFluid::Outlets[1], BoundaryFluid::Outlets[2],
                    BoundaryFluid::Outlets[3], BoundaryFluid::Outlets[4],
                    BoundaryFluid::Outlets[5], BoundaryFluid::FSIRing);

    // Solid Robin data, per reference area:
    //   J_a [ rVC (dState - dPred) + alpha vPred - alpha u_f^{k-1} ],
    // with u_f sampled from the projected uWall (latest fluid iterate).
    auto robinInterfaceData = VectorFunction(dim, [&](const Point &xs) {
      const Point xf = forwardSolidPointToFluid(xs, meshFluid, interfaceMap);
      const Real Ja = arealStretchAt(xs);

      Math::SpatialVector<Real> value(dim);
      value.setZero();

      // Guarded cross-mesh fluid velocity sample (see fluidStress note).
      Math::SpatialVector<Real> uf(dim);
      uf.setZero();
      {
        const auto ufRaw = uWall(xf);
        for (Index i = 0; i < static_cast<Index>(dim); ++i)
          uf(i) = ufRaw(i);
        if (!isFiniteVec(uf)) {
          static bool reported = false;
          if (!reported) {
            reported = true;
            const auto &pc = xs.getPhysicalCoordinates();
            Alert::Warning()
                << "robinInterfaceData: non-finite cross-mesh velocity sample "
                << "at xs=(" << pc.transpose()
                << "); using zero there (warned once)." << Alert::Raise;
          }
          uf.setZero();
        }
      }
      const auto vp = vPred(xs);
      const auto dS = dState(xs);
      const auto dP = dPred(xs);
      for (Index i = 0; i < static_cast<Index>(dim); ++i)
        value(i) = Ja * (robinVelocityCoeff * (dS(i) - dP(i)) +
                         robinAlpha * vp(i) - robinAlpha * uf(i));

      return value;
    });

    // ----------------------------------------------------------------------
    // Solid Newmark / NeoHookean problem (nonlinear: SNES on the increment
    // etaState, with dState = dOld + etaState).
    // ----------------------------------------------------------------------
    Problem solid(d, w);
    solid =
        solidMass * Integral(d, w) + solidTangent +
        solidMass * Integral(dState, w) - solidMass * Integral(dPred, w)
        + solidInternal
        + DirichletBC(d, zero).on(BoundarySolid::Inlet, BoundarySolid::Outlets[0], BoundarySolid::Outlets[1], BoundarySolid::Outlets[2], BoundarySolid::Outlets[3], BoundarySolid::Outlets[4], BoundarySolid::Outlets[5])
        + DirichletBC(d, zero).on(BoundarySolid::FSIRing)
        // Robin-Robin transmission (solid side), per current area (J_a):
        //   sigma_s n_s + alpha d_dot = alpha u_f^{lag} + t_f^{lag},
        //   d_dot = vPred + solidVelocityCoeff (dState - dPred)  (Newmark).
        + robinVelocityCoeff *
              BoundaryIntegral(arealStretch * Dot(d, w))
                  .over(BoundarySolid::FSI) +
        BoundaryIntegral(robinInterfaceData, w).over(BoundarySolid::FSI)
        // Neumann load: the traction exerted by the fluid on the solid is
        //   t_f = sigma_f n_s = -sigma_f n_f = fluidStress  (p>0 pushes the
        // wall OUTWARD).
        - BoundaryIntegral(fluidStress, w)
            .over(BoundarySolid::FSI);

    solid.assemble();
    Solver::KSP kspSolid(solid);
    Solver::SNES snes(kspSolid);
    snes.setTolerances(1.0e-10, 1.0e-8, 1.0e-10, 50, 10000);
    snes.setStateUpdate([&](const PETSc::Math::Vector &state) {
      etaState.setData(state, 0);
      dState = dOld;
      dState += etaState;
    });

    // ----------------------------------------------------------------------
    // Quasi-static wall prestress (cfg.prestressSteps > 0): ramp a uniform
    // static NeoHookean increments.  Lumen pressure pushes the wall OUTWARD:
    // on the inner (FSI) surface the solid outward normal n_s points INTO the
    // lumen, so the traction is t = -p n_s, entering the residual as
    //   -int t.w = + int p (n_s . w).
    // ----------------------------------------------------------------------
    if (cfg.prestressSteps > 0) {
      // Prestress to the FULL initial arterial pressure so the static wall
      // equilibrium is (approximately) balanced by the transferred fluid
      // traction at the dynamic start (pCur = par), giving a small first
      // residual.  (The previous 0.9 factor left a 0.1*par deficit that kicked
      // the wall on step 1 -> ALE mesh-velocity spike -> the "awful" fluid
      // velocity field reported at the static->dynamic handoff.)
      const Real p0 = model.getState().par;
      Real prestressPressure = 0.0;

      PETSc::Variational::TrialFunction dPre(dh);
      PETSc::Variational::TestFunction wPre(dh);
      Solid::MaterialTangent preTangent(law, dPre, wPre, dState);
      Solid::InternalForce preInternal(law, wPre, dState);

      // Follower pressure (exact deformed-surface load + consistent
      // tangent): full quadratic Newton up to total pressure.
      Solid::FollowerPressureForce preLoad(prestressPressure, wPre, dState);
      preLoad.over(BoundarySolid::FSI);
      Solid::FollowerPressureTangent preLoadK(prestressPressure, dPre, wPre,
                                              dState);
      preLoadK.over(BoundarySolid::FSI);

      Problem prestress(dPre, wPre);
      prestress =
          preTangent + preInternal + preLoadK + preLoad +
          DirichletBC(dPre, zero).on(
              BoundarySolid::Inlet, BoundarySolid::Outlets[0],
              BoundarySolid::Outlets[1], BoundarySolid::Outlets[2],
              BoundarySolid::Outlets[3], BoundarySolid::Outlets[4],
              BoundarySolid::Outlets[5], BoundarySolid::FSIRing);

      prestress.assemble();
      Solver::KSP kspPre(prestress);
      Solver::SNES snesPre(kspPre);
      snesPre.setTolerances(1.0e-10, 1.0e-8, 1.0e-10, 50, 10000);
      snesPre.setStateUpdate([&](const PETSc::Math::Vector &state) {
        etaState.setData(state, 0);
        dState = dOld;   // dOld stays zero: etaState accumulates the full lift
        dState += etaState;
      });

      for (size_t k = 1; k <= cfg.prestressSteps; ++k) {
        prestressPressure =
            (static_cast<Real>(k) / static_cast<Real>(cfg.prestressSteps)) * p0;
        snesPre.solve();
        if (!snesPre.converged()) {
          if (isRoot)
            std::cerr << "Prestress SNES failed at increment " << k << " / "
                      << cfg.prestressSteps << "; continuing with the last "
                      << "converged (partial) prestress state.\n";
          break;
        }
      }
      if (isRoot)
        Alert::Info() << "Prestressed wall to " << prestressPressure
                      << " Pa in " << cfg.prestressSteps << " increment(s)"
                      << Alert::Raise;

      // Commit the prestressed state as the dynamic initial condition (zero
      // velocity / acceleration are already set) and lift the fluid mesh.
      dOld.setData(dState.getData());
      dIter.setData(dState.getData());
      restoreMeshToReference(meshFluid, referenceVertices);
      ale.assemble();
      Solver::KSP(ale).solve();
      aleDisp.setData(dMove.getSolution().getData());
      aleDispOld.setData(aleDisp.getData());
      moveMeshWithVertexDisplacement(meshFluid, referenceVertices, uh,
                                     aleDisp);

      // Pressurized fluid start: seed the fluid pressure at par so the
      // transferred wall traction starts in (approximate) balance with the
      // prestressed wall equilibrium (no follower level any more).
      pOld = model.getState().par;
      pCur = model.getState().par;
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
    //   (needed when the added-mass effect destabilizes the loose scheme).
    // ======================================================================

    for (size_t step = 1; step <= cfg.nsteps; ++step) {
      const auto rep = model.step(dt);
      if (!rep.converged) {
        if (isRoot)
          std::cerr << "0D model failed at step " << step << '\n';
        break;
      }

      const auto &s = model.getState();

      // Startup load ramp ("rampage").  With the prestress the wall is already
      // in pressurized equilibrium, so the loads are imposed at full strength
      // from step 1 (loadRamp == 1).  WITHOUT the prestress (prestressSteps ==
      // 0) the lumen is inflated from the unloaded reference state by climbing
      // the imposed inlet/outlet pressures 0 -> 1 over rampSteps steps.
      //
      // A SMOOTHSTEP profile (C^1, zero slope at both ends) is used instead of
      // a linear ramp: it keeps d/dt of the load small, so the compliant wall
      // inflates quasi-statically and the fluid never builds the large inflow
      // velocity (-> Bernoulli suction / negative interior pressure) that a
      // fast ramp produces.  Inflating a compliant vessel from rest is an
      // intrinsically violent transient; ramp SLOWLY (rampSteps large) or,
      // better, use the prestress, which inflates the wall quasi-statically
      // BEFORE the dynamics so this transient never happens.
      Real loadRamp = 1.0;
      if (cfg.prestressSteps == 0 && cfg.rampSteps > 0) {
        const Real sRamp =
            std::min(Real(1.0),
                     static_cast<Real>(step) / static_cast<Real>(cfg.rampSteps));
        loadRamp = sRamp * sRamp * (3.0 - 2.0 * sRamp); // smoothstep
      }

      // Physiological pressures (ramped on a no-prestress start).
      pinValue = loadRamp * s.par;
      for (const auto &[tag, bc] : wk) {
        // Full windkessel pressure pout = pc + Rp*Q (lagged), imposed as a
        // Neumann traction -- identical to CoupledLV0DCoronary3D, the proven
        // 0D-3D outlet treatment.  Stability requires a PHYSIOLOGICAL Rp (see
        // RCR::Rp); an over-stiff Rp makes pout hypersensitive to Q and the
        // outlet oscillates.  The same loadRamp scales the outlet level so the
        // inlet-to-outlet pressure gradient inflates consistently.
        outletPressureValue[tag] =
            loadRamp * (s.par - cfg.pressureDropScale * (s.par - bc.pout));
      }
      // No direct wall load: the wall is driven purely by the fluid traction,
      // which itself follows the ramped inlet/outlet pressures, so the wall
      // still inflates gradually on a no-prestress start.

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

      if (isRoot) {
        Alert::Info() << "Coronary explicit ALE FSI step " << step << " / "
                      << cfg.nsteps
                      << "  (coupling iterations: " << cfg.couplingIterations
                      << ")" << Alert::Raise;
      }

      std::map<Attribute, Real> qOut;
      for (const Attribute outlet : BoundaryFluid::Outlets)
        qOut[outlet] = 0.0;
      Real qOutSum = 0.0;
      Real qIn = 0.0;
      bool stepFailed = false;
      Real lastRel = 0.0;     // final relative interface change (coupling)
      size_t couplesDone = 0; // coupling iterates actually performed

      // Omega^n and u^n are fixed within the step: assemble massOld once.
      moveMeshWithVertexDisplacement(meshFluid, referenceVertices, uh,
                                     aleDispOld);
      massOld.assemble();

      {
        for (size_t couple = 1; couple <= cfg.couplingIterations; ++couple) {
          // Solid Newmark / NeoHookean solve on the reference mesh.
          snes.solve();
          if (!snes.converged()) {
            if (isRoot)
              std::cerr << "Solid SNES failed at step " << step
                        << " (coupling iterate " << couple << ") after "
                        << snes.getIterationNumber() << " iterations.\n";
            stepFailed = true;
            break;
          }

          // Full (un-relaxed) interface update: dIter <- dState.  NO Aitken
          // relaxation -- loose Robin-Robin coupling (couplingIterations == 1
          // by default; for > 1 this is a plain Picard fixed point).  rel is
          // kept purely as a per-step diagnostic.
          auto delta = dState;
          delta -= dIter;
          PetscReal deltaNorm = 0.0;
          PetscReal stateNorm = 0.0;
          VecNorm(delta.getData(), NORM_2, &deltaNorm);
          VecNorm(dState.getData(), NORM_2, &stateNorm);
          const Real rel = (stateNorm > 0.0)
                               ? (static_cast<Real>(deltaNorm) /
                                  static_cast<Real>(stateNorm))
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

          if (isRoot) {
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
          // Re-assemble on the reference mesh so the FSI DirichletBC picks up
          // the CURRENT solid displacement iterate 'dIter' (see note above).
          ale.assemble();
          Solver::KSP(ale).solve();
          aleDisp.setData(dMove.getSolution().getData());


          // ALE mesh velocity.
          meshVelocity = aleDisp;
          meshVelocity -= aleDispOld;
          meshVelocity *= 1.0 / dt;

          moveMeshWithVertexDisplacement(meshFluid, referenceVertices, uh,
                                         aleDisp);

          // Refresh the ALE convecting velocity uConv = u^n - w and re-project
          // the VMS fields (projected convection, tau, dynamic subscale) on the
          // just-moved mesh, so the stabilization is consistent with the Oseen
          // convection assembled below.  Order matters: up -> tau -> subscale
          // (the subscale update reads both up and tau).
          uConv = uOld;
          uConv -= meshVelocity;
          vmsL2Conv.assemble();
          Solver::KSP(vmsL2Conv).solve();
          vmsTauProj.assemble();
          Solver::KSP(vmsTauProj).solve();
          vmsSubProj.assemble();
          Solver::KSP(vmsSubProj).solve();
          // Orthogonal grad-div IMEX projections.  Project sqrt(tau_C^{n+1})
          // and form the bilinear coefficient as its EXACT pointwise square
          // tau_C = (sqrt tau_C)^2 (so the two are consistent).  Then the
          // lagged pi~^n = Pi( sqrt(tau_C^n) div u^n ) (depends only on time-n
          // data: uOld, vmsSqrtTauCOld; re-projected on the moved mesh).
          vmsSqrtTauCProj.assemble();
          Solver::KSP(vmsSqrtTauCProj).solve();
          VecPointwiseMult(vmsTauC.getSolution().getData(),
                           vmsSqrtTauC.getSolution().getData(),
                           vmsSqrtTauC.getSolution().getData());
          vmsPiTildeProj.assemble();
          Solver::KSP(vmsPiTildeProj).solve();

          flow.assemble().setFieldSplits();

          // Inject (rho/dt)(u^n, v)|_{Omega^n} into the velocity block.
          {
            ::Vec b = flow.getLinearSystem().getVector();
            const PetscInt vOff = static_cast<PetscInt>(
                flow.getTestOffsets()[0]); // velocity block
            const ::Vec &mOld = massOld.getVector();
            PetscInt lo = 0, hi = 0;
            VecGetOwnershipRange(mOld, &lo, &hi);
            const PetscScalar *arr = nullptr;
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
          tractionTransfer.project(Region::Faces, tractionFSI,
                                   BoundaryFluid::FSI);
          uWall.project(Region::Faces, uCur, BoundaryFluid::FSI);

        }
      }

      if (stepFailed)
        break;

      // Converged interface displacement and consistent kinematics.
      dState = dIter;
      solidAcceleration = dState;
      solidAcceleration -= dPred;
      solidAcceleration *= 1.0 / (betaN * dt * dt);
      solidVelocity = vPred;
      tmp = solidAcceleration;
      tmp *= gammaN * dt;
      solidVelocity += tmp;

      flux = BoundaryIntegral(Dot(uCur, normalFluid), qFlux)
                 .over(BoundaryFluid::Inlet);
      flux.assemble();
      qIn = flux(one);
      qOutSum = 0.0;
      for (const Attribute outlet : BoundaryFluid::Outlets) {
        flux = BoundaryIntegral(Dot(uCur, normalFluid), qFlux).over(outlet);
        flux.assemble();
        const Real qo = flux(one);
        qOut[outlet] = qo;
        qOutSum += qo;
      }

      // ---- Interface power (Robin-Robin energy consistency) -------------
      // E_Gamma = int_FSI (sigma_f . n_f) . (u_f - d_dot_s) dGamma, the rate
      // of work the fluid does against the kinematic mismatch on Gamma_FSI.
      // For an exactly energy-conserving (and conservatively-coupled) FSI this
      // is 0; with the loose Robin-Robin scheme AND the NON-conservative
      // cross-mesh NODAL transfer it is nonzero.  Comparing it between the
      // P2/P1 and P1/P1 cases isolates how much of that energy comes from the
      // interpolation.  sigma_f . n_f = -tractionFSI = muFsi(grad u+grad u^T)n
      // - p n.  Computed as two integrals to avoid GridFunction-Function
      // arithmetic inside the expression.
      const auto sigmaFn =
          muFsi * Mult(strainRateFsi, normalFluid) - pCur * normalFluid;
      flux =
          BoundaryIntegral(Dot(sigmaFn, uCur), qFlux).over(BoundaryFluid::FSI);
      flux.assemble();
      const Real ePowerFluid = flux(one);
      flux = BoundaryIntegral(Dot(sigmaFn, interfaceSolidVelocity), qFlux)
                 .over(BoundaryFluid::FSI);
      flux.assemble();
      const Real ePowerSolid = flux(one);
      const Real eInterface = ePowerFluid - ePowerSolid;

      // ---- Startup/stability diagnostics --------------------------------
      // p range: a strongly NEGATIVE interior pressure means the wall is
      // inflating faster than the inflow can fill -> Bernoulli/inertial
      // suction (ramp too fast, or use the prestress).
      // mass = qIn + qOutSum = -d|lumen|/dt (no-slip wall): < 0 while the wall
      // inflates, -> 0 once the lumen reaches equilibrium.  If |mass| does NOT
      // decay toward 0 as loadRamp -> 1, suspect a GCL / mass-conservation
      // problem rather than a physical transient.
      {
        PetscReal pMin = 0.0, pMax = 0.0;
        VecMin(pCur.getData(), PETSC_NULLPTR, &pMin);
        VecMax(pCur.getData(), PETSC_NULLPTR, &pMax);
        if (isRoot)
          Alert::Info() << "  [diag] p in [" << pMin << ", " << pMax << "] Pa"
                        << "  mass(qIn+qOut) = " << (qIn + qOutSum)
                        << "  E_iface = " << eInterface << " W"
                        << "  | coupling: iters = " << couplesDone << "/"
                        << cfg.couplingIterations
                        << "  interface change = " << lastRel << Alert::Raise;
      }

      // RCR / Windkessel update from the converged interface flux.
      for (const Attribute outlet : BoundaryFluid::Outlets)
        updateRCRNonNew(cfg, outlet, model, wk[outlet], qOut[outlet], dt);

      // Commit the step-(n+1) state.
      uOld.setData(uCur.getData());
      pOld.setData(pCur.getData());
      dOld.setData(dState.getData());
      solidVelocityOld.setData(solidVelocity.getData());
      solidAccelerationOld.setData(solidAcceleration.getData());
      // Carry the dynamic VMS subscale u'^{n+1} -> u'^n for the next step.
      vmsSubOld.setData(vmsSub.getSolution().getData());
      // Carry sqrt(tau_C^{n+1}) -> sqrt(tau_C^n) for the lagged pi~ projection.
      vmsSqrtTauCOld.setData(vmsSqrtTauC.getSolution().getData());


      aleDispOld.setData(aleDisp.getData());


      xdmf_fluid.write(s.t).flush();
      xdmf_solid.write(s.t).flush();

      if (isRoot && csv) {
        csv << s.t << ',' << s.y << ',' << s.v << ',' << s.pv << ',' << s.par
            << ',' << s.pd << ',' << qIn << ',' << qOutSum;
        for (const Attribute outlet : BoundaryFluid::Outlets)
          csv << ',' << qOut[outlet] << ',' << wk[outlet].pout;
        csv << ',' << eInterface << ',' << lastRel << '\n';
        csv.flush();
      }
    }

    xdmf_fluid.close();
    xdmf_solid.close();
    if (isRoot && csv)
      csv.close();
  } catch (const std::exception &e) {
    std::cerr << "CoronaryArtery_FSI_Explicit_PETSc_Seq failed: " << e.what()
              << '\n';
    PetscFinalize();
    return 1;
  }

  PetscFinalize();
  return 0;
}
