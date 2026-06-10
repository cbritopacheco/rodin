/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * Explicit (staggered) PETSc SEQUENTIAL coronary ALE FSI prototype.
 *
 * Single-process (non-MPI) sibling of CoronaryArtery_FSI_Explicit_PETSc_MPI.cpp:
 * identical physics and numerics, but the mesh lives in Context::Local (no
 * boost::mpi environment, no Sharder/partitioner, no distributed interface).
 * PETSc is still used for the linear/nonlinear solves (KSP/SNES) and runs in
 * serial, exactly as in examples/PETSc/PDEs/Seq_BDF1_ALE_FSI.cpp.  Because the
 * fluid<->solid interface coupling was only ever valid on a single rank, this
 * version drops the MPI scaffolding outright instead of guarding against it.
 *
 * Same physics (coronary Carreau-Yasuda fluid, NeoHookean/Newmark
 * solid, harmonic ALE, 0D CCMLC2014 heart model + RCR outlets) but a
 * Robin-Neumann staggered coupling instead of a monolithic SNES.
 *
 * Staggered step n -> n+1 (per coupling sub-iteration, default 1 = loose):
 *   1. Advance the 0D heart model; read inlet/outlet pressures.
 *   2. Newmark predictors dPred, vPred.
 *   3. SOLID solve (total Lagrangian, reference config): NeoHookean/Newmark
 *      SNES with the fluid interface traction as Neumann load plus a Robin
 *      transmission term on Gamma_FSI.
 *   4. Harmonic ALE extension of the solid displacement iterate over the fluid;
 *      meshVelocity = (aleDisp - aleDispOld) / dt; move the fluid mesh.
 *   5. FLUID Oseen ALE solve (u, p) on the moved mesh, with the Robin
 *      transmission of the solid interface velocity on Gamma_FSI.
 *   6. Under-relax the interface displacement, test convergence.
 *   7. Commit state, compute inlet/outlet fluxes, update the RCR models.
 *
 * Boundary/interface attributes: 2 = FSI wall/interface, plus inlet/outlet caps
 * and a one-element FSI ring band (attr 99) clamped at the cut ends.
 *
 * IMPORTANT - added-mass instability: loosely coupled FSI is only conditionally
 * stable and worst when fluidDensity == solidDensity (here both 1060). For
 * production increase -coronary_coupling_iterations with
 * -coronary_coupling_relaxation < 1, and/or reduce -coronary_traction_scale.
 */

#include "Rodin/Variational/BoundaryIntegral.h"
#include "Rodin/Variational/ForwardDecls.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cmath>
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

#include "CoronaryArtery/VMSConvectionIntegrator.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Math;
using namespace Rodin::Solver;
using namespace Rodin::Variational;
using namespace Rodin::Examples::Heart;

namespace {
using Model = Rodin::Heart::CCMLC2014T<>;
using MeshType = Geometry::Mesh<Context::Local>;

struct BoundaryFluid {
  static constexpr Attribute FSI = 2;
  static constexpr Attribute Inlet = 4;
  static constexpr std::array<Attribute, 6> Outlets{{7, 8, 9, 10, 14, 15}};
  // One-element-wide band of FSI faces touching the inlet/outlet caps. Promoted
  // at startup from FSI so the fluid velocity can be strongly pinned to zero on
  // the cap rings while the bulk FSI keeps its Robin transmission.
  static constexpr Attribute FSIRing = 99;
};

struct BoundarySolid {
  static constexpr Attribute FSI = 2;
  static constexpr Attribute Inlet = 17;
  static constexpr std::array<Attribute, 6> Outlets{{18, 19, 20, 21, 22, 31}};
  static constexpr Attribute FSIRing = 99;
};

struct RCR {
  Real Rp = 5.0e9;
  Real C = 1.0e-12;
  Real Rd = 5.0e10;
  Real pd = 500.0;
  Real pc = 10600.0;
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
  std::string xdmfBasename = "CoronaryArtery_FSI_Explicit";
  std::string csvPath = "CoronaryArtery_FSI_Explicit.csv";
  Real meshScale = 1.0e-3;

  Real dt = 1.0e-3;
  size_t nsteps = 3 * static_cast<int>(0.85 / 1.0e-3);

  // Startup load ramp: smoothly raise the imposed inlet/outlet pressures from
  // 0 to their full value over this window so the un-prestressed wall eases
  // into the transmural load instead of being step-loaded in a single dt
  // (which makes meshVelocity = (aleDisp - aleDispOld)/dt blow up).  With the
  // Robin-Robin wall (lagged-traction data, see the flow problem) the domain
  // pressurizes over the first few steps as sigma^{lag} accumulates, so the
  // ramp's only job is to ease the wall load; a window of ~50-100 steps keeps
  // the loose coupling well-behaved.
  Real pressureRampTime = 0.05;

  // Scale applied to the inlet->outlet pressure drop.  The 0D model produces a
  // transmural gradient of ~1000 Pa between the inlet (arterial) pressure and
  // the outlet (Windkessel) pressures; this compresses each outlet toward the
  // inlet so the driving gradient is pressureDropScale * original (0.5 -> ~500
  // Pa).  1.0 leaves the gradient untouched.
  Real pressureDropScale = 0.5;

  Real fluidDensity = 1060.0;
  Real pressurePenalty = 1.0e-12;
  Real inletBackflowStabilization = 1.0;
  Real outletBackflowStabilization = 1.0;
  CarreauYasuda viscosity;
  // Geometry + nonlinear-solve parameters for the non-Newtonian outlet flow
  // law used by updateRCRNonNew (per-outlet proximal/distal surrogate vessels).
  OutletFlowLaw outletFlowLaw;

  Real inletImpedance = 1.e3;
  Real inletTangentialDamping = 1.e3;

  // Quasi-static wall prestress: before the dynamic loop, ramp a uniform
  // lumen pressure 0 -> par(0) on the solid FSI surface in this many static
  // NeoHookean increments (0 = disabled).  The wall then enters the time loop
  // already in equilibrium with the transmural load (no inlet ballooning),
  // the fluid mesh is lifted to match, the fluid pressure is initialized so
  // the first traction exchange is balanced, and the startup pressure ramp
  // is skipped.  NOTE: this is forward pressurization (the geometry ends
  // slightly inflated w.r.t. the imaged one), not an inverse prestress.
  size_t prestressSteps = 60;

  Real solidDensity = 1060.0;
  Real solidYoungModulus = 5.0e5;
  Real solidPoissonRatio = 0.4;
  Real newmarkBeta = 0.25;
  Real newmarkGamma = 0.5;

  // Explicit-coupling specific options.
  Real tractionScale = 1.0; // overall scale of the transferred traction
  Real tractionPressureScale =
      1.0;                         // scale of the pressure part (Seq used 0.05)
  Real tractionViscousScale = 1.0; // scale of the viscous part
  size_t couplingIterations = 2;   // 1 = loosely coupled; >1 = strong coupling
  Real couplingRelaxation = 1.0;   // displacement under-relaxation for k>1
  Real couplingTolerance = 1.0e-6; // relative interface-displacement tolerance

  Real robinAlpha = 1e4;
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

// Linear (Newtonian) RCR update; superseded by updateRCRNonNew below but kept
// for reference / fallback.
[[maybe_unused]] static void updateRCR(const Model &model, RCR &bc, Real Q,
                                       Real dt) {
  const Real cap = bc.C / dt;
  const auto &s = model.getState();

  bc.pc = (cap * bc.pc + Q + s.pv / bc.Rd) / (cap + 1.0 / bc.Rd);
  bc.qd = (bc.pc - bc.pd) / bc.Rd;
  bc.pout = bc.pc + bc.Rp * Q;
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
  // No reconcile(): that synchronizes shared entities across MPI ranks and is
  // only defined on the distributed mesh.  A local mesh needs no reconciliation.

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
// the remaining (bulk) FSI faces.  A face is in the band iff it shares at least
// one vertex with a cap face, i.e. it touches the ring curve.  Done identically
// on the fluid and solid reference meshes (assumed conforming on the interface)
// so the bulk FSI faces still pair up one-to-one in buildInterfaceMap.
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
// In this sequential build the whole interface lives on one process, so the
// match is exact and complete (no cross-rank concerns).
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

  Point q(p);
  q.setPolytope(*solidFace);

  return q;
}

static Point forwardSolidPointToFluid(const Point &p, const MeshType &fluidMesh,
                                      const InterfaceMap &map) {
  const auto found = map.solidToFluid.find(p.getPolytope().getIndex());

  if (found == map.solidToFluid.end())
    throw std::runtime_error("Solid point is not on a mapped FSI face.");

  auto fluidFace = fluidMesh.getFace(found->second);

  Point q(p);
  q.setPolytope(*fluidFace);

  return q;
}

static std::string lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return s;
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

  PetscReal tractionScale = cfg.tractionScale;
  PetscBool tractionScaleSet = PETSC_FALSE;
  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_traction_scale",
                      &tractionScale, &tractionScaleSet);
  if (tractionScaleSet)
    cfg.tractionScale = tractionScale;

  PetscReal tractionPressureScale = cfg.tractionPressureScale;
  PetscBool tractionPressureScaleSet = PETSC_FALSE;
  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
                      "-coronary_traction_pressure_scale",
                      &tractionPressureScale, &tractionPressureScaleSet);
  if (tractionPressureScaleSet)
    cfg.tractionPressureScale = tractionPressureScale;

  PetscReal tractionViscousScale = cfg.tractionViscousScale;
  PetscBool tractionViscousScaleSet = PETSC_FALSE;
  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR,
                      "-coronary_traction_viscous_scale", &tractionViscousScale,
                      &tractionViscousScaleSet);
  if (tractionViscousScaleSet)
    cfg.tractionViscousScale = tractionViscousScale;

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
}

// Non-Newtonian (Carreau-Yasuda) RCR outlet update.  Replaces the linear
// 'updateRCR': the proximal and distal surrogate vessels each carry a
// shear-thinning Poiseuille-type flow law (wall-shear root solve + RK4 WRMS
// integral), and the distal capacitor pressure pc is advanced implicitly by a
// Newton solve.  'tag' selects the per-outlet geometry from
// cfg.outletFlowLaw.geometricParam.
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


    VelocityFES uph(std::integral_constant<size_t, 2>{}, meshFluid, dimFluid);
    PressureFES tauh(std::integral_constant<size_t, 1>{}, meshFluid);

    PETSc::Variational::TrialFunction upTrial(uph);  // Pi[(grad u^n) u^n]
    PETSc::Variational::TestFunction vp(uph);
    PETSc::Variational::TrialFunction subTrial(uph); // dynamic subscale u'^{n+1}
    PETSc::Variational::TrialFunction tauTrial(tauh); // projected tau_K
    PETSc::Variational::TestFunction tp(tauh);
    PETSc::Variational::TrialFunction muTrial(tauh); // projected mu(u^n)
    PETSc::Variational::TestFunction mp(tauh);

    PETSc::Variational::GridFunction subOld(uph); // u'^n

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
    PETSc::Variational::GridFunction dIter(dh); // relaxed coupling iterate
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
    subOld = zero;

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
      csv << '\n';
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

    const Real robinVelocityCoeff = cfg.robinAlpha * solidVelocityCoeff;

    const Real tractionPScale = cfg.tractionScale * cfg.tractionPressureScale;
    const Real tractionVScale = cfg.tractionScale * cfg.tractionViscousScale;

    const auto &cy = cfg.viscosity;
    const Real gammaReg = cy.gammaRegularization;
    const Real deltaMu = cy.mu0 - cy.muInf;

    const auto normalFluid = BoundaryNormal(meshFluid);

    Solid::NeoHookean law(solidLambda, solidMu);
    Solid::MaterialTangent solidTangent(law, d, w, dState);
    Solid::InternalForce solidInternal(law, w, dState);

    // Current fluid solution (read back after each Oseen solve).
    PETSc::Variational::GridFunction uCur(uh);
    PETSc::Variational::GridFunction pCur(ph);
    uCur = zero;
    pCur = 0.0;

    const auto gradUfsi = Jacobian(uCur);
    const auto strainRateFsi = gradUfsi + Transpose(gradUfsi);
    const auto symUfsi = 0.5 * strainRateFsi;
    const auto shearFsi =
        Sqrt(gammaReg * gammaReg + 2.0 * Dot(symUfsi, symUfsi));
    const auto muFsi =
        cy.muInf + deltaMu * Pow(1.0 + Pow(cy.lambda * shearFsi, cy.yasuda),
                                 (cy.n - 1.0) / cy.yasuda);

    // PHYSICAL (unscaled) fluid traction t_f = -sigma_f n_f: the fluid's OWN
    // lagged Robin datum sigma_f^{lag} n and the traction output field.
    const auto tractionFSI =
        (1.0 * pCur) * normalFluid -
        (1.0 * muFsi) * Mult(strainRateFsi, normalFluid);

    // Scaled traction TRANSFERRED to the solid
    const auto tractionToSolid =
        (tractionPScale * pCur) * normalFluid -
        (tractionVScale * muFsi) * Mult(strainRateFsi, normalFluid);

    // Areal stretch J_a = A_t / A_0 of a solid FSI face, evaluated at the
    // CURRENT coupling iterate 'dIter' (not the lagged dOld): converts
    // per-current-area (Cauchy) interface data to the per-reference-area
    // measure of the total-Lagrangian solid boundary integrals.  Shared by
    // the fluid traction pullback ('fluidStress') and the solid Robin
    // alpha-terms, so the transmission parameter alpha acts on the SAME
    // (current) area measure on both sides of the interface.  'dIter' is
    // frozen during each SNES solve, so the stretch is explicit data and
    // contributes nothing to the Jacobian.
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
      const auto force = tractionToSolid(xf);
      for (Index i = 0; i < static_cast<Index>(dim); ++i)
        value(i) = force(i);

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

      return value;
    });

    const auto transportLag = uOld - meshVelocity;
    const auto convU = Mult(Jacobian(u), transportLag);
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

    // NOTE: the BDF1 time term has TWO parts that live on DIFFERENT meshes:
    //   (rho/dt) [ \int_{Omega^{n+1}} u . v  -  \int_{Omega^n} u^n . v ].
    // The implicit part (over Omega^{n+1}) stays in 'flow' below.
    Problem flow(u, p, v, q);
    flow =
        (cfg.fluidDensity / dt) * Integral(u, v) +
        cfg.fluidDensity * Integral(Dot(convU, v)) +
        0.5 * cfg.fluidDensity * Integral(divGeomTemam * Dot(u, v)) +
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

        + cfg.inletImpedance *
              BoundaryIntegral(Dot(Dot(u, normalFluid) * normalFluid, v))
                  .over(BoundaryFluid::Inlet)
        + cfg.inletTangentialDamping *
              BoundaryIntegral(Dot(duTangential, v)).over(BoundaryFluid::Inlet)

        + cfg.robinAlpha * BoundaryIntegral(u,v).over(BoundaryFluid::FSI)
        - cfg.robinAlpha * BoundaryIntegral(interfaceSolidVelocity,v).over(BoundaryFluid::FSI)
        + BoundaryIntegral(tractionFSI,v).over(BoundaryFluid::FSI)

        + DirichletBC(u, zero).on(BoundaryFluid::FSIRing);
        // Projected dynamic-VMS convective stabilization (lagged): bilinear
        // streamline term plus its projection/subscale source.  tauTrial,
        // upTrial and subTrial are L2-projected each step just before assembly.
       // + VMSConvectionBilinearIntegrator(u, v, uOld, tauTrial.getSolution(),
         //                                 cfg.fluidDensity)
        //- VMSConvectionLinearIntegrator(v, subTrial.getSolution(), uOld,
          //                              upTrial.getSolution(),
            //                            tauTrial.getSolution(), cfg.fluidDensity,
              //                          dt);

    PETSc::Variational::TestFunction vMass(uh);
    LinearForm<VelocityFES, ::Vec> massOld(vMass);
    massOld = (cfg.fluidDensity / dt) * Integral(uOld, vMass);

    // ---- VMS projection forms (mass-matrix L2 projections, reassembled each
    // step on the current configuration before the fluid solve).
    //   convectionTarget = (grad u^n) u^n  (lagged convective acceleration).
    const auto convectionTarget = Mult(Jacobian(uOld), uOld);

    // Projected lagged viscosity (feeds tau); muLag is the Carreau-Yasuda
    // viscosity evaluated at u^n, already defined for the Oseen operator above.
    Problem muProj(muTrial, mp);
    muProj = Integral(muTrial, mp) - Integral(muLag, mp);

    // Projected convective acceleration Pi[(grad u^n) u^n].
    Problem l2ConvU(upTrial, vp);
    l2ConvU = Integral(upTrial, vp) - Integral(convectionTarget, vp);

    // Stabilization parameter tau_K, projected onto the P1 scalar space.
    RealFunction tauVMS = [&](const Point &pt) -> Real {
      const auto uo = uOld.getValue(pt);
      const Real muv = muTrial.getSolution().getValue(pt);
      const Real hK = std::pow(pt.getPolytope().getMeasure(),
                               1.0 / pt.getPolytope().getDimension());
      const Real order = 2.0;
      const Real speed = std::sqrt(Math::dot(uo, uo));
      const Real Tau =
          1.0 / (4.0 * std::pow(order, 4.0) * muv /
                     (cfg.fluidDensity * hK * hK) +
                 2.0 * order * speed / hK);
      return 1.0 / (cfg.fluidDensity / dt + cfg.fluidDensity / Tau);
    };
    Problem tauProj(tauTrial, tp);
    tauProj = Integral(tauTrial, tp) - Integral(tauVMS, tp);

    // Dynamic subscale u'^{n+1} = tau rho ( u'^n / dt - (conv - Pi[conv]) ).
    auto subUpdate = VectorFunction(dim, [&](const Point &pt) {
      const auto conv = convectionTarget.getValue(pt);
      const auto proj = upTrial.getSolution().getValue(pt);
      const auto old = subOld.getValue(pt);
      const Real tv = tauTrial.getSolution().getValue(pt);
      Math::SpatialVector<Real> out(dim);
      for (Index c = 0; c < static_cast<Index>(dim); ++c)
        out(c) = tv * cfg.fluidDensity *
                 (1.0 / dt * old(c) - (conv(c) - proj(c)));
      return out;
    });
    Problem subProj(subTrial, vp);
    subProj = Integral(subTrial, vp) - Integral(subUpdate, vp);


    auto interfaceSolidDisplacement = VectorFunction(dim, [&](const Point& xf)
    {
      const Point xs =
        forwardFluidPointToSolid(xf, meshSolid, interfaceMap);

      Math::SpatialVector<Real> value(dim);
      value.setZero();

      const auto ds = dIter(xs);
      for (Index i = 0; i < static_cast<Index>(dim); ++i)
        value(i) = ds(i);

      return value;
    });

    // ----------------------------------------------------------------------
    // Harmonic ALE mesh-extension problem.
    //   * harmonic in the fluid,
    //   * equal to the solid displacement iterate in the solid block,
    //   * equal to the solid displacement iterate on the FSI interface,
    //   * fixed (zero) at the inlet / outlet rings.
    // ----------------------------------------------------------------------
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


    auto robinInterfaceData = VectorFunction(dim, [&](const Point &xs) {
      const Point xf = forwardSolidPointToFluid(xs, meshFluid, interfaceMap);
      const Real Ja = arealStretchAt(xs);

      Math::SpatialVector<Real> value(dim);
      value.setZero();

      const auto uf = uCur(xf);
      const auto vp = vPred(xs);
      const auto dS = dState(xs);
      const auto dP = dPred(xs);
      for (Index i = 0; i < static_cast<Index>(dim); ++i)
        value(i) = Ja * (robinVelocityCoeff * (dS(i) - dP(i)) +
                         cfg.robinAlpha * vp(i) - cfg.robinAlpha * uf(i));

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

        + robinVelocityCoeff *
              BoundaryIntegral(arealStretch * Dot(d, w))
                  .over(BoundarySolid::FSI) +
        BoundaryIntegral(robinInterfaceData, w).over(BoundarySolid::FSI)
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
    // lumen pressure 0 -> tractionPScale * par(0) on the solid FSI surface in
    // static NeoHookean increments.
    // ----------------------------------------------------------------------
    if (cfg.prestressSteps > 0) {
      const Real p0 = tractionPScale * 11000.;
      Real prestressPressure = 0.0;
      const auto normalSolid = BoundaryNormal(meshSolid);
      auto pPre = RealFunction([&](const Point &) { return prestressPressure; });

      PETSc::Variational::TrialFunction dPre(dh);
      PETSc::Variational::TestFunction wPre(dh);
      Solid::MaterialTangent preTangent(law, dPre, wPre, dState);
      Solid::InternalForce preInternal(law, wPre, dState);

      Problem prestress(dPre, wPre);
      prestress =
          preTangent + preInternal +
          DirichletBC(dPre, zero).on(
              BoundarySolid::Inlet, BoundarySolid::Outlets[0],
              BoundarySolid::Outlets[1], BoundarySolid::Outlets[2],
              BoundarySolid::Outlets[3], BoundarySolid::Outlets[4],
              BoundarySolid::Outlets[5], BoundarySolid::FSIRing) +
          BoundaryIntegral(pPre * Dot(wPre, normalSolid))
              .over(BoundarySolid::FSI);

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

      dOld.setData(dState.getData());
      dIter.setData(dState.getData());

      restoreMeshToReference(meshFluid, referenceVertices);
      ale.assemble();
      Solver::KSP(ale).solve();
      aleDisp.setData(dMove.getSolution().getData());
      aleDispOld.setData(aleDisp.getData());
      moveMeshWithVertexDisplacement(meshFluid, referenceVertices, uh, aleDisp);

      pOld = model.getState().par;
      pCur = model.getState().par;
      cfg.pressureRampTime = 0.0;
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
      // Smooth cosine ramp of the imposed pressures over the startup window so
      // the wall is loaded gradually instead of being step-loaded in one dt.
      constexpr Real kPi = 3.14159265358979323846;
      const Real ramp =
          (cfg.pressureRampTime > 0.0)
              ? 0.5 * (1.0 - std::cos(kPi * std::min(s.t / cfg.pressureRampTime,
                                                     1.0)))
              : 1.0;
      pinValue = ramp * s.par;
      for (const auto &[tag, bc] : wk) {
        const Real outlet = ramp * bc.pout;
        // Compress the inlet->outlet drop so the driving gradient is
        // pressureDropScale * original (default 0.5 -> ~500 Pa instead of
        // ~1000 Pa).  pressureDropScale = 1.0 recovers the raw 0D pressures.
        outletPressureValue[tag] =
            pinValue - cfg.pressureDropScale * (pinValue - outlet);
      }

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

          // Under-relaxed interface displacement update + convergence test.
          auto delta = dState;
          delta -= dIter;
          PetscReal deltaNorm = 0.0;
          PetscReal stateNorm = 0.0;
          VecNorm(delta.getData(), NORM_2, &deltaNorm);
          delta *= cfg.couplingRelaxation;
          dIter += delta;
          VecNorm(dState.getData(), NORM_2, &stateNorm);
          const Real rel = (stateNorm > 0.0)
                               ? (static_cast<Real>(deltaNorm) /
                                  static_cast<Real>(stateNorm))
                               : static_cast<Real>(deltaNorm);

          solidAcceleration = dIter;
          solidAcceleration -= dPred;
          solidAcceleration *= 1.0 / (betaN * dt * dt);
          solidVelocity = vPred;
          tmp = solidAcceleration;
          tmp *= gammaN * dt;
          solidVelocity += tmp;

          if (cfg.couplingIterations > 1 && isRoot) {
            Alert::Info() << "  coupling iterate " << couple << " / "
                          << cfg.couplingIterations
                          << "  relative interface change = " << rel
                          << Alert::Raise;
          }

          if (couple > 1 && rel < cfg.couplingTolerance)
            break;

          restoreMeshToReference(meshFluid, referenceVertices);
          ale.assemble();
          Solver::KSP(ale).solve();
          aleDisp.setData(dMove.getSolution().getData());

          // ALE mesh velocity.
          meshVelocity = aleDisp;
          meshVelocity -= aleDispOld;
          meshVelocity *= 1.0 / dt;

          moveMeshWithVertexDisplacement(meshFluid, referenceVertices, uh,
                                         aleDisp);

          // VMS projections on the current configuration (order matters: tau
          // and the subscale read the projected convection / viscosity).
          muProj.assemble();
          Solver::KSP(muProj).solve();
          l2ConvU.assemble();
          Solver::KSP(l2ConvU).solve();
          tauProj.assemble();
          Solver::KSP(tauProj).solve();
          subProj.assemble();
          Solver::KSP(subProj).solve();

          flow.assemble().setFieldSplits();

          // Inject (rho/dt)(u^n, v)|_{Omega^n} into the velocity block of RHS.
          // massOld's vector ordering coincides with the velocity-block
          // ordering of the monolithic RHS (global index = vOff + i), which is
          // exact in this single-process build (VecGetOwnershipRange spans the
          // whole vector).
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

      // RCR / Windkessel update from the converged interface flux.
      for (const Attribute outlet : BoundaryFluid::Outlets)
        updateRCRNonNew(cfg, outlet, model, wk[outlet], qOut[outlet], dt);

      // Commit the step-(n+1) state.
      uOld.setData(uCur.getData());
      pOld.setData(pCur.getData());
      dOld.setData(dState.getData());
      solidVelocityOld.setData(solidVelocity.getData());
      solidAccelerationOld.setData(solidAcceleration.getData());
      // VMS dynamic-subscale history u'^n <- u'^{n+1}.
      subOld.setData(subTrial.getSolution().getData());


      aleDispOld.setData(aleDisp.getData());


      xdmf_fluid.write(s.t).flush();
      xdmf_solid.write(s.t).flush();

      if (isRoot && csv) {
        csv << s.t << ',' << s.y << ',' << s.v << ',' << s.pv << ',' << s.par
            << ',' << s.pd << ',' << qIn << ',' << qOutSum;
        for (const Attribute outlet : BoundaryFluid::Outlets)
          csv << ',' << qOut[outlet] << ',' << wk[outlet].pout;
        csv << '\n';
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
