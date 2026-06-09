/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * Explicit (staggered) PETSc/MPI coronary ALE FSI prototype.
 *
 * This file is the partitioned/explicit counterpart of
 * examples/Heart/CoronaryArtery_FSI_PETSc_MPI.cpp (monolithic). It keeps the
 * exact same physics (coronary Carreau-Yasuda fluid, NeoHookean/Newmark solid,
 * harmonic ALE, 0D CCMLC2014 heart model + RCR outlets) but replaces the single
 * monolithic SNES by a Dirichlet-Neumann staggered coupling, following the
 * proven explicit scheme of examples/PETSc/PDEs/Seq_BDF1_ALE_FSI.cpp.
 *
 * Mesh assumption (identical to the monolithic example):
 *   A single conforming FSI mesh with cell attributes
 *
 *     1: fluid
 *     2: solid
 *
 *   and boundary/interface attributes
 *
 *     2: coronary FSI wall/interface   (INTERNAL interface between cells 1 and 2)
 *     3: inlet
 *     4..9: outlets
 *     100: solid ring/support
 *
 * Staggered Dirichlet-Neumann step n -> n+1:
 *
 *   1. Advance the 0D heart model; read inlet/outlet pressures.
 *   2. Newmark predictors dPred, vPred for the solid displacement/velocity.
 *   3. (coupling sub-iteration, default 1 = loosely coupled / pure explicit)
 *      a. Harmonic ALE extension of the current solid-displacement iterate over
 *         the fluid: aleDisp = harmonic(dIter). This yields the FULL mesh-motion
 *         field (harmonic in the fluid, equal to dIter in the solid).
 *      b. meshVelocity = (aleDisp - aleDispOld) / dt   (GENUINE, nonzero).
 *      c. Move the (shared) mesh to the aleDisp configuration.
 *      d. FLUID: linear Oseen ALE solve for (u, p). The fluid receives the mesh
 *         velocity as a Dirichlet condition on the FSI interface
 *         (no-slip + Geometric Conservation Law consistent), and natural
 *         pressure on inlet/outlets.
 *      e. Capture the fluid Cauchy traction on the FSI interface as a dead field
 *         t_f = -(p I - 2 mu sym(grad u)) . n_f, evaluated from the FLUID side
 *         (traceOf) on the CURRENT mesh.
 *      f. Restore the mesh to the reference configuration.
 *      g. SOLID: nonlinear NeoHookean/Newmark SNES solve on the REFERENCE
 *         configuration (total Lagrangian) with the captured traction applied as
 *         an explicit Neumann load on the FSI interface.
 *      h. Optional Aitken-free constant under-relaxation of the displacement and
 *         convergence check on the interface displacement increment.
 *   4. Re-solve the harmonic ALE with the final solid displacement, move the
 *      mesh to the consistent output geometry, store aleDispOld for the next
 *      step, and commit (uOld, pOld, dOld, solidVelocityOld, ...).
 *   5. Compute inlet/outlet flow rates and update the RCR/Windkessel models.
 *
 * Why a separate file:
 *   The monolithic version enforces traction continuity automatically. The
 *   staggered version must (i) transfer the fluid traction explicitly, (ii)
 *   transfer the solid velocity explicitly as a Dirichlet condition for the
 *   fluid, and (iii) assemble the solid on the reference configuration while the
 *   fluid is assembled on the moved configuration. All three are handled here.
 *
 * IMPORTANT - added-mass instability:
 *   Loosely coupled (Dirichlet-Neumann, 1 sub-iteration) FSI is only
 *   conditionally stable and becomes unstable when the fluid and solid densities
 *   are comparable on slender domains. The coronary case has
 *   fluidDensity == solidDensity == 1060, which is the worst case. For
 *   production use either the monolithic example, or increase
 *   -coronary_coupling_iterations (strong coupling) with
 *   -coronary_coupling_relaxation < 1, and/or reduce -coronary_traction_scale.
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
#include <map>
#include <numbers>
#include <unordered_map>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <Rodin/Alert.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/MPI.h>
#include <Rodin/MPI/Context/MPI.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/PETSc.h>
#include <Rodin/Solid.h>
#include <Rodin/Solver.h>
#include <Rodin/Variational.h>

#ifdef RODIN_USE_SCOTCH
#include <Rodin/Scotch/MeshPartitioner.h>
#endif

#include "Rodin/Heart/CCMLC2014.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Math;
using namespace Rodin::Solver;
using namespace Rodin::Variational;

namespace
{
  using Model = Rodin::Heart::CCMLC2014T<>;
  using MeshType = Geometry::Mesh<Context::MPI>;

  constexpr int RootRank = 0;

  struct Volume
  {
    static constexpr Attribute Fluid = 1;
    static constexpr Attribute Solid = 2;
  };

  struct Boundary
  {
    static constexpr Attribute FSI = 2;
    static constexpr Attribute Inlet = 3;
    static constexpr std::array<Attribute, 6> Outlets{{4, 5, 6, 7, 8, 9}};
    static constexpr Attribute SolidRing = 100;
  };

  struct RCR
  {
    Real Rp = 1.0e8;
    Real C = 1.0e-11;
    Real Rd = 1.0e9;
    Real pd = 500.0;
    Real pc = 10500.0;
    Real pout = 11000.0;
    Real qd = 0.0;
  };

  struct CarreauYasuda
  {
    Real mu0 = 0.0186058;
    Real muInf = 0.0042963;
    Real lambda = 0.2435;
    Real n = 0.2079;
    Real yasuda = 1.5410;
    Real gammaRegularization = 1.0e-3;
  };

  struct Config
  {
    std::string xdmfBasename = "CoronaryArtery_FSI_Explicit";
    std::string csvPath = "CoronaryArtery_FSI_Explicit.csv";
    Real meshScale = 1.0e-3;

    Real dt = 1.0e-3;
    size_t nsteps = 3 * static_cast<int>(0.85 / 1.0e-3);

    Real fluidDensity = 1060.0;
    Real pressurePenalty = 1.0e-12;
    Real inletBackflowStabilization = 1.0;
    Real outletBackflowStabilization = 1.0;
    CarreauYasuda viscosity;

    Real solidDensity = 1060.0;
    Real solidYoungModulus = 5.0e5;
    Real solidPoissonRatio = 0.3;
    Real newmarkBeta = 0.25;
    Real newmarkGamma = 0.5;

    Real ringStiffness = 8.0e4;
    Real ringDamping = 4.0e4;
    Real inactiveRegularization = 1.0e-10;

    // Explicit-coupling specific options.
    Real tractionScale = 1.0;          // overall scale of the transferred traction
    Real tractionPressureScale = 1.0;  // scale of the pressure part (Seq used 0.05)
    Real tractionViscousScale = 1.0;   // scale of the viscous part
    size_t couplingIterations = 1;     // 1 = loosely coupled; >1 = strong coupling
    Real couplingRelaxation = 1.0;     // displacement under-relaxation for k>1
    Real couplingTolerance = 1.0e-6;   // relative interface-displacement tolerance
    Real robinAlpha = 1e4; // Value Robin Parameter for coupling
  };

  static Real periodic_activation(Real t)
  {
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

    Real alpha = 0.0;
    Real value = minValue;
    if (tau < t1)
    {
      alpha = -(tau - t1) / t1;
      value = alpha * minValue + (1.0 - alpha) * maxValue;
    }
    else if (tau < t2)
    {
      value = maxValue;
    }
    else if (tau < t3)
    {
      alpha = -(tau - t3) / (t3 - t2);
      value = alpha * maxValue + (1.0 - alpha) * minValue;
    }
    else if (tau < t4)
    {
      alpha = -(tau - t4) / (t4 - t3);
      value = alpha * minValue + (1.0 - alpha) * secondThreshold;
    }
    else if (tau < t5)
    {
      value = secondThreshold;
    }
    else if (tau < t6)
    {
      alpha = -(tau - t6) / (t6 - t5);
      value = alpha * secondThreshold + (1.0 - alpha) * minValue;
    }
    return value;
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

  static void updateRCR(const Model& model, RCR& bc, Real Q, Real dt)
  {
    const Real cap = bc.C / dt;
    const auto& s = model.getState();

    bc.pc = (cap * bc.pc + Q + s.pv / bc.Rd) / (cap + 1.0 / bc.Rd);
    bc.qd = (bc.pc - bc.pd) / bc.Rd;
    bc.pout = bc.pc + bc.Rp * Q;
  }

  static MeshType makeMesh(const Context::MPI& context, const Config& cfg, const std::string& meshPath)
  {
    const auto& comm = context.getCommunicator();

    Rodin::MPI::Sharder sharder(context);
    if (comm.rank() == RootRank)
    {
      Geometry::Mesh<Context::Local> local;
      local.load(meshPath, IO::FileFormat::MEDIT);

      if (local.getSpaceDimension() != 3)
        throw std::runtime_error("Expected a 3D coronary FSI mesh.");

      const size_t D = local.getDimension();
      local.getConnectivity().compute(D, D);
      local.getConnectivity().compute(D, 0);
      local.getConnectivity().compute(D, D - 1);
      local.getConnectivity().compute(D - 1, D);
      local.getConnectivity().compute(D - 1, 0);
      local.getConnectivity().compute(D - 1, 1);
      local.getConnectivity().compute(1, 0);

#ifdef RODIN_USE_SCOTCH
      Scotch::Partitioner partitioner(local);
#else
      Geometry::BalancedCompactPartitioner partitioner(local);
#endif
      partitioner.partition(static_cast<size_t>(comm.size()));
      sharder.shard(partitioner);
      sharder.scatter(RootRank);
    }

    MeshType mesh = sharder.gather(RootRank);
    mesh.scale(cfg.meshScale);

    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);
    mesh.getConnectivity().compute(D, 0);
    mesh.getConnectivity().compute(D, D - 1);
    mesh.getConnectivity().compute(D - 1, D);
    mesh.getConnectivity().compute(D - 1, 0);
    mesh.getConnectivity().compute(D - 1, 1);
    mesh.getConnectivity().compute(1, 0);
    mesh.reconcile(2);
    mesh.reconcile(1);

    return mesh;
  }

  template <class MeshT>
  static void saveReferenceVertices(
      const MeshT& mesh,
      std::vector<Math::SpatialPoint>& vertices)
  {
    vertices.resize(mesh.getVertexCount());
    for (auto it = mesh.getVertex(); it; ++it)
      vertices[it->getIndex()] = mesh.getVertexCoordinates(it->getIndex());
  }

  template <class MeshT>
  static void restoreMeshToReference(
      MeshT& mesh,
      const std::vector<Math::SpatialPoint>& referenceVertices)
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
  static void moveMeshWithVertexDisplacement(
      MeshT& mesh,
      const std::vector<Math::SpatialPoint>& referenceVertices,
      const FESType& displacementFES,
      const GridFunctionType& displacement)
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
  //
  // Both meshes are loaded from the SAME mesh file and partitioned with the SAME
  // deterministic partitioner, so every FSI face on one mesh has a geometric
  // twin on the other with identical centroid coordinates.  The map is therefore
  // built PER RANK: each process matches only its own local FSI faces.
  //
  // The original 2D-strip helper keyed the match on the x-coordinate of the face
  // centroid alone.  That is ambiguous in 3D (many faces share an x) and is not
  // MPI-safe.  Here the key combines ALL spatial coordinates of the centroid, so
  // the fluid<->solid pairing is unambiguous in space on every rank.
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

  // Quantize a centroid into an integer key on ALL spatial coordinates so that
  // fluid and solid FSI faces are matched by their true location in space.
  static std::array<long long, 3> centroidKey(const Math::SpatialPoint& c)
  {
    constexpr Real scale = 1e9;
    std::array<long long, 3> k{ 0, 0, 0 };
    const Index n =
      std::min<Index>(3, static_cast<Index>(c.size()));
    for (Index i = 0; i < n; ++i)
      k[i] = static_cast<long long>(std::llround(scale * c(i)));
    return k;
  }

  static InterfaceMap buildInterfaceMap(
      const MeshType& fluidReferenceMesh,
      const MeshType& solidReferenceMesh)
  {
    InterfaceMap map;
    std::map<std::array<long long, 3>, Index> solidByCentroid;

    // Index every LOCAL solid FSI face by its full (3D) centroid key.
    for (auto it = solidReferenceMesh.getBoundary(); it; ++it)
    {
      if (it->getAttribute() != Boundary::FSI)
        continue;

      solidByCentroid.emplace(
          centroidKey(centroid(solidReferenceMesh, *it)),
          it->getIndex());
    }

    // Match every LOCAL fluid FSI face to its solid twin on the same rank.
    for (auto it = fluidReferenceMesh.getBoundary(); it; ++it)
    {
      if (it->getAttribute() != Boundary::FSI)
        continue;

      const auto found =
        solidByCentroid.find(centroidKey(centroid(fluidReferenceMesh, *it)));

      if (found == solidByCentroid.end())
        throw std::runtime_error(
            "Failed to match a fluid FSI face to a solid FSI face on this rank.");

      const Index fluidFace = it->getIndex();
      const Index solidFace = found->second;

      map.fluidToSolid.emplace(fluidFace, solidFace);
      map.solidToFluid.emplace(solidFace, fluidFace);
      map.segments.push_back({ fluidFace, solidFace });
    }

    return map;
  }

  [[maybe_unused]] static Point forwardFluidPointToSolid(
      const Point& p,
      const MeshType& solidReferenceMesh,
      const InterfaceMap& map)
  {
    const auto found = map.fluidToSolid.find(p.getPolytope().getIndex());

    if (found == map.fluidToSolid.end())
      throw std::runtime_error("Fluid point is not on a mapped FSI face.");

    auto solidFace = solidReferenceMesh.getFace(found->second);

    Point q(p);
    q.setPolytope(*solidFace);

    return q;
  }

  static Point forwardSolidPointToFluid(
      const Point& p,
      const MeshType& fluidMesh,
      const InterfaceMap& map)
  {
    const auto found = map.solidToFluid.find(p.getPolytope().getIndex());

    if (found == map.solidToFluid.end())
      throw std::runtime_error("Solid point is not on a mapped FSI face.");

    auto fluidFace = fluidMesh.getFace(found->second);

    Point q(p);
    q.setPolytope(*fluidFace);

    return q;
  }

  static std::string lower(std::string s)
  {
    std::transform(
        s.begin(), s.end(), s.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
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
    PetscOptionsGetString(
        PETSC_NULLPTR, PETSC_NULLPTR,
        "-coronary_fsi_mesh",
        meshPath,
        sizeof(meshPath),
        &meshPathSet);
    if (meshPathSet)
      cfg.meshPath = meshPath;

    PetscReal dt = cfg.dt;
    PetscBool dtSet = PETSC_FALSE;
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_dt", &dt, &dtSet);
    if (dtSet)
      cfg.dt = dt;

    PetscInt nsteps = static_cast<PetscInt>(cfg.nsteps);
    PetscBool nstepsSet = PETSC_FALSE;
    PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_nsteps", &nsteps, &nstepsSet);
    if (nstepsSet)
      cfg.nsteps = static_cast<size_t>(std::max<PetscInt>(0, nsteps));

    PetscReal tractionScale = cfg.tractionScale;
    PetscBool tractionScaleSet = PETSC_FALSE;
    PetscOptionsGetReal(
        PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_traction_scale",
        &tractionScale, &tractionScaleSet);
    if (tractionScaleSet)
      cfg.tractionScale = tractionScale;

    PetscReal tractionPressureScale = cfg.tractionPressureScale;
    PetscBool tractionPressureScaleSet = PETSC_FALSE;
    PetscOptionsGetReal(
        PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_traction_pressure_scale",
        &tractionPressureScale, &tractionPressureScaleSet);
    if (tractionPressureScaleSet)
      cfg.tractionPressureScale = tractionPressureScale;

    PetscReal tractionViscousScale = cfg.tractionViscousScale;
    PetscBool tractionViscousScaleSet = PETSC_FALSE;
    PetscOptionsGetReal(
        PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_traction_viscous_scale",
        &tractionViscousScale, &tractionViscousScaleSet);
    if (tractionViscousScaleSet)
      cfg.tractionViscousScale = tractionViscousScale;

    PetscInt couplingIterations = static_cast<PetscInt>(cfg.couplingIterations);
    PetscBool couplingIterationsSet = PETSC_FALSE;
    PetscOptionsGetInt(
        PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_coupling_iterations",
        &couplingIterations, &couplingIterationsSet);
    if (couplingIterationsSet)
      cfg.couplingIterations =
        static_cast<size_t>(std::max<PetscInt>(1, couplingIterations));

    PetscReal couplingRelaxation = cfg.couplingRelaxation;
    PetscBool couplingRelaxationSet = PETSC_FALSE;
    PetscOptionsGetReal(
        PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_coupling_relaxation",
        &couplingRelaxation, &couplingRelaxationSet);
    if (couplingRelaxationSet)
      cfg.couplingRelaxation = couplingRelaxation;

    PetscReal couplingTolerance = cfg.couplingTolerance;
    PetscBool couplingToleranceSet = PETSC_FALSE;
    PetscOptionsGetReal(
        PETSC_NULLPTR, PETSC_NULLPTR, "-coronary_coupling_tolerance",
        &couplingTolerance, &couplingToleranceSet);
    if (couplingToleranceSet)
      cfg.couplingTolerance = couplingTolerance;
  }
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  // Robust direct sub-solvers by default (overridable on the command line).
  setPETScDefault("-ksp_type", "preonly");
  setPETScDefault("-pc_type", "lu");
  setPETScDefault("-pc_factor_mat_solver_type", "mumps");

  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world(PETSC_COMM_WORLD, boost::mpi::comm_attach);
  Rodin::Context::MPI context(env, world);
  const auto& comm = context.getCommunicator();
  const bool isRoot = comm.rank() == RootRank;

  try
  {
    Config cfg;
    readOptions(cfg);

    Model::Input modelInput = makeModelInput();
    Model model(modelInput);
    initializeModel(model, modelInput);

    const std::string fluidMesh = "../resources/examples/Heart/CoronaryArtery_FSI.medit.mesh";
    MeshType meshFluid = makeMesh(context, cfg, fluidMesh);
    const size_t dimFluid = meshFluid.getSpaceDimension();

    const std::string solidMesh = "../resources/examples/Heart/CoronaryArtery_FSI.medit.mesh";
    MeshType meshSolid= makeMesh(context, cfg, solidMesh);
    const size_t dimSolid = meshSolid.getSpaceDimension();

    // Spatial dimension shared by both blocks (used by the coupling fields).
    const size_t dim = dimFluid;

    std::vector<Math::SpatialPoint> referenceVertices;
    saveReferenceVertices(meshFluid, referenceVertices);

    const InterfaceMap interfaceMap = buildInterfaceMap(meshFluid, meshSolid);

    using VelocityFES =
      H1<2, Math::SpatialVector<Real>, MeshType>;
    using PressureFES =
      H1<1, Real, MeshType>;
    using DisplacementFES =
      H1<1, Math::SpatialVector<Real>, MeshType>;
    // ALE mesh-motion space lives on the fluid mesh; same structure as the solid
    // displacement space but a distinct (vector) field on meshFluid.
    using DisplacementFluidFES =
      H1<1, Math::SpatialVector<Real>, MeshType>;

    VelocityFES uh(std::integral_constant<size_t, 2>{}, meshFluid, dimFluid);
    PressureFES ph(std::integral_constant<size_t, 1>{}, meshFluid);
    DisplacementFluidFES dfh(std::integral_constant<size_t, 1>{}, meshFluid, dimFluid);

    DisplacementFES dh(std::integral_constant<size_t, 1>{}, meshSolid, dimSolid);

    // Fluid trial/test (velocity-pressure).
    PETSc::Variational::TrialFunction u(uh);
    PETSc::Variational::TrialFunction p(ph);
    PETSc::Variational::TestFunction  v(uh);
    PETSc::Variational::TestFunction  q(ph);

    // Harmonic ALE trial/test (full mesh-motion field).
    PETSc::Variational::TrialFunction dAle(dfh);
    PETSc::Variational::TestFunction  vAle(dfh);

    // Solid (displacement increment) trial/test.
    PETSc::Variational::TrialFunction d(dh);
    PETSc::Variational::TestFunction  w(dh);

    u.setName("u");
    p.setName("p");
    dAle.setName("dAle");
    d.setName("disp");

    // Fluid state (uOld/pOld hold the most recent fluid solution).
    PETSc::Variational::GridFunction uOld(uh);
    PETSc::Variational::GridFunction pOld(ph);
    PETSc::Variational::GridFunction one(ph);

    // ALE problem
    PETSc::Variational::GridFunction aleDisp(dfh);     // full mesh displacement (this step)
    PETSc::Variational::GridFunction aleDispOld(dfh);  // full mesh displacement (previous step)
    PETSc::Variational::GridFunction meshVelocity(dfh);// (aleDisp - aleDispOld) / dt

    // Solid / mesh-motion state.
    PETSc::Variational::GridFunction dState(dh);  // total solid displacement (current)
    PETSc::Variational::GridFunction dOld(dh);    // total solid displacement (previous step)
    PETSc::Variational::GridFunction dIter(dh);   // relaxed coupling iterate
    PETSc::Variational::GridFunction etaState(dh);// solid displacement increment (SNES state)
    PETSc::Variational::GridFunction dPred(dh);   // Newmark displacement predictor
    PETSc::Variational::GridFunction vPred(dh);   // Newmark velocity predictor
    PETSc::Variational::GridFunction solidVelocity(dh);
    PETSc::Variational::GridFunction solidAcceleration(dh);
    PETSc::Variational::GridFunction solidVelocityOld(dh);
    PETSc::Variational::GridFunction solidAccelerationOld(dh);
    PETSc::Variational::GridFunction solidFluidTraction(dh); // projected FSI traction

    auto zero = VectorFunction(dim, [&](const Point&)
    {
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
    solidFluidTraction = zero;

    uOld.setName("FluidVelocity");
    pOld.setName("FluidPressure");
    dState.setName("Displacement");
    solidVelocity.setName("SolidVelocity");
    meshVelocity.setName("ALEMeshVelocity");
    solidFluidTraction.setName("FluidTraction");

    IO::XDMF xdmf(comm, cfg.xdmfBasename, RootRank);
    xdmf.setMesh(meshFluid);
    xdmf.add("FluidVelocity", uOld);
    xdmf.add("FluidPressure", pOld);
    xdmf.add("Displacement", dState);
    xdmf.add("SolidVelocity", solidVelocity);
    xdmf.add("ALEMeshVelocity", meshVelocity);
    xdmf.add("FluidTraction", solidFluidTraction);
    xdmf.write(0.0).flush();

    std::ofstream csv;
    if (isRoot)
    {
      csv.open(cfg.csvPath);
      csv << "t,lv_y,lv_v,lv_pv,lv_par,lv_pd,q_in,q_out_total";
      for (const Attribute outlet : Boundary::Outlets)
        csv << ",q_out_" << outlet << ",p_out_" << outlet;
      csv << '\n';
    }

    std::map<Attribute, RCR> wk;
    for (const Attribute outlet : Boundary::Outlets)
      wk.emplace(outlet, RCR{});

    Real pinValue = model.getState().par;
    std::map<Attribute, Real> outletPressureValue;
    for (const Attribute outlet : Boundary::Outlets)
      outletPressureValue[outlet] = wk.at(outlet).pout;

    auto pin = RealFunction([&](const Point&) { return pinValue; });
    auto pout0 = RealFunction([&](const Point&) { return outletPressureValue[Boundary::Outlets[0]]; });
    auto pout1 = RealFunction([&](const Point&) { return outletPressureValue[Boundary::Outlets[1]]; });
    auto pout2 = RealFunction([&](const Point&) { return outletPressureValue[Boundary::Outlets[2]]; });
    auto pout3 = RealFunction([&](const Point&) { return outletPressureValue[Boundary::Outlets[3]]; });
    auto pout4 = RealFunction([&](const Point&) { return outletPressureValue[Boundary::Outlets[4]]; });
    auto pout5 = RealFunction([&](const Point&) { return outletPressureValue[Boundary::Outlets[5]]; });

    // ----------------------------------------------------------------------
    // Constitutive and time-integration constants.
    // ----------------------------------------------------------------------
    const Real solidLambda =
      cfg.solidYoungModulus * cfg.solidPoissonRatio
      / ((1.0 + cfg.solidPoissonRatio) * (1.0 - 2.0 * cfg.solidPoissonRatio));
    const Real solidMu =
      cfg.solidYoungModulus / (2.0 * (1.0 + cfg.solidPoissonRatio));

    const Real dt = cfg.dt;
    const Real betaN = cfg.newmarkBeta;
    const Real gammaN = cfg.newmarkGamma;
    const Real solidMass = cfg.solidDensity / (betaN * dt * dt);
    const Real solidVelocityCoeff = gammaN / (betaN * dt);
    const Real ringVelocityCoeff = cfg.ringDamping * solidVelocityCoeff;

    // Robin coefficient that multiplies the (implicit) Newmark interface
    // velocity on the solid side:  alpha * d(d_s)/dt with
    // d(d_s)/dt = vPred + solidVelocityCoeff * (dState - dPred).
    const Real robinVelocityCoeff = cfg.robinAlpha * solidVelocityCoeff;

    const Real tractionPScale = cfg.tractionScale * cfg.tractionPressureScale;
    const Real tractionVScale = cfg.tractionScale * cfg.tractionViscousScale;

    const auto& cy = cfg.viscosity;
    const Real gammaReg = cy.gammaRegularization;
    const Real deltaMu = cy.mu0 - cy.muInf;

    const auto normal = BoundaryNormal(meshFluid);

    Solid::NeoHookean law(solidLambda, solidMu);
    Solid::MaterialTangent solidTangent(law, d, w, dState);
    Solid::InternalForce   solidInternal(law, w, dState);

    // Current fluid solution (read back after each Oseen solve).
    PETSc::Variational::GridFunction uCur(uh);
    PETSc::Variational::GridFunction pCur(ph);
    uCur = zero;
    pCur = 0.0;

    // Interface velocity imposed on the fluid at the FSI interface (no-slip).
    PETSc::Variational::GridFunction interfaceVelocity(dh);
    interfaceVelocity = zero;

    // ----------------------------------------------------------------------
    // Fluid Oseen / BDF1-ALE problem (linear: transport and viscosity are
    // lagged from the previous fluid solution, so a single KSP solve is
    // enough).  The ALE convective velocity uses the genuine mesh velocity.
    // ----------------------------------------------------------------------
    const auto transportLag = uOld - meshVelocity;
    const auto convU = Mult(Jacobian(u), transportLag);
    const auto divTransportLag = Div(uOld) - Div(meshVelocity);

    const auto symU = 0.5 * (Jacobian(u) + Transpose(Jacobian(u)));
    const auto symV = 0.5 * (Jacobian(v) + Transpose(Jacobian(v)));
    const auto symLag = 0.5 * (Jacobian(uOld) + Transpose(Jacobian(uOld)));
    const auto shearLag = Sqrt(gammaReg * gammaReg + 2.0 * Dot(symLag, symLag));
    const auto muLag =
      cy.muInf
      + deltaMu * Pow(1.0 + Pow(cy.lambda * shearLag, cy.yasuda),
                      (cy.n - 1.0) / cy.yasuda);

    const auto outletBeta = Max(-Dot(transportLag, normal), 0.0);
    const auto inletBeta = Max(Dot(transportLag, normal), 0.0);
    const auto outletBackflow =
      0.5 * cfg.outletBackflowStabilization * cfg.fluidDensity * outletBeta;
    const auto inletBackflow =
      0.5 * cfg.inletBackflowStabilization * cfg.fluidDensity * inletBeta;

    // NOTE: the BDF1 time term has TWO parts that live on DIFFERENT meshes:
    //   (rho/dt) [ \int_{Omega^{n+1}} u . v  -  \int_{Omega^n} u^n . v ].
    // The implicit part (over Omega^{n+1}) stays in 'flow' below.  The explicit
    // u^n part must be integrated on the PREVIOUS configuration Omega^n; it is
    // therefore NOT included here and is assembled separately as 'massOld'
    // (see below) and injected into the fluid RHS at solve time.
    Problem flow(u, p, v, q);
    flow =
        (cfg.fluidDensity / dt) * Integral(u, v)
      + cfg.fluidDensity * Integral(Dot(convU, v))
      + 0.5 * cfg.fluidDensity
        * Integral(divTransportLag * Dot(u, v))
      + 2.0 * Integral(muLag * symU, symV)
      - Integral(p, Div(v))
      + Integral(Div(u), q)
      + cfg.pressurePenalty * Integral(p, q)
      + BoundaryIntegral(inletBackflow * Dot(u, v)).over(Boundary::Inlet)
      + BoundaryIntegral(outletBackflow * Dot(u, v)).over(
          Boundary::Outlets[0], Boundary::Outlets[1], Boundary::Outlets[2],
          Boundary::Outlets[3], Boundary::Outlets[4], Boundary::Outlets[5])
      + BoundaryIntegral(pin * Dot(v, normal)).over(Boundary::Inlet)
      + BoundaryIntegral(pout0 * Dot(v, normal)).over(Boundary::Outlets[0])
      + BoundaryIntegral(pout1 * Dot(v, normal)).over(Boundary::Outlets[1])
      + BoundaryIntegral(pout2 * Dot(v, normal)).over(Boundary::Outlets[2])
      + BoundaryIntegral(pout3 * Dot(v, normal)).over(Boundary::Outlets[3])
      + BoundaryIntegral(pout4 * Dot(v, normal)).over(Boundary::Outlets[4])
      + BoundaryIntegral(pout5 * Dot(v, normal)).over(Boundary::Outlets[5])
      // ----- Robin transmission condition on the FSI interface (fluid side) ---
      //   sigma_f(u, p) n_f + alpha * u_f = alpha * u_s        on Gamma_FSI,
      // where u_s = interfaceVelocity is the (just-computed) solid interface
      // velocity.  The implicit alpha*u_f term goes to the system matrix; the
      // data alpha*u_s goes to the RHS.  This is a Robin condition for the
      // fluid; together with the solid Robin condition below it forms the
      // Robin-Robin scheme.
      + cfg.robinAlpha * BoundaryIntegral(u, v).over(Boundary::FSI)
      - cfg.robinAlpha * BoundaryIntegral(interfaceVelocity, v).over(Boundary::FSI);
    // Standalone load that carries the explicit BDF1 mass term integrated on
    // the PREVIOUS mesh Omega^n:  (rho/dt) (u^n, v).  Assembled each step while
    // the fluid mesh is moved to Omega^n, then added to the fluid RHS.
    PETSc::Variational::TestFunction vMass(uh);
    LinearForm<VelocityFES, ::Vec> massOld(vMass);
    massOld = (cfg.fluidDensity / dt) * Integral(uOld, vMass);

    // ----------------------------------------------------------------------
    // Harmonic ALE mesh-extension problem.
    //   * harmonic in the fluid,
    //   * equal to the solid displacement iterate in the solid block,
    //   * equal to the solid displacement iterate on the FSI interface,
    //   * fixed (zero) at the inlet / outlet rings.
    // ----------------------------------------------------------------------
    Problem ale(dAle, vAle);
    ale =
        Integral(Jacobian(dAle), Jacobian(vAle))
      + DirichletBC(dAle, dIter).on(Boundary::FSI)
      + DirichletBC(dAle, zero).on(
          Boundary::Inlet,
          Boundary::Outlets[0], Boundary::Outlets[1], Boundary::Outlets[2],
          Boundary::Outlets[3], Boundary::Outlets[4], Boundary::Outlets[5]);

    // ----------------------------------------------------------------------
    // Fluid traction on the FSI interface, evaluated from the FLUID side of
    // the (deformed) interface:  t = -sigma_f n_f, with sigma_f the Cauchy
    // stress of the Carreau-Yasuda fluid and n_f the fluid outward normal.
    //   sigma_f = -p I + mu (grad u + grad u^T)
    //    t       =  p n_f - mu (grad u + grad u^T) n_f

    // ----------------------------------------------------------------------
    const auto fsiNormal = FaceNormal(meshFluid).traceOf(Volume::Fluid);
    const auto gradUfsi = Jacobian(uCur).traceOf(Volume::Fluid);
    const auto strainRateFsi = gradUfsi + Transpose(gradUfsi);
    const auto symUfsi = 0.5 * strainRateFsi;
    const auto shearFsi =
      Sqrt(gammaReg * gammaReg + 2.0 * Dot(symUfsi, symUfsi));
    const auto muFsi =
      cy.muInf
      + deltaMu * Pow(1.0 + Pow(cy.lambda * shearFsi, cy.yasuda),
                      (cy.n - 1.0) / cy.yasuda);
    const auto tractionFSI =
        (tractionPScale * pCur) * fsiNormal
      - (tractionVScale * muFsi) * Mult(strainRateFsi, fsiNormal);

    // Fluid interface velocity u_f^n evaluated at a SOLID interface point, used
    // as the lagged Robin data on the solid side.  We map the solid point to its
    // matching fluid face and sample the previous fluid velocity uOld there.
    auto interfaceFluidVelocity = VectorFunction(dim, [&](const Point& xs)
    {
      const Point xf =
        forwardSolidPointToFluid(xs, meshFluid, interfaceMap);

      Math::SpatialVector<Real> value(dim);
      value.setZero();

      const auto uf = uOld(xf);
      for (Index i = 0; i < static_cast<Index>(dim); ++i)
        value(i) = uf(i);

      return value;
    });


    // ----------------------------------------------------------------------
    // Solid Newmark / NeoHookean problem (nonlinear: SNES on the increment
    // etaState, with dState = dOld + etaState).  The transferred fluid
    // traction enters as an explicit Neumann load on the FSI faces.  The
    // increment field is regularized over the fluid block (where it carries
    // no equation of its own) to keep the global matrix nonsingular.
    // ----------------------------------------------------------------------
    Problem solid(d, w);
    solid =
        solidMass * Integral(d, w)
      + solidTangent
      + solidMass * Integral(dState, w)
      - solidMass * Integral(dPred, w)
      + solidInternal
      + (cfg.ringStiffness + ringVelocityCoeff)
        * BoundaryIntegral(d, w).over(Boundary::SolidRing)
      + (cfg.ringStiffness + ringVelocityCoeff)
        * BoundaryIntegral(dState, w).over(Boundary::SolidRing)
      - ringVelocityCoeff
        * BoundaryIntegral(dPred, w).over(Boundary::SolidRing)
      + cfg.ringDamping
        * BoundaryIntegral(vPred, w).over(Boundary::SolidRing)
      // ----- Robin transmission condition on the FSI interface (solid side) ---
      //   sigma_s n_s + alpha * d(d_s)/dt = alpha * u_f^n + t_f^n   on Gamma_FSI,
      // with the Newmark interface velocity
      //   d(d_s)/dt = vPred + solidVelocityCoeff (dState - dPred).
      // Following the same SNES residual/tangent convention used for the ring
      // (a 'd' term feeds the Jacobian, the matching 'dState' term the residual):
      //   alpha * d(d_s)/dt  ->  robinVelocityCoeff (d  and  dState)
      //                          - robinVelocityCoeff dPred + alpha vPred
      //   - alpha u_f^n       ->  lagged fluid interface velocity (RHS data)
      //   - t_f^n             ->  transferred fluid Cauchy traction (RHS data)
      + robinVelocityCoeff * BoundaryIntegral(d, w).over(Boundary::FSI)
      + robinVelocityCoeff * BoundaryIntegral(dState, w).over(Boundary::FSI)
      - robinVelocityCoeff * BoundaryIntegral(dPred, w).over(Boundary::FSI)
      + cfg.robinAlpha     * BoundaryIntegral(vPred, w).over(Boundary::FSI)
      - cfg.robinAlpha     * BoundaryIntegral(interfaceFluidVelocity, w).over(Boundary::FSI)
      - FaceIntegral(solidFluidTraction, w).over(Boundary::FSI);

    solid.assemble();
    Solver::KSP kspSolid(solid);
    Solver::SNES snes(kspSolid);
    snes.setTolerances(1.0e-10, 1.0e-8, 1.0e-10, 50, 10000);
    snes.setStateUpdate([&](const PETSc::Math::Vector& state)
    {
      etaState.setData(state, 0);
      dState = dOld;
      dState += etaState;
    });

    // Interface flux functional q_flux = \int_Gamma (u . n) for the RCR/0D
    // Windkessel coupling at the inlet and outlets.
    using FluxLinearForm = LinearForm<PressureFES, ::Vec>;
    PETSc::Variational::TestFunction qFlux(ph);
    FluxLinearForm flux(qFlux);

    // ======================================================================
    // Time loop: staggered Dirichlet-Neumann FSI.
    //
    //   For each step:
    //     0D heart model advance  ->  inlet / outlet pressures.
    //     Newmark predictors      ->  dPred, vPred.
    //     Coupling iteration k = 1 .. couplingIterations:
    //       (a) harmonic ALE extension from the current interface iterate;
    //       (b) mesh velocity (aleDisp - aleDispOld) / dt;
    //       (c) move the mesh to the ALE configuration;
    //       (d) fluid Oseen solve (interface velocity = solid velocity);
    //       (e) interface flux for the RCR update;
    //       (f) transfer the fluid traction to the FSI faces;
    //       (g) solid Newmark / NeoHookean solve on the reference mesh;
    //       (h) update solid velocity / acceleration;
    //       (i) under-relax the interface displacement and test convergence;
    //       (j) refresh the interface velocity for the next iterate.
    //     Commit the new state, update the RCR, write output.
    //
    //   couplingIterations == 1 reproduces the classical loosely-coupled
    //   (single pass) scheme; couplingIterations > 1 yields strong coupling
    //   (needed when the added-mass effect destabilizes the loose scheme).
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
      pinValue = s.par;
      for (const auto& [tag, bc] : wk)
        outletPressureValue[tag] = bc.pout;

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

      // Initialize the coupling iterate with the predictor and the interface
      // velocity with the previous solid velocity (loose-coupling default).
      dIter = dPred;
      interfaceVelocity = solidVelocityOld;

      if (isRoot)
      {
        Alert::Info()
          << "Coronary explicit ALE FSI step " << step << " / " << cfg.nsteps
          << "  (coupling iterations: " << cfg.couplingIterations << ")"
          << Alert::Raise;
      }

      std::map<Attribute, Real> qOut;
      for (const Attribute outlet : Boundary::Outlets)
        qOut[outlet] = 0.0;
      Real qOutSum = 0.0;
      Real qIn = 0.0;
      bool stepFailed = false;

      for (size_t couple = 1; couple <= cfg.couplingIterations; ++couple)
      {

          solidFluidTraction.project(Region::Faces, tractionFSI, Boundary::FSI);

          // Solid Newmark / NeoHookean solve on the reference mesh.
          snes.solve();
          if (!snes.converged())
          {
            if (isRoot)
              std::cerr << "Solid SNES failed at step " << step
                        << " (coupling iterate " << couple << ") after "
                        << snes.getIterationNumber() << " iterations.\n";
            stepFailed = true;
            break;
          }

          // (h) Solid velocity / acceleration consistent with the new state.
          solidAcceleration = dState;
          solidAcceleration -= dPred;
          solidAcceleration *= 1.0 / (betaN * dt * dt);
          solidVelocity = vPred;
          tmp = solidAcceleration;
          tmp *= gammaN * dt;
          solidVelocity += tmp;

          // (i) Under-relaxed interface displacement update + convergence test.
          auto delta = dState;
          delta -= dIter;
          PetscReal deltaNorm = 0.0;
          PetscReal stateNorm = 0.0;
          VecNorm(delta.getData(), NORM_2, &deltaNorm);
          delta *= cfg.couplingRelaxation;
          dIter += delta;
          VecNorm(dState.getData(), NORM_2, &stateNorm);
          const Real rel =
            (stateNorm > 0.0)
              ? (static_cast<Real>(deltaNorm) / static_cast<Real>(stateNorm))
              : static_cast<Real>(deltaNorm);

        // (a) Harmonic ALE extension on the reference mesh.
        restoreMeshToReference(meshFluid, referenceVertices);
        ale.assemble();
        Solver::KSP(ale).solve();
        aleDisp.setData(dAle.getSolution().getData());

        // (b) Genuine ALE mesh velocity .
        meshVelocity = aleDisp;
        meshVelocity -= aleDispOld;
        meshVelocity *= 1.0 / dt;

        // Solid is solved before the fluid in this staggered order, so the
        // fluid Robin condition uses the just-computed solid interface velocity.
        interfaceVelocity.setData(solidVelocity.getData());

        // Assemble the explicit BDF1 mass term (rho/dt)(u^n, v) on the
        //      PREVIOUS mesh Omega^n.  Move the fluid mesh to the previous ALE
        //      configuration (reference + aleDispOld), assemble massOld there.
        moveMeshWithVertexDisplacement(meshFluid, referenceVertices, dfh, aleDispOld);
        massOld.assemble();

        // (c2) Move the mesh to the current ALE configuration Omega^{n+1}.
        moveMeshWithVertexDisplacement(meshFluid, referenceVertices, dfh, aleDisp);

        // (d) Fluid Oseen solve on the moved (current) mesh.  Assemble all
        //     Omega^{n+1} terms first, then inject the Omega^n mass vector.
        flow.assemble().setFieldSplits();

        // Add (rho/dt)(u^n, v)|_{Omega^n} into the velocity block of the fluid
        // RHS.
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

        // (e) Interface flux for the RCR / Windkessel update (moved mesh).
        flux = BoundaryIntegral(Dot(uCur, normal), qFlux).over(Boundary::Inlet);
        flux.assemble();
        qIn = flux(one);
        qOutSum = 0.0;
        for (const Attribute outlet : Boundary::Outlets)
        {
          flux = BoundaryIntegral(Dot(uCur, normal), qFlux).over(outlet);
          flux.assemble();
          const Real qo = flux(one);
          qOut[outlet] = qo;
          qOutSum += qo;
        }


        if (cfg.couplingIterations > 1 && isRoot)
        {
          Alert::Info()
            << "  coupling iterate " << couple << " / "
            << cfg.couplingIterations
            << "  relative interface change = " << rel
            << Alert::Raise;
        }

        if (rel < cfg.couplingTolerance)
          break;
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

      // RCR / Windkessel update from the converged interface flux.
      for (const Attribute outlet : Boundary::Outlets)
        updateRCR(model, wk[outlet], qOut[outlet], dt);

      // Commit the step-(n+1) state.
      uOld.setData(uCur.getData());
      pOld.setData(pCur.getData());
      dOld.setData(dState.getData());
      solidVelocityOld.setData(solidVelocity.getData());
      solidAccelerationOld.setData(solidAcceleration.getData());


      restoreMeshToReference(meshFluid, referenceVertices);
      ale.assemble();
      Solver::KSP(ale).solve();
      aleDisp.setData(dAle.getSolution().getData());
      meshVelocity = aleDisp;
      meshVelocity -= aleDispOld;
      meshVelocity *= 1.0 / dt;
      aleDispOld.setData(aleDisp.getData());
      moveMeshWithVertexDisplacement(meshFluid, referenceVertices, dfh, aleDisp);

      xdmf.write(s.t).flush();

      if (isRoot && csv)
      {
        csv << s.t << ','
            << s.y << ','
            << s.v << ','
            << s.pv << ','
            << s.par << ','
            << s.pd << ','
            << qIn << ','
            << qOutSum;
        for (const Attribute outlet : Boundary::Outlets)
          csv << ',' << qOut[outlet] << ',' << wk[outlet].pout;
        csv << '\n';
        csv.flush();
      }
    }

    xdmf.close();
    if (isRoot && csv)
      csv.close();
  }
  catch (const std::exception& e)
  {
    if (comm.rank() == RootRank)
      std::cerr << "CoronaryArtery_FSI_Explicit_PETSc_MPI failed: "
                << e.what() << '\n';
    PetscFinalize();
    return 1;
  }

  PetscFinalize();
  return 0;
}
