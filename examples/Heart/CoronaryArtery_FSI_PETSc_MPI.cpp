/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * Monolithic PETSc/MPI coronary ALE FSI prototype.
 *
 * This file follows the monolithic pattern from
 * examples/PETSc/PDEs/Seq_BDF1_ALE_FSI_Monolithic.cpp, but uses the coronary
 * fluid residual from examples/Heart/CoronaryArtery and the nonlinear
 * NeoHookean/Newmark solid residual from CoronaryArterySolid_PETSc_MPI.cpp.
 *
 * Mesh assumption:
 *   The input mesh must be a single conforming FSI mesh with cell attributes
 *
 *     1: fluid
 *     2: solid
 *
 *   and boundary/interface attributes
 *
 *     2: coronary FSI wall/interface
 *     3: inlet
 *     4..9: outlets
 *     100: solid ring/support
 *
 *   The current repository only ships the fluid coronary mesh, so this example
 *   needs a conforming coronary FSI mesh before it can run on the real case.
 *
 * Unknowns in the PETSc SNES state are correction variables:
 *
 *   du     fluid velocity correction
 *   dp     fluid pressure correction
 *   deta   solid/fluid-mesh displacement-increment correction
 *
 * The nonlinear state is updated by SNES through:
 *
 *   uState   = current fluid velocity
 *   pState   = current pressure
 *   etaState = current displacement increment d^{n+1} - d^n
 *   dState   = dOld + etaState
 *
 * The FSI kinematic row is imposed at correction level:
 *
 *   du = cV * deta + (vSolidState - uState) on Gamma_FSI
 *
 * where cV = gamma / (beta dt) is the Newmark derivative of solid velocity
 * with respect to displacement.
 *
 * Flow mode:
 *   - newton: full ALE Navier-Stokes tangent, including current transport
 *     derivatives with respect to du and deta, plus Carreau-Yasuda viscosity
 *     derivative with respect to du.
 *   - oseen: lagged ALE transport and lagged viscosity in the fluid tangent.
 *
 * In both modes the coupled system is solved with SNES because the solid is
 * always nonlinear.
 */

#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

#include <mpi.h>
#include <petscmat.h>
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

  static double beginTimer(MPI_Comm comm, bool isRoot, const std::string& label)
  {
    MPI_Barrier(comm);
    if (isRoot)
      std::cout << "[timer] begin " << label << std::endl;
    return MPI_Wtime();
  }

  static void endTimer(
      MPI_Comm comm,
      bool isRoot,
      const std::string& label,
      double start)
  {
    const double local = MPI_Wtime() - start;
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, comm);

    if (isRoot)
    {
      std::cout << std::fixed << std::setprecision(6)
                << "[timer] end   " << label
                << " : " << global << " s max/rank"
                << std::endl;
    }
  }

  static void printLinearSystemInfo(
      MPI_Comm comm,
      bool isRoot,
      const char* label,
      ::Mat A,
      ::Vec x,
      ::Vec b)
  {
    (void) isRoot;

    const char* xtype = PETSC_NULLPTR;
    const char* btype = PETSC_NULLPTR;
    const char* atype = PETSC_NULLPTR;

    PetscErrorCode ierr = VecGetType(x, &xtype);
    if (ierr != PETSC_SUCCESS)
      xtype = "(VecGetType failed)";

    ierr = VecGetType(b, &btype);
    if (ierr != PETSC_SUCCESS)
      btype = "(VecGetType failed)";

    ierr = MatGetType(A, &atype);
    if (ierr != PETSC_SUCCESS)
      atype = "(MatGetType failed)";

    PetscInt xLocal = 0, xGlobal = 0;
    PetscInt bLocal = 0, bGlobal = 0;
    PetscInt mLocal = 0, nLocal = 0;
    PetscInt mGlobal = 0, nGlobal = 0;

    ierr = VecGetLocalSize(x, &xLocal);
    if (ierr != PETSC_SUCCESS) xLocal = -1;
    ierr = VecGetSize(x, &xGlobal);
    if (ierr != PETSC_SUCCESS) xGlobal = -1;

    ierr = VecGetLocalSize(b, &bLocal);
    if (ierr != PETSC_SUCCESS) bLocal = -1;
    ierr = VecGetSize(b, &bGlobal);
    if (ierr != PETSC_SUCCESS) bGlobal = -1;

    ierr = MatGetLocalSize(A, &mLocal, &nLocal);
    if (ierr != PETSC_SUCCESS) { mLocal = -1; nLocal = -1; }
    ierr = MatGetSize(A, &mGlobal, &nGlobal);
    if (ierr != PETSC_SUCCESS) { mGlobal = -1; nGlobal = -1; }

    PetscPrintf(
        comm,
        "[petsc] %s\n"
        "        A type=%s local=%" PetscInt_FMT "x%" PetscInt_FMT
        " global=%" PetscInt_FMT "x%" PetscInt_FMT "\n"
        "        x type=%s local=%" PetscInt_FMT " global=%" PetscInt_FMT "\n"
        "        b type=%s local=%" PetscInt_FMT " global=%" PetscInt_FMT "\n",
        label,
        atype ? atype : "(null)",
        mLocal, nLocal, mGlobal, nGlobal,
        xtype ? xtype : "(null)",
        xLocal, xGlobal,
        btype ? btype : "(null)",
        bLocal, bGlobal);
  }


  struct Volume
  {
    static constexpr Attribute Fluid = 1;
    static constexpr Attribute Solid = 2;
  };

  struct Boundary
  {
    static constexpr Attribute FSI = 2;
    static constexpr Attribute Inlet = 4;
    static constexpr std::array<Attribute, 6> Outlets{{10, 15, 9, 14, 8, 7}};
    static constexpr Attribute SolidRing = 17;
  };

  enum class FlowMode
  {
    Newton,
    Oseen
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
    std::string meshPath = "malla_merge.mesh";
    std::string xdmfBasename = "CoronaryArtery_FSI";
    std::string csvPath = "CoronaryArtery_FSI.csv";
    Real meshScale = 1.0e-3;

    Real dt = 1.0e-3;
    size_t nsteps = 3 * static_cast<int>(0.85 / 1.0e-3);

    Real fluidDensity = 1060.0;
    Real pressurePenalty = 1.0e-12;
    Real inletBackflowStabilization = 1.0;
    Real outletBackflowStabilization = 1.0;
    FlowMode flowMode = FlowMode::Oseen;
    CarreauYasuda viscosity;

    Real solidDensity = 1060.0;
    Real solidYoungModulus = 5.0e4;
    Real solidPoissonRatio = 0.3;
    Real newmarkBeta = 0.25;
    Real newmarkGamma = 0.5;

    Real ringStiffness = 8.0e4;
    Real ringDamping = 0;
    Real inactiveRegularization = 1.0e-10;
    PetscBool moveMeshDuringNewton = PETSC_FALSE;
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
      Rodin::Heart::CCMLC2014::Model::WindkesselRheology::Quemada;

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

  static MeshType makeMesh(const Context::MPI& context, const Config& cfg)
  {
    const auto& comm = context.getCommunicator();

    Rodin::MPI::Sharder sharder(context);
    if (comm.rank() == RootRank)
    {
      Geometry::Mesh<Context::Local> local;
      local.load(cfg.meshPath, IO::FileFormat::MEDIT);

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
      std::cout << "Using Scotch for mesh partitioning.\n";
      Scotch::Partitioner partitioner(local);
#else
      std::cout << "Using BalancedCompactPartitioner for mesh partitioning.\n";
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

  static std::string lower(std::string s)
  {
    std::transform(
        s.begin(), s.end(), s.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
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

    char flowMode[32] = {};
    PetscBool flowModeSet = PETSC_FALSE;
    PetscOptionsGetString(
        PETSC_NULLPTR, PETSC_NULLPTR,
        "-coronary_flow_mode",
        flowMode,
        sizeof(flowMode),
        &flowModeSet);
    if (flowModeSet)
    {
      const std::string mode = lower(flowMode);
      if (mode == "newton")
        cfg.flowMode = FlowMode::Newton;
      else if (mode == "oseen")
        cfg.flowMode = FlowMode::Oseen;
      else
        throw std::runtime_error("Invalid -coronary_flow_mode. Expected newton or oseen.");
    }

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

    PetscOptionsGetBool(
        PETSC_NULLPTR, PETSC_NULLPTR,
        "-coronary_fsi_move_mesh_during_newton",
        &cfg.moveMeshDuringNewton,
        PETSC_NULLPTR);
  }
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world(PETSC_COMM_WORLD, boost::mpi::comm_attach);
  Rodin::Context::MPI context(env, world);
  const auto& comm = context.getCommunicator();
  const bool isRoot = comm.rank() == RootRank;

  if (isRoot)
    std::cout << "Starting CoronaryArtery_FSI_PETSc_MPI on " << comm.size() << " processes." << std::endl;

  try
  {
    Config cfg;
    readOptions(cfg);

    Model::Input modelInput = makeModelInput();
    Model model(modelInput);
    initializeModel(model, modelInput);

    const double makeMeshStart = beginTimer(PETSC_COMM_WORLD, isRoot, "makeMesh");
    MeshType mesh = makeMesh(context, cfg);
    endTimer(PETSC_COMM_WORLD, isRoot, "makeMesh", makeMeshStart);

    const size_t dim = mesh.getSpaceDimension();

    std::vector<Math::SpatialPoint> referenceVertices;
    {
      const double t0 = beginTimer(PETSC_COMM_WORLD, isRoot, "saveReferenceVertices");
      saveReferenceVertices(mesh, referenceVertices);
      endTimer(PETSC_COMM_WORLD, isRoot, "saveReferenceVertices", t0);
    }

    using VelocityFES =
      H1<2, Math::SpatialVector<Real>, MeshType>;
    using PressureFES =
      H1<1, Real, MeshType>;
    using DisplacementFES =
      H1<1, Math::SpatialVector<Real>, MeshType>;

    const double fesStart = beginTimer(PETSC_COMM_WORLD, isRoot, "construct finite element spaces");
    VelocityFES uh(std::integral_constant<size_t, 2>{}, mesh, dim);
    PressureFES ph(std::integral_constant<size_t, 1>{}, mesh);
    DisplacementFES dh(std::integral_constant<size_t, 1>{}, mesh, dim);
    endTimer(PETSC_COMM_WORLD, isRoot, "construct finite element spaces", fesStart);

    if (isRoot)
    {
      std::cout << "[dofs] velocity uh = " << uh.getSize() << std::endl;
      std::cout << "[dofs] pressure ph = " << ph.getSize() << std::endl;
      std::cout << "[dofs] displacement dh = " << dh.getSize() << std::endl;
      std::cout << "[dofs] monolithic total = "
                << (uh.getSize() + ph.getSize() + dh.getSize())
                << std::endl;
    }

    PETSc::Variational::TrialFunction du(uh);
    PETSc::Variational::TrialFunction dp(ph);
    PETSc::Variational::TrialFunction deta(dh);
    PETSc::Variational::TestFunction v(uh);
    PETSc::Variational::TestFunction q(ph);
    PETSc::Variational::TestFunction z(dh);

    du.setName("du");
    dp.setName("dp");
    deta.setName("deta");

    PETSc::Variational::GridFunction uState(uh);
    PETSc::Variational::GridFunction pState(ph);
    PETSc::Variational::GridFunction etaState(dh);
    PETSc::Variational::GridFunction dState(dh);
    PETSc::Variational::GridFunction uOld(uh);
    PETSc::Variational::GridFunction pOld(ph);
    PETSc::Variational::GridFunction etaOld(dh);
    PETSc::Variational::GridFunction dOld(dh);

    PETSc::Variational::GridFunction solidVelocity(dh);
    PETSc::Variational::GridFunction solidAcceleration(dh);
    PETSc::Variational::GridFunction solidVelocityOld(dh);
    PETSc::Variational::GridFunction solidAccelerationOld(dh);
    PETSc::Variational::GridFunction dPred(dh);
    PETSc::Variational::GridFunction vPred(dh);
    PETSc::Variational::GridFunction meshVelocity(dh);
    PETSc::Variational::GridFunction one(ph);

    auto zero = VectorFunction(dim, [&](const Point&)
    {
      Math::SpatialVector<Real> value(dim);
      value.setZero();
      return value;
    });

    {
      const double t0 = beginTimer(PETSC_COMM_WORLD, isRoot, "initialize grid functions");

    uState = zero;
    pState = 0.0;
    etaState = zero;
    dState = zero;
    uOld = zero;
    pOld = 0.0;
    etaOld = zero;
    dOld = zero;
    solidVelocity = zero;
    solidAcceleration = zero;
    solidVelocityOld = zero;
    solidAccelerationOld = zero;
    dPred = zero;
    vPred = zero;
    meshVelocity = zero;
    one = 1.0;

      endTimer(PETSC_COMM_WORLD, isRoot, "initialize grid functions", t0);
    }

    uState.setName("FluidVelocity");
    pState.setName("FluidPressure");
    etaState.setName("DisplacementIncrement");
    dState.setName("Displacement");
    solidVelocity.setName("SolidVelocity");
    solidAcceleration.setName("SolidAcceleration");
    meshVelocity.setName("ALEMeshVelocity");

    size_t nv = mesh.getVertexCount();
    size_t nc = mesh.getCellCount();
    if (isRoot)
    {
      std::cout << "Mesh loaded with " << nv << " vertices and " << nc << " cells." << std::endl;
    }

    IO::XDMF xdmf(comm, cfg.xdmfBasename, RootRank);
    xdmf.setMesh(mesh);
    xdmf.add("FluidVelocity", uState);
    xdmf.add("FluidPressure", pState);
    xdmf.add("Displacement", dState);
    xdmf.add("SolidVelocity", solidVelocity);
    xdmf.add("ALEMeshVelocity", meshVelocity);
    {
      const double t0 = beginTimer(PETSC_COMM_WORLD, isRoot, "initial XDMF write");
      xdmf.write(0.0).flush();
      endTimer(PETSC_COMM_WORLD, isRoot, "initial XDMF write", t0);
    }

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

    const Real solidLambda =
      cfg.solidYoungModulus * cfg.solidPoissonRatio
      / ((1.0 + cfg.solidPoissonRatio) * (1.0 - 2.0 * cfg.solidPoissonRatio));
    const Real solidMu =
      cfg.solidYoungModulus / (2.0 * (1.0 + cfg.solidPoissonRatio));
    Solid::NeoHookean law(solidLambda, solidMu);
    Solid::MaterialTangent solidTangent(law, deta, z, dState);
    Solid::InternalForce solidInternal(law, z, dState);

    const Real dt = cfg.dt;
    const Real betaN = cfg.newmarkBeta;
    const Real gammaN = cfg.newmarkGamma;
    const Real solidMass = cfg.solidDensity / (betaN * dt * dt);
    const Real solidVelocityCoeff = gammaN / (betaN * dt);
    const Real ringVelocityCoeff = cfg.ringDamping * solidVelocityCoeff;

    const auto normal = BoundaryNormal(mesh);

    const auto meshVelocityState = (1.0 / dt) * etaState;
    const auto meshVelocityLag = (1.0 / dt) * etaOld;
    const auto transportState = uState - meshVelocityState;
    const auto transportLag = uOld - meshVelocityLag;

    const auto gradDU = Jacobian(du);
    const auto gradU = Jacobian(uState);
    const auto gradLag = Jacobian(uOld);

    const auto convState = Mult(gradU, transportState);
    const auto convTangentVelocity = Mult(gradDU, transportState);
    const auto convTangentTransportU = Mult(gradU, du);
    const auto convTangentTransportEta = Mult(gradU, deta);
    const auto oseenConvTangent = Mult(gradDU, transportLag);
    const auto oseenConvState = Mult(gradU, transportLag);

    const auto divTransportState = Div(uState) - (1.0 / dt) * Div(etaState);
    const auto divTransportLag = Div(uOld) - (1.0 / dt) * Div(etaOld);

    const auto outletBeta = Max(-Dot(transportLag, normal), 0.0);
    const auto inletBeta = Max(Dot(transportLag, normal), 0.0);
    const auto outletBackflow =
      0.5 * cfg.outletBackflowStabilization * cfg.fluidDensity * outletBeta;
    const auto inletBackflow =
      0.5 * cfg.inletBackflowStabilization * cfg.fluidDensity * inletBeta;

    const auto symDU = 0.5 * (Jacobian(du) + Transpose(Jacobian(du)));
    const auto symV = 0.5 * (Jacobian(v) + Transpose(Jacobian(v)));
    const auto symU = 0.5 * (Jacobian(uState) + Transpose(Jacobian(uState)));
    const auto symLag = 0.5 * (gradLag + Transpose(gradLag));

    const auto& cy = cfg.viscosity;
    const Real gammaReg = cy.gammaRegularization;
    const Real deltaMu = cy.mu0 - cy.muInf;
    const auto shear =
      Sqrt(gammaReg * gammaReg + 2.0 * Dot(symU, symU));
    const auto shearLag =
      Sqrt(gammaReg * gammaReg + 2.0 * Dot(symLag, symLag));
    const auto carreauBase =
      1.0 + Pow(cy.lambda * shear, cy.yasuda);
    const auto carreauBaseLag =
      1.0 + Pow(cy.lambda * shearLag, cy.yasuda);
    const auto mu =
      cy.muInf + deltaMu * Pow(carreauBase, (cy.n - 1.0) / cy.yasuda);
    const auto muLag =
      cy.muInf + deltaMu * Pow(carreauBaseLag, (cy.n - 1.0) / cy.yasuda);
    const auto dshear =
      2.0 * Dot(symU, symDU) / shear;
    const auto dmu =
      deltaMu
      * (cy.n - 1.0)
      * Pow(carreauBase, (cy.n - 1.0 - cy.yasuda) / cy.yasuda)
      * std::pow(cy.lambda, cy.yasuda)
      * Pow(shear, cy.yasuda - 1.0)
      * dshear;

    const auto solidVelocityState =
      vPred + solidVelocityCoeff * (dState - dPred);

    Problem fsi(du, dp, deta, v, q, z);

    {
      const double t0 = beginTimer(PETSC_COMM_WORLD, isRoot, "build variational FSI form expression");

    if (cfg.flowMode == FlowMode::Newton)
    {
      fsi =
          /* Newton ALE fluid tangent. */
            (cfg.fluidDensity / dt) * Integral(du, v).over(Volume::Fluid)
          + cfg.fluidDensity * Integral(Dot(convTangentVelocity, v)).over(Volume::Fluid)
          + cfg.fluidDensity * Integral(Dot(convTangentTransportU, v)).over(Volume::Fluid)
          - (cfg.fluidDensity / dt) * Integral(Dot(convTangentTransportEta, v)).over(Volume::Fluid)
          + 0.5 * cfg.fluidDensity * Integral(divTransportState * Dot(du, v)).over(Volume::Fluid)
          + 0.5 * cfg.fluidDensity * Integral(Dot(Div(du) * uState, v)).over(Volume::Fluid)
          - (0.5 * cfg.fluidDensity / dt) * Integral(Dot(Div(deta) * uState, v)).over(Volume::Fluid)
          + 2.0 * Integral(mu * symDU, symV).over(Volume::Fluid)
          + 2.0 * Integral(dmu * symU, symV).over(Volume::Fluid)
          - Integral(dp, Div(v)).over(Volume::Fluid)
          + Integral(Div(du), q).over(Volume::Fluid)
          + cfg.pressurePenalty * Integral(dp, q).over(Volume::Fluid)
          + BoundaryIntegral(inletBackflow * Dot(du, v)).over(Boundary::Inlet)
          + BoundaryIntegral(outletBackflow * Dot(du, v)).over(
              Boundary::Outlets[0], Boundary::Outlets[1], Boundary::Outlets[2],
              Boundary::Outlets[3], Boundary::Outlets[4], Boundary::Outlets[5])

          /* Newton ALE fluid residual. */
          + (cfg.fluidDensity / dt) * Integral(uState - uOld, v).over(Volume::Fluid)
          + cfg.fluidDensity * Integral(Dot(convState, v)).over(Volume::Fluid)
          + 0.5 * cfg.fluidDensity * Integral(Dot(divTransportState * uState, v)).over(Volume::Fluid)
          + 2.0 * Integral(mu * symU, symV).over(Volume::Fluid)
          - Integral(pState, Div(v)).over(Volume::Fluid)
          + Integral(Div(uState), q).over(Volume::Fluid)
          + cfg.pressurePenalty * Integral(pState, q).over(Volume::Fluid)
          + BoundaryIntegral(pin * Dot(v, normal)).over(Boundary::Inlet)
          + BoundaryIntegral(pout0 * Dot(v, normal)).over(Boundary::Outlets[0])
          + BoundaryIntegral(pout1 * Dot(v, normal)).over(Boundary::Outlets[1])
          + BoundaryIntegral(pout2 * Dot(v, normal)).over(Boundary::Outlets[2])
          + BoundaryIntegral(pout3 * Dot(v, normal)).over(Boundary::Outlets[3])
          + BoundaryIntegral(pout4 * Dot(v, normal)).over(Boundary::Outlets[4])
          + BoundaryIntegral(pout5 * Dot(v, normal)).over(Boundary::Outlets[5])
          + BoundaryIntegral(inletBackflow * Dot(uState, v)).over(Boundary::Inlet)
          + BoundaryIntegral(outletBackflow * Dot(uState, v)).over(
              Boundary::Outlets[0], Boundary::Outlets[1], Boundary::Outlets[2],
              Boundary::Outlets[3], Boundary::Outlets[4], Boundary::Outlets[5])

          /* Harmonic ALE displacement increment in the fluid. */
          + Integral(Jacobian(deta), Jacobian(z)).over(Volume::Fluid)
          + Integral(Jacobian(etaState), Jacobian(z)).over(Volume::Fluid)

          /* Nonlinear Newmark solid block over the solid cells. */
          + solidMass * Integral(deta, z).over(Volume::Solid)
          + solidTangent.over(Volume::Solid)
          + solidMass * Integral(dState, z).over(Volume::Solid)
          - solidMass * Integral(dPred, z).over(Volume::Solid)
          + solidInternal.over(Volume::Solid)

          /* Coronary solid ring support from CoronaryArterySolid_PETSc_MPI.cpp. */
          + (cfg.ringStiffness + ringVelocityCoeff) *
            BoundaryIntegral(deta, z).over(Boundary::SolidRing)
          + (cfg.ringStiffness + ringVelocityCoeff) *
            BoundaryIntegral(dState, z).over(Boundary::SolidRing)
          - ringVelocityCoeff *
            BoundaryIntegral(dPred, z).over(Boundary::SolidRing)
          + cfg.ringDamping *
            BoundaryIntegral(vPred, z).over(Boundary::SolidRing)

          /* Inactive regularization for globally defined fields. */
          + cfg.inactiveRegularization * Integral(du, v).over(Volume::Solid)
          + cfg.inactiveRegularization * Integral(uState, v).over(Volume::Solid)
          + cfg.inactiveRegularization * Integral(dp, q).over(Volume::Solid)
          + cfg.inactiveRegularization * Integral(pState, q).over(Volume::Solid)

          /* Essential and FSI kinematic constraints. */
          + DirichletBC(deta, -etaState).on(
              Boundary::Inlet,
              Boundary::Outlets[0], Boundary::Outlets[1], Boundary::Outlets[2],
              Boundary::Outlets[3], Boundary::Outlets[4], Boundary::Outlets[5])
          + DirichletBC(
              du,
              solidVelocityCoeff * deta,
              solidVelocityState - uState).on(Boundary::FSI);
    }
    else
    {
const Real fsiPenalty = 1e2; // tune

      fsi =
          /* Oseen ALE fluid tangent with lagged transport and viscosity. */
            (cfg.fluidDensity / dt) * Integral(du, v).over(Volume::Fluid)
          + cfg.fluidDensity * Integral(Dot(oseenConvTangent, v)).over(Volume::Fluid)
          + 0.5 * cfg.fluidDensity * Integral(divTransportLag * Dot(du, v)).over(Volume::Fluid)
          + 2.0 * Integral(muLag * symDU, symV).over(Volume::Fluid)
          - Integral(dp, Div(v)).over(Volume::Fluid)
          + Integral(Div(du), q).over(Volume::Fluid)
          + cfg.pressurePenalty * Integral(dp, q).over(Volume::Fluid)
          + BoundaryIntegral(inletBackflow * Dot(du, v)).over(Boundary::Inlet)
          + BoundaryIntegral(outletBackflow * Dot(du, v)).over(
              Boundary::Outlets[0], Boundary::Outlets[1], Boundary::Outlets[2],
              Boundary::Outlets[3], Boundary::Outlets[4], Boundary::Outlets[5])

          /* Oseen residual, still embedded in global SNES for the solid. */
          + (cfg.fluidDensity / dt) * Integral(uState - uOld, v).over(Volume::Fluid)
          + cfg.fluidDensity * Integral(Dot(oseenConvState, v)).over(Volume::Fluid)
          + 0.5 * cfg.fluidDensity * Integral(Dot(divTransportLag * uState, v)).over(Volume::Fluid)
          + 2.0 * Integral(muLag * symU, symV).over(Volume::Fluid)
          - Integral(pState, Div(v)).over(Volume::Fluid)
          + Integral(Div(uState), q).over(Volume::Fluid)
          + cfg.pressurePenalty * Integral(pState, q).over(Volume::Fluid)
          + BoundaryIntegral(pin * Dot(v, normal)).over(Boundary::Inlet)
          + BoundaryIntegral(pout0 * Dot(v, normal)).over(Boundary::Outlets[0])
          + BoundaryIntegral(pout1 * Dot(v, normal)).over(Boundary::Outlets[1])
          + BoundaryIntegral(pout2 * Dot(v, normal)).over(Boundary::Outlets[2])
          + BoundaryIntegral(pout3 * Dot(v, normal)).over(Boundary::Outlets[3])
          + BoundaryIntegral(pout4 * Dot(v, normal)).over(Boundary::Outlets[4])
          + BoundaryIntegral(pout5 * Dot(v, normal)).over(Boundary::Outlets[5])

          /* Harmonic ALE displacement increment in the fluid. */
          + Integral(Jacobian(deta), Jacobian(z)).over(Volume::Fluid)
          + Integral(Jacobian(etaState), Jacobian(z)).over(Volume::Fluid)

          /* Nonlinear Newmark solid block over the solid cells. */
          + solidMass * Integral(deta, z).over(Volume::Solid)
          + solidTangent.over(Volume::Solid)
          + solidMass * Integral(dState, z).over(Volume::Solid)
          - solidMass * Integral(dPred, z).over(Volume::Solid)
          + solidInternal.over(Volume::Solid)

          /* Coronary solid ring support from CoronaryArterySolid_PETSc_MPI.cpp. */
          + (cfg.ringStiffness + ringVelocityCoeff) *
            BoundaryIntegral(deta, z).over(Boundary::SolidRing)
          + (cfg.ringStiffness + ringVelocityCoeff) *
            BoundaryIntegral(dState, z).over(Boundary::SolidRing)
          - ringVelocityCoeff *
            BoundaryIntegral(dPred, z).over(Boundary::SolidRing)
          + cfg.ringDamping *
            BoundaryIntegral(vPred, z).over(Boundary::SolidRing)

          /* Inactive regularization for globally defined fields. */
          + cfg.inactiveRegularization * Integral(du, v).over(Volume::Solid)
          + cfg.inactiveRegularization * Integral(uState, v).over(Volume::Solid)
          + cfg.inactiveRegularization * Integral(dp, q).over(Volume::Solid)
          + cfg.inactiveRegularization * Integral(pState, q).over(Volume::Solid)

+ fsiPenalty * BoundaryIntegral(Dot(du, v)).over(Boundary::FSI)

- fsiPenalty * solidVelocityCoeff
    * BoundaryIntegral(Dot(deta, v)).over(Boundary::FSI)

- fsiPenalty
    * BoundaryIntegral(Dot(solidVelocityState - uState, v)).over(Boundary::FSI)

- fsiPenalty * solidVelocityCoeff
    * BoundaryIntegral(Dot(du, z)).over(Boundary::FSI)

+ fsiPenalty * solidVelocityCoeff * solidVelocityCoeff
    * BoundaryIntegral(Dot(deta, z)).over(Boundary::FSI)

+ fsiPenalty * solidVelocityCoeff
    * BoundaryIntegral(Dot(solidVelocityState - uState, z)).over(Boundary::FSI)

          /* Essential and FSI kinematic constraints. */
          + DirichletBC(deta, -etaState).on(
              Boundary::Inlet,
              Boundary::Outlets[0], Boundary::Outlets[1], Boundary::Outlets[2],
              Boundary::Outlets[3], Boundary::Outlets[4], Boundary::Outlets[5]);
    }

      endTimer(PETSC_COMM_WORLD, isRoot, "build variational FSI form expression", t0);
    }

    if (isRoot)
      std::cout << "Setting up solver...\n";

    {
      const double t0 = beginTimer(PETSC_COMM_WORLD, isRoot, "initial fsi.assemble");
      fsi.assemble(); // .setFieldSplits();
      endTimer(PETSC_COMM_WORLD, isRoot, "initial fsi.assemble", t0);

      auto& system = fsi.getLinearSystem();
      printLinearSystemInfo(
          PETSC_COMM_WORLD,
          isRoot,
          "after initial fsi.assemble",
          system.getOperator(),
          system.getSolution(),
          system.getVector());
    }

    if (isRoot)
      std::cout << "Configuring SNES solver...\n";

    const double solverConfigStart =
      beginTimer(PETSC_COMM_WORLD, isRoot, "construct/configure KSP and SNES");
    Solver::KSP ksp(fsi);
    Solver::SNES snes(ksp);
    snes.setTolerances(1.0e-10, 1.0e-8, 1.0e-10, 50, 10000);
    endTimer(
        PETSC_COMM_WORLD,
        isRoot,
        "construct/configure KSP and SNES",
        solverConfigStart);

    auto moveMeshToCurrentDisplacement = [&]()
    {
      moveMeshWithVertexDisplacement(mesh, referenceVertices, dh, dState);
    };

    snes.setStateUpdate([&](const PETSc::Math::Vector& state)
    {
      const size_t uOffset = 0;
      const size_t pOffset = uOffset + uh.getSize();
      const size_t etaOffset = pOffset + ph.getSize();

      uState.setData(state, uOffset);
      pState.setData(state, pOffset);
      etaState.setData(state, etaOffset);

      dState = dOld + etaState;
      solidAcceleration = dState;
      solidAcceleration -= dPred;
      solidAcceleration *= 1.0 / (betaN * dt * dt);

      solidVelocity = vPred;
      auto tmp = solidAcceleration;
      tmp *= gammaN * dt;
      solidVelocity += tmp;

      meshVelocity = etaState;
      meshVelocity *= 1.0 / dt;

      if (cfg.moveMeshDuringNewton == PETSC_TRUE)
        moveMeshToCurrentDisplacement();
    });

    using FluxLinearForm =
      LinearForm<PressureFES, ::Vec>;
    PETSc::Variational::TestFunction qFlux(ph);
    FluxLinearForm flux(qFlux);

    for (size_t step = 1; step <= cfg.nsteps; ++step)
    {
      if (isRoot)
        std::cout << "Starting step " << step << " / " << cfg.nsteps << std::endl;

      const double modelStepStart =
        beginTimer(PETSC_COMM_WORLD, isRoot, "0D model.step");
      const auto rep = model.step(dt);
      endTimer(PETSC_COMM_WORLD, isRoot, "0D model.step", modelStepStart);
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

      const double predictorStart =
        beginTimer(PETSC_COMM_WORLD, isRoot, "Newmark predictor/state reset");
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

      uState = uOld;
      pState = pOld;
      etaState = dPred;
      etaState -= dOld;
      dState = dPred;
      meshVelocity = etaState;
      meshVelocity *= 1.0 / dt;
      moveMeshToCurrentDisplacement();
      endTimer(
          PETSC_COMM_WORLD,
          isRoot,
          "Newmark predictor/state reset",
          predictorStart);

      if (isRoot)
      {
        Alert::Info()
          << "Coronary monolithic ALE FSI step " << step
          << " / " << cfg.nsteps
          << "  flow="
          << (cfg.flowMode == FlowMode::Newton ? "newton" : "oseen")
          << Alert::Raise;
      }

      {
        auto& system = fsi.getLinearSystem();
        printLinearSystemInfo(
            PETSC_COMM_WORLD,
            isRoot,
            "before snes.solve",
            system.getOperator(),
            system.getSolution(),
            system.getVector());

        const double t0 = beginTimer(PETSC_COMM_WORLD, isRoot, "snes.solve");
        snes.solve();
        endTimer(PETSC_COMM_WORLD, isRoot, "snes.solve", t0);

        printLinearSystemInfo(
            PETSC_COMM_WORLD,
            isRoot,
            "after snes.solve",
            system.getOperator(),
            system.getSolution(),
            system.getVector());
      }

      if (!snes.converged())
      {
        if (isRoot)
          std::cerr << "SNES failed at step " << step
                    << " after " << snes.getIterationNumber()
                    << " iterations.\n";
        break;
      }

      const double fluxStart =
        beginTimer(PETSC_COMM_WORLD, isRoot, "flux assembly/evaluation");
      std::map<Attribute, Real> qOut;
      Real qOutSum = 0.0;
      Real qIn = 0.0;

      flux = BoundaryIntegral(Dot(uState, normal), qFlux).over(Boundary::Inlet);
      flux.assemble();
      qIn = flux(one);

      for (const Attribute outlet : Boundary::Outlets)
      {
        flux = BoundaryIntegral(Dot(uState, normal), qFlux).over(outlet);
        flux.assemble();
        const Real q = flux(one);
        qOut[outlet] = q;
        qOutSum += q;
      }
      endTimer(PETSC_COMM_WORLD, isRoot, "flux assembly/evaluation", fluxStart);

      const double rcrStart =
        beginTimer(PETSC_COMM_WORLD, isRoot, "RCR update");
      for (const Attribute outlet : Boundary::Outlets)
        updateRCR(model, wk[outlet], qOut[outlet], dt);
      endTimer(PETSC_COMM_WORLD, isRoot, "RCR update", rcrStart);

      const double historyStart =
        beginTimer(PETSC_COMM_WORLD, isRoot, "history update");
      uOld.setData(uState.getData());
      pOld.setData(pState.getData());
      etaOld.setData(etaState.getData());
      dOld.setData(dState.getData());
      solidVelocityOld.setData(solidVelocity.getData());
      solidAccelerationOld.setData(solidAcceleration.getData());
      endTimer(PETSC_COMM_WORLD, isRoot, "history update", historyStart);

      const double outputStart =
        beginTimer(PETSC_COMM_WORLD, isRoot, "move mesh + XDMF write");
      moveMeshToCurrentDisplacement();
      xdmf.write(s.t).flush();
      endTimer(PETSC_COMM_WORLD, isRoot, "move mesh + XDMF write", outputStart);

      if (isRoot && csv)
      {
        const double csvStart = MPI_Wtime();
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
        std::cout << std::fixed << std::setprecision(6)
                  << "[timer] end   CSV write : "
                  << (MPI_Wtime() - csvStart)
                  << " s rank/root"
                  << std::endl;
      }
    }

    xdmf.close();
    if (isRoot && csv)
      csv.close();
  }
  catch (const std::exception& e)
  {
    if (comm.rank() == RootRank)
      std::cerr << "CoronaryArtery_FSI_PETSc_MPI failed: " << e.what() << '\n';
    PetscFinalize();
    return 1;
  }

  PetscFinalize();
  return 0;
}
