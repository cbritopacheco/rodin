/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * Distributed (MPI) PETSc port of the Newmark-beta quasi-compressible
 * NeoHookean coronary solid.
 *
 *   rho * u_tt - div(P(F(u))) = 0
 *
 * CoronaryArterySolid_PETSc.cpp:
 *
 *  1. boost::mpi environment + Rodin::Context::MPI wrap the run; PetscInitialize
 *     is still called (PETSc objects deduce their communicator from the mesh).
 *  2. The MEDIT mesh is loaded on rank 0 as Mesh<Context::Local>, partitioned
 *     (BalancedCompactPartitioner / Scotch), sharded and scattered, then
 *     gathered into a Mesh<Context::MPI> and reconciled for assembly.
 *  3. The 0D model runs redundantly on every rank (deterministic, cheap), so
 *     the 0D->3D coupling value is identical everywhere with no broadcast.
 *     CSV / stdout output is guarded to rank 0.
 *  4. Newmark predictor/corrector vector algebra uses the ghost-aware
 *     GridFunction operators (operator=, +=, -=, *= scalar), NOT raw VecAXPY,
 *     so the ghost layer stays consistent before each (re)assembly.
 *  5. No explicit predictor seeding of the SNES solution vector: in MPI the
 *     ghosted GridFunction vector and the distributed system vector do not
 *     share a layout. SNES warm-starts from the previous step; the predictor
 *     still enters the residual via the -aMass*(uPred, w) term.
 *
 * Run example:
 *   VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
 *   mpirun -n 8 ./CoronaryArterySolid_PETSc_MPI \
 *     -snes_monitor -snes_converged_reason -ksp_converged_reason \
 *     -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps
 */
#include "Rodin/IO/ForwardDecls.h"
#include "Rodin/Variational/ForwardDecls.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <numbers>
#include <stdexcept>

#include <petscsys.h>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <Rodin/PETSc.h>
#include <Rodin/Alert.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/MPI.h>
#include <Rodin/MPI/Context/MPI.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/Solid.h>
#include <Rodin/Solver.h>
#include <Rodin/Variational.h>
#include <type_traits>
#ifdef RODIN_USE_SCOTCH
#include <Rodin/Scotch/MeshPartitioner.h>
#endif
#include "Rodin/Heart/CCMLC2014.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

using Model = Rodin::Heart::CCMLC2014T<>;

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
  // Piecewise-linear approximation of Caruel et al. Fig. 7.
  const Real low_ec = 0.0;
  const Real high_ec = 2.0;
  const Real low_value = 1.6;
  const Real high_value = 1.0;
  if (ec <= low_ec)
    return low_value;
  if (ec >= high_ec)
    return high_value;
  const Real s = (ec - low_ec) / (high_ec - low_ec);
  return (1.0 - s) * low_value + s * high_value;
}
static Real load_dependent_relaxation_dm0(Real ec) {
  const Real low_ec = 0.0;
  const Real high_ec = 2.0;
  const Real low_value = 1.6;
  const Real high_value = 1.0;
  if (ec <= low_ec || ec >= high_ec)
    return 0.0;
  return (high_value - low_value) / (high_ec - low_ec);
}
static Real atrial_pressure(Real t) {
  const Real T = 0.85;
  const Real tau = t - T * std::floor(t / T);
  const Real min_value = 500.0;
  const Real max_value = 1000.0;
  const Real second_threshold = 1250.0;
  const Real t1 = 0.02;
  const Real t2 = 0.15;
  const Real t3 = 0.17;
  const Real t4 = 0.56;
  const Real t5 = 0.62;
  const Real t6 = 0.85;
  Real alpha = 0.0;
  Real value = min_value;
  if (tau < t1) {
    alpha = -(tau - t1) / t1;
    value = alpha * min_value + (1.0 - alpha) * max_value;
  } else if (tau < t2) {
    value = max_value;
  } else if (tau < t3) {
    alpha = -(tau - t3) / (t3 - t2);
    value = alpha * max_value + (1.0 - alpha) * min_value;
  } else if (tau < t4) {
    alpha = -(tau - t4) / (t4 - t3);
    value = alpha * min_value + (1.0 - alpha) * second_threshold;
  } else if (tau < t5) {
    value = second_threshold;
  } else if (tau < t6) {
    alpha = -(tau - t6) / (t6 - t5);
    value = alpha * second_threshold + (1.0 - alpha) * min_value;
  } else {
    value = min_value;
  }
  return value;
}

int main(int argc, char** argv) {

  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  Rodin::Context::MPI context(env, world);
  const auto& comm = context.getCommunicator();
  constexpr int RootRank = 0;
  const bool isRoot = (comm.rank() == RootRank);

  // ----------------------------------------------------------------------
  // Initialization 0D model (maybe to run it in just one rank?).
  // ----------------------------------------------------------------------
  Model::Input in;
  // Geometry / inertia
  in.rho = 1.0e3;
  in.R0 = 2.36e-2;
  in.d0 = 1.42e-2;
  // Active law parameters
  in.Es = 3.0e7;
  in.mu = 70.0;
  in.eta = 70.0;
  in.alpha = 1.5;
  in.alphaR = 0.12;
  in.k0 = 1.0e5;
  in.sigma0 = 1.24e5;
  // Windkessel
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
  // Valve parameters
  in.Kat = 9.0e-6;
  in.Kp = 5.0e-10;
  in.Kar = 1.3e-5;
  in.cavityCapacity = 5.0e-12;
  // Local active solve
  in.localTolerance = 1e-12;
  in.localMaxIterations = 50;
  in.localDamping = 1.0;
  in.absRegularization = 1e-14;
  // Initial local active state
  in.initFibDef = 0.0;
  in.initActiveStiffness = 0.0;
  in.initActiveStress = 0.0;
  in.pSv = [](Real) { return 1.0e3; };
  in.pAt = atrial_pressure;
  in.u = periodic_activation;
  in.m0 = load_dependent_relaxation_m0;
  in.dm0 = load_dependent_relaxation_dm0;
  {
    using PassiveEnergy = std::decay_t<decltype(in.passiveEnergy)>;
    typename PassiveEnergy::Parameters hp;
    hp.mu1 = 0.0;
    hp.mu2 = 0.0;
    hp.C0 = 1.9e3;
    hp.C1 = 1.1e-1;
    hp.C2 = 1.9e3;
    hp.C3 = 1.1e-1;
    in.passiveEnergy = PassiveEnergy(hp);
  }
  Model model(in);
  model.setMaxIterations(200)
      .setAbsoluteTolerance(1e-8)
      .setRelativeTolerance(1e-8)
      .setStepTolerance(1e-10)
      .setDampingFactor(1.0);
  Model::State s0;
  s0.t = 0.0;
  s0.y = 0.0;
  s0.v = 0.0;
  s0.pv = in.pAt(0.0) - 100.0;
  s0.par = 11000.0;
  s0.pd = 10000.0;
  // Initial local active state
  s0.ec = in.initFibDef;
  s0.gamma = std::sqrt(std::max<Real>(in.initActiveStiffness, 0.0));
  s0.beta = (s0.gamma > 0.0) ? (in.initActiveStress / s0.gamma) : 0.0;
  s0.kc = s0.gamma * s0.gamma;
  s0.tauc = s0.gamma * s0.beta;
  s0.w = in.m0(s0.ec);
  model.initialize(s0);
  if (isRoot) {
    const auto &s = model.getState();
    std::cout << "Initial state:\n"
              << "  y     = " << s.y << '\n'
              << "  v     = " << s.v << '\n'
              << "  pv    = " << s.pv << '\n'
              << "  par   = " << s.par << '\n'
              << "  pd    = " << s.pd << '\n'
              << "  ec    = " << s.ec << '\n'
              << "  gamma = " << s.gamma << '\n'
              << "  beta  = " << s.beta << '\n'
              << "  w     = " << s.w << '\n'
              << "  kc    = " << s.kc << '\n'
              << "  tauc  = " << s.tauc << '\n';
  }
  std::ofstream out;
  if (isRoot) {
    out.open("ccmlc2014_0d_cycle.csv");
    out << "t,y,v,pv,par,pd,ec,gamma,beta,w,kc,tauc,V,Q,pat\n";
  }


  const char *meshFile = "../resources/examples/Heart/malla_solido.mesh";

  Rodin::MPI::Sharder sharder(context);
  if (comm.rank() == RootRank)
  {
    Geometry::Mesh<Context::Local> local;
    local.load(meshFile, IO::FileFormat::MEDIT);

    if (local.getSpaceDimension() != 3)
      throw std::runtime_error("Expected a 3D coronary solid mesh.");

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

  Geometry::Mesh<Context::MPI> mesh = sharder.gather(RootRank);
  mesh.scale(1.e-3);

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


  std::array<Attribute, 6> SolidOutlets{{18, 19, 20, 21, 22, 31}};
  Attribute SolidRing = 17;
  Attribute SolidExt = 102;



  // ---- Finite-element space -----------------------------------------------
  const size_t dim = mesh.getSpaceDimension();
  H1 Vh(std::integral_constant<size_t,1>{}, mesh, dim);   // vector P1 on the distributed mesh

  // ---- Laplacian "map" problem (PETSc KSP / CG) ---------------------------
  //H1 V_lh(mesh);
  //PETSc::Variational::TrialFunction u_l(V_lh);
  //PETSc::Variational::TestFunction  v_l(V_lh);
  //auto zero_lap = RealFunction{Zero()};

  //Problem Laplacian(u_l, v_l);
  //Laplacian = Integral(Grad(u_l), Grad(v_l))
    //        + DirichletBC(u_l, zero_lap).on(GammaRing);
 // {
   // Laplacian.assemble();
   // Solver::KSP lapSolver(Laplacian);
   // lapSolver.setType(KSPCG);
   // Laplacian.solve(lapSolver);
   // }

  //IO::XDMF xdmf_lap("Laplacian");
  //xdmf_lap.setMesh(mesh);
  //xdmf_lap.add("map", u_l.getSolution());
  //xdmf_lap.write(0.0);

  // ---- Solid Material -----------------------------------------------------------
  const Real E = 5e5;
  const Real nu = 0.3;
  const Real lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  const Real mu = E / (2.0 * (1.0 + nu));
  Solid::NeoHookean law(lambda, mu);

  // ---- dynamics parameters ------------------------------------------------
  const Real rho = 1060.;
  const Real dt = 1e-3;
  const int nSteps = 3 * static_cast<int>(0.85 / dt);
  // Newmark average acceleration
  const Real beta = 0.25;
  const Real gamma = 0.5;
  // Effective coefficients
  const Real aMass = rho / (beta * dt * dt);

  // ---- fields ---------------------------------
  PETSc::Variational::GridFunction u(Vh);   // displacement
  PETSc::Variational::GridFunction vel(Vh); // velocity
  PETSc::Variational::GridFunction acc(Vh); // acceleration
  auto zero = VectorFunction{Zero(), Zero(), Zero()};
  u   = zero;
  vel = zero;
  acc = zero;

  PETSc::Variational::GridFunction uPred(Vh);
  PETSc::Variational::GridFunction vPred(Vh);
  PETSc::Variational::GridFunction aNew(Vh);
  PETSc::Variational::GridFunction tmp(Vh);
  uPred = zero;
  vPred = zero;
  aNew  = zero;
  tmp   = zero;

  // ---- output -------------------------------------------------------------
  // Distributed (MPI) XDMF: each rank writes its own HDF5 piece and the root
  // rank writes the master .xdmf that references them. Using the serial
  // constructor on an MPI mesh makes every rank fight over the same file
  // (HDF5 lock errors) and produces broken geometry in ParaView.
  IO::XDMF xdmf(comm, "CoronarySolid", RootRank);
  xdmf.setMesh(mesh);
  xdmf.add("Displacement", u);
  xdmf.add("Velocity", vel);
  xdmf.add("Acceleration", acc);
  xdmf.write(0.0).flush();

  PETSc::Variational::TrialFunction du(Vh);
  PETSc::Variational::TestFunction  w(Vh);

  const auto normal = BoundaryNormal(mesh);

  const Real k = 1e5;
  const Real a = 1e5;
  const Real b = 5e4;
  const Real aVel = b * gamma / (beta * dt);


  Real disp0DValue = 0.0;
  auto disp_0D = RealFunction([&](const Geometry::Point&) { return disp0DValue; });

  Solid::MaterialTangent tangent(law, du, w, u);
  Solid::InternalForce   internal(law, w, u);

  Problem newton(du, w);
  newton = tangent + aMass * Integral(du, w) + internal
         + aMass * Integral(u, w)
         + k * BoundaryIntegral(Dot(du, normal), Dot(w, normal)).over(SolidExt)
         + k * BoundaryIntegral(Dot(u,  normal), Dot(w, normal)).over(SolidExt)
         - k * BoundaryIntegral(disp_0D, Dot(w, normal)).over(SolidExt)
         - aMass * Integral(uPred, w)
         + DirichletBC(du, zero).on(SolidRing);
         + a   * BoundaryIntegral(du, w).over(SolidOutlets[0], SolidOutlets[1], SolidOutlets[2], SolidOutlets[3], SolidOutlets[4], SolidOutlets[5])
         + aVel * BoundaryIntegral(du, w).over(SolidOutlets[0], SolidOutlets[1], SolidOutlets[2], SolidOutlets[3], SolidOutlets[4], SolidOutlets[5])
         + a    * BoundaryIntegral(u, w).over(SolidOutlets[0], SolidOutlets[1], SolidOutlets[2], SolidOutlets[3], SolidOutlets[4], SolidOutlets[5])
         + aVel * BoundaryIntegral(u, w).over(SolidOutlets[0], SolidOutlets[1], SolidOutlets[2], SolidOutlets[3], SolidOutlets[4], SolidOutlets[5])
         - aVel * BoundaryIntegral(uPred, w).over(SolidOutlets[0], SolidOutlets[1], SolidOutlets[2], SolidOutlets[3], SolidOutlets[4], SolidOutlets[5])
         + b    * BoundaryIntegral(vPred, w).over(SolidOutlets[0], SolidOutlets[1], SolidOutlets[2], SolidOutlets[3], SolidOutlets[4], SolidOutlets[5]);

  newton.assemble();

  Solver::KSP ksp(newton);
  Solver::SNES snes(ksp);
  // setTolerances(abstol, rtol, stol, maxIt, maxF)
  snes.setTolerances(1e-10, 1e-8, 1e-10, 50, 10000)
      .setStateUpdate([&](const PETSc::Math::Vector& x)
      {
        // setData refreshes ghosts on MPI grid functions.
        u.setData(x, 0);
      });

  for (size_t step = 1; step <= static_cast<size_t>(nSteps); ++step) {
    const Real t = step * dt;

    // ---- 0D step (redundant on every rank) --------------------------------
    const auto rep = model.step(dt);
    if (isRoot) {
      std::cout << "  0D Newton: "
                << (rep.converged ? "converged" : "not converged")
                << ", iterations = " << rep.iterations
                << ", final residual = " << rep.finalResidual
                << ", final step norm = " << rep.finalStepNorm << '\n';
    }
    if (!rep.converged) {
      if (isRoot)
        std::cerr << "0D solver failed to converge at step " << step
                  << ", t = " << model.getState().t << "\n";
      break;
    }
    const auto &s = model.getState();
    if (isRoot) {
      const Real R = in.R0 + s.y;
      const Real V = (4.0 / 3.0) * std::numbers::pi_v<Real> * R * R * R;
      const Real Q = 4.0 * std::numbers::pi_v<Real> * R * R * s.v;
      const Real pat = in.pAt(s.t);
      out << s.t << "," << s.y << "," << s.v << "," << s.pv << "," << s.par
          << "," << s.pd << "," << s.ec << "," << s.gamma << "," << s.beta
          << "," << s.w << "," << s.kc << "," << s.tauc << "," << V << "," << Q
          << "," << pat << "\n";
    }

    // 0D -> 3D imposed displacement for this step (same value on all ranks).
    disp0DValue = s.y;

    // ---- Newmark predictors (ghost-aware GridFunction algebra) -------------
    // uPred = u + dt*vel + dt*dt*(0.5 - beta)*acc
    uPred = u;
    tmp = vel; tmp *= dt;                         uPred += tmp;
    tmp = acc; tmp *= dt * dt * (0.5 - beta);     uPred += tmp;
    // vPred = vel + dt*(1 - gamma)*acc
    vPred = vel;
    tmp = acc; tmp *= dt * (1.0 - gamma);         vPred += tmp;

    // Initial guess for u^{n+1}: the predictor (state only; SNES x warm-starts
    // from the previous step on its own).
    u = uPred;

    // ---- nonlinear solve with SNES ----------------------------------------
    snes.solve();

    if (!snes.converged()) {
      if (isRoot)
        std::cerr << "SNES failed to converge at step " << step
                  << " after " << snes.getIterationNumber() << " iterations.\n";
      break;
    }
    if (isRoot)
      std::cout << "Step " << step << ", time " << t
                << ", SNES iterations = " << snes.getIterationNumber()
                << std::endl;

    // ---- Newmark correctors ----------------------------------------------
    // aNew = (1/(beta*dt^2)) * (u - uPred)
    aNew = u; aNew -= uPred; aNew *= 1.0 / (beta * dt * dt);
    // vel = vPred + gamma*dt*aNew
    vel = vPred; tmp = aNew; tmp *= gamma * dt; vel += tmp;
    // acc = aNew
    acc = aNew;

    xdmf.write(t).flush();
  }

  xdmf.close();
  if (isRoot)
    out.close();

  PetscFinalize();
  return 0;
}
