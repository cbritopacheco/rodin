/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * The Coronary Artery is modeled as a quasi compressible NeoHookean solid with
 * inertia and simple mass-proportional viscous damping:
 *
 *   rho * u_tt - div(P(F(u))) = 0
 *
 * Time integration is performed with the Newmark-beta method using the
 * average-acceleration choice beta = 1/4, gamma = 1/2.
 *
 * Output is written to XDMF for visualization in ParaView.
 *
 * This version also instruments a discrete derivative check for the assembled
 * nonlinear residual R(u) and tangent J(u):
 *
 *   FD(eps) = (R(u + eps * eta) - R(u)) / eps
 *
 * and compares FD(eps) against J(u) * eta.
 *
 * If the tangent is consistent, ||FD(eps) - J(u) eta|| should scale like O(eps)
 * for sufficiently small eps, until roundoff dominates.
 */
#include "Rodin/IO/ForwardDecls.h"
#include "Rodin/Variational/ForwardDecls.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>

#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Solid.h>
#include <Rodin/Solver/CG.h>
#include <Rodin/Solver/NewtonSolver.h>
#include <Rodin/Solver/SparseLU.h>
#include <Rodin/Variational.h>

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

int main(int, char **) {

  // Initialization 0D model

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

  {
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

  std::ofstream out("ccmlc2014_0d_cycle.csv");
  out << "t,y,v,pv,par,pd,ec,gamma,beta,w,kc,tauc,V,Q,pat\n";

  // Initialization 3D Model

  const char *meshFile = "../resources/examples/Heart/CoronaryArterySolid.mesh";

  // Define boundary attributes
  Attribute GammaRing = 100, Gamma1 = 101, Gamma2 = 105;

  // Load mesh
  Mesh mesh;
  mesh.load(meshFile, IO::FileFormat::MEDIT);
  mesh.scale(1.e-3);
  mesh.getConnectivity().compute(2, 3);

  // ---- Finite-element space -----------------------------------------------
  const size_t dim = mesh.getSpaceDimension();
  P1 Vh(mesh, dim);

  // Laplacian problem
  P1 V_lh(mesh);
  // Fe functions
  TrialFunction u_l(V_lh);
  TestFunction v_l(V_lh);
  auto zero_lap = RealFunction{Zero()};
  //Fields
  Problem Laplacian(u_l, v_l);

  Laplacian = Integral(Grad(u_l), Grad(v_l))
            + DirichletBC(u_l, zero_lap).on(GammaRing);

  CG(Laplacian).solve();

  // Save solution
  IO::XDMF xdmf_lap("Laplacian");
  xdmf_lap.grid().setMesh(mesh).add("map", u_l.getSolution());
  xdmf_lap.write();



  // ---- material -----------------------------------------------------------
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

  // ---- fields -------------------------------------------------------------
  GridFunction u(Vh);   // displacement
  GridFunction vel(Vh); // velocity
  GridFunction acc(Vh); // acceleration

  u.setName("Displacement");
  vel.setName("Velocity");
  acc.setName("Acceleration");

  auto zero = VectorFunction{Zero(), Zero(), Zero()};
  u = zero;
  vel = zero;
  acc = zero;

  GridFunction uPred(Vh);
  GridFunction vPred(Vh);
  GridFunction aNew(Vh);

  uPred.setName("PredictedDisplacement");
  vPred.setName("PredictedVelocity");
  aNew.setName("AccelerationNew");

  // ---- output -------------------------------------------------------------
  IO::XDMF xdmf("CoronarySolid");
  auto grid = xdmf.grid();
  grid.setMesh(mesh);
  grid.add(u);
  grid.add(vel);
  grid.add(acc);

  xdmf.write(0.0);

  TrialFunction du(Vh);
  TestFunction w(Vh);

  for (size_t step = 1; step <= nSteps; ++step) {

    const Real t = step * dt;

    const auto rep = model.step(dt);
    std::cout << "  Newton step: "
              << (rep.converged ? "converged" : "not converged")
              << ", iterations = " << rep.iterations
              << ", final residual = " << rep.finalResidual
              << ", final step norm = " << rep.finalStepNorm << '\n';

    if (!rep.converged) {
      std::cerr << "Solver failed to converge at step " << step
                << ", t = " << model.getState().t << "\n";
      break;
    }

    const auto &s = model.getState();
    const Real R = in.R0 + s.y;
    const Real V = (4.0 / 3.0) * std::numbers::pi_v<Real> * R * R * R;
    const Real Q = 4.0 * std::numbers::pi_v<Real> * R * R * s.v;
    const Real pat = in.pAt(s.t);

    out << s.t << "," << s.y << "," << s.v << "," << s.pv << "," << s.par << ","
        << s.pd << "," << s.ec << "," << s.gamma << "," << s.beta << "," << s.w
        << "," << s.kc << "," << s.tauc << "," << V << "," << Q << "," << pat
        << "\n";


    // Passing the displacement to be imposed in the 3D model
    const Real disp = s.y;
    auto disp_0D = RealFunction{disp};
    const auto normal = BoundaryNormal(mesh);
    const Real k = 1.e4;

    // ---- Newmark predictors -----------------------------------------------
    uPred = u + dt * vel + (dt * dt * (0.5 - beta)) * acc;
    vPred = vel + (dt * (1.0 - gamma)) * acc;

    // Use predictor as initial guess for u^{n+1}
    u = uPred;

    // ---- nonlinear solid operators ----------------------------------------
    Solid::MaterialTangent tangent(law, du, w, u);

    Solid::InternalForce internal(law, w, u);

    // Effective nonlinear problem:
    //
    //   R(u; w)
    //   = F_int(u; w)
    //   + aMass * (u, w)
    //   - F_ext(w)
    //   - aMass * (uPred, w)
    //
    // Newton increment equation:
    //
    //   J(u^k)[du] + R(u^k) = 0
    //
    Problem newton(du, w);
    newton = tangent + aMass * Integral(du, w) + internal +
             aMass * Integral(u, w) -
             + k * BoundaryIntegral(Dot(du, normal), Dot(w,normal)).over(Gamma1,Gamma2)
             + k * BoundaryIntegral(Dot(u, normal), Dot(w,normal)).over(Gamma1,Gamma2)
             - k * BoundaryIntegral(disp_0D, Dot(w,normal)).over(Gamma1,Gamma2)
             - aMass * Integral(uPred, w) + DirichletBC(du, zero).on(GammaRing);

    SparseLU linearSolver(newton);
    NewtonSolver solver(linearSolver);
    solver.setMaxIterations(50)
        .setDampingFactor(1.0)
        .setAbsoluteTolerance(1e-10)
        .setRelativeTolerance(1e-8);
    solver.solve(u);

    std::cout << "Step " << step << ", time " << t << std::endl;

    // ---- Newmark correctors -----------------------------------------------
    aNew = (1.0 / (beta * dt * dt)) * (u - uPred);
    vel = vPred + (gamma * dt) * aNew;
    acc = aNew;

    xdmf.write(t).flush();
  }

  xdmf.close();
  return 0;
}
