/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file AtriumFlow.cpp
 * @brief Pressure-driven Newtonian flow in the left atrium on a fixed domain,
 *        assembled and solved on the device.
 *
 * First step of the port of examples/Heart/Atrium/Viscoelastic
 * (model/Atrium) to the Kokkos module: the openings, the lumped ventricle and
 * the time loop of that driver, with a linear Oseen step per time step. See
 * README.md of this folder for what is and is not ported.
 *
 * Each step solves, for P1/P1 velocity and pressure and w = u^n,
 *
 *   rho/dt (u, v) + rho/2 [((grad u) w, v) - (u, (grad v) w)] + mu (grad u, grad v)
 *   - (p, div v) + (div u, q)
 *   + rho (tau (grad u) w, (grad v) w) + (tau/rho) (grad p, grad q)
 *   + <p_pv n, v>_veins + <p_out n, v>_mitral
 *   + Z_pv <u, v>_veins + Z_m <u.n, v.n>_mitral  =  rho/dt (w, v),
 *
 * with u = 0 on the wall. The convection is skew-symmetric, so it neither
 * produces nor removes kinetic energy whatever the direction of the flow
 * through an opening; this replaces the backflow term of Atrium/Fluid. tau is
 * the cellwise parameter [(2/dt)^2 + (2|w|/h)^2 + (4 nu/h^2)^2]^(-1/2),
 * evaluated on the device from w and the cell volume. The streamline term is
 * not the residual-based stabilisation of the CPU drivers: it has no
 * orthogonal projection and no dynamic subscale.
 *
 * Openings, as in Atrium/Fluid (labels of the Atrium meshes): veins 1, 2, 3,
 * 5 with p_pv = mean + amplitude sin(2 pi t/T) and impedance Z_pv; mitral
 * plane 6 with the lumped ventricle Rodin::Heart::CCMLC2014, p_out = p_v and
 * Z_m = |Gamma|/K, K the conductance of its valve law, open while the mean
 * pressure of the plane exceeds p_v; wall 4. The ventricle runs on the host:
 * it is a system of five scalar equations.
 *
 * The mesh is a tetrahedral MEDIT mesh of the fluid in millimetres
 * (resources/examples/Heart/Atrium/Coarse.mesh), scaled to metres.
 *
 * Usage:
 *   KokkosAtriumFlow mesh.mesh [key=value ...]
 * Keys: cycles (2), steps (850, per cycle), period (0.85), rho (1060),
 * mu (0.0049), p_in (1200), p_in_amplitude (200), ramp (0.25, cycles),
 * z_in (1000), lv (1), p_out (0), output (AtriumFlow), output_every (17),
 * tolerance (1e-8), restart (200), scale (1e-3).
 * Output: <output>.xdmf (velocity, pressure) and <output>.csv, one row per step.
 * With cuDSS (RODIN_USE_CUDSS) or RocSolverRF the step is solved by a sparse
 * direct solver on the device; otherwise by GMRES.
 */
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Heart/CCMLC2014.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Kokkos.h>
#include <Rodin/Kokkos/Solver.h>
#include <Rodin/Solver.h>
#include <Rodin/Variational.h>

#if defined(KOKKOS_ENABLE_CUDA)
using Backend = Rodin::Context::CUDA;
#elif defined(KOKKOS_ENABLE_HIP)
using Backend = Rodin::Context::HIP;
#elif defined(KOKKOS_ENABLE_SYCL)
using Backend = Rodin::Context::SYCL;
#else
using Backend = ::Kokkos::DefaultHostExecutionSpace;
#endif

// A global using directive on Rodin would make Kokkos ambiguous in the Kokkos headers.
using DeviceContext = Rodin::Context::Kokkos<Rodin::Context::Local, Backend>;
using Ventricle = Rodin::Heart::CCMLC2014T<>;

namespace
{
  using Real = Rodin::Real;

  /// @brief Run parameters, set by key=value arguments.
  struct Parameters
  {
    Real scale = 1e-3;
    Real period = 0.85;
    size_t cycles = 2;
    size_t steps = 850;
    Real rho = 1060;
    Real mu = 0.0049;
    Real inletPressureMean = 1200;
    Real inletPressureAmplitude = 200;
    Real rampCycles = 0.25;
    Real inletImpedance = 1000;
    bool ventricle = true;
    Real outletPressure = 0;
    std::string output = "AtriumFlow";
    size_t outputEvery = 17;
    Real tolerance = 1e-8;
    size_t restart = 200;
  };

  /// @brief Smooth step on [0, 1].
  Real smooth(Real s)
  {
    s = s < 0 ? 0 : (s > 1 ? 1 : s);
    return s * s * (3 - 2 * s);
  }

  /// @brief Ventricular activation of Atrium/Fluid, periodic.
  Real activation(Real period, Real t)
  {
    constexpr Real rampStart = 0.15, rampEnd = 0.21, plateauEnd = 0.36, relaxEnd = 0.45,
                   negativeEnd = 0.6, positive = 35, negative = -20;
    const Real tau = t - period * std::floor(t / period);
    if (tau < rampStart)
      return 0;
    if (tau < rampEnd)
      return positive * smooth((tau - rampStart) / (rampEnd - rampStart));
    if (tau < plateauEnd)
      return positive;
    if (tau < relaxEnd)
    {
      return positive +
        (negative - positive) * smooth((tau - plateauEnd) / (relaxEnd - plateauEnd));
    }
    if (tau < negativeEnd)
      return negative;
    return negative * (1 - smooth((tau - negativeEnd) / (period - negativeEnd)));
  }

  /**
   * @brief Lumped ventricle of Atrium/Fluid, its parameters unchanged.
   * @param period Cardiac period.
   * @param atrialPressure Mean pressure of the mitral plane, read at every step.
   */
  Ventricle::Input ventricleInput(Real period, const Real* atrialPressure)
  {
    Ventricle::Input input;
    input.rho = 1.0e3;
    input.R0 = 2.45e-2;
    input.d0 = 1.4e-2;
    input.Es = 3.0e6;
    input.mu = 70.0;
    input.eta = 70.0;
    input.alpha = 1.5;
    input.alphaR = 0.12;
    input.k0 = 1.0e5;
    input.sigma0 = 1.2e5;
    input.Rp = 5.0e7;
    input.Cp = 5.0e-9;
    input.Rd = 1.0e8;
    input.Cd = 5.0e-10;
    input.proximalRadius = 0.0125;
    input.proximalLength = 0.4;
    input.distalRadius = 0.00175;
    input.distalLength = 0.2;
    input.mu_0 = 0.301;
    input.mu_Inf = 0.0055;
    input.lambda = 16.152;
    input.n = 0.21;
    input.yasuda = 0.77;
    input.Kat = 1.0e-6;
    input.Kp = 5.0e-10;
    input.Kar = 2.0e-7;
    input.cavityCapacity = 5.0e-12;
    input.localTolerance = 1.0e-12;
    input.localMaxIterations = 50;
    input.localDamping = 1.0;
    input.absRegularization = 1.0e-14;
    input.pSv = [](Real) { return Real(1.0e3); };
    input.pAt = [atrialPressure](Real) { return *atrialPressure; };
    input.u = [period](Real t) { return activation(period, t); };
    input.m0 = [](Real ec) {
      constexpr Real low = 1.6, high = 1.0, lowEc = 0.0, highEc = 2.0;
      if (ec >= highEc)
        return high;
      if (ec <= lowEc)
        return low;
      const Real s = (ec - lowEc) / (highEc - lowEc);
      return (1 - s) * low + s * high;
    };
    input.dm0 = [](Real ec) {
      constexpr Real low = 1.6, high = 1.0, lowEc = 0.0, highEc = 2.0;
      if (ec <= lowEc || ec >= highEc)
        return Real(0);
      return (high - low) / (highEc - lowEc);
    };
    using PassiveEnergy = std::decay_t<decltype(input.passiveEnergy)>;
    typename PassiveEnergy::Parameters hp;
    hp.C0 = 1.9e3;
    hp.C1 = 1.1e-1;
    hp.C2 = 1.9e3;
    hp.C3 = 1.1e-1;
    input.passiveEnergy = PassiveEnergy(hp);
    return input;
  }

  /// @brief Reads key=value arguments into the parameters.
  Parameters parse(int argc, char** argv)
  {
    Parameters P;
    std::map<std::string, std::string> kv;
    for (int i = 2; i < argc; ++i)
    {
      const std::string a = argv[i];
      const auto eq = a.find('=');
      if (a.rfind("--", 0) == 0 || eq == std::string::npos)
        continue; // Kokkos runtime arguments.
      kv[a.substr(0, eq)] = a.substr(eq + 1);
    }
    const auto real = [&](const char* key, Real& out) {
      if (kv.count(key))
        out = std::stod(kv[key]);
    };
    const auto count = [&](const char* key, size_t& out) {
      if (kv.count(key))
        out = std::stoul(kv[key]);
    };
    real("scale", P.scale);
    real("period", P.period);
    count("cycles", P.cycles);
    count("steps", P.steps);
    real("rho", P.rho);
    real("mu", P.mu);
    real("p_in", P.inletPressureMean);
    real("p_in_amplitude", P.inletPressureAmplitude);
    real("ramp", P.rampCycles);
    real("z_in", P.inletImpedance);
    real("p_out", P.outletPressure);
    real("tolerance", P.tolerance);
    count("restart", P.restart);
    count("output_every", P.outputEvery);
    if (kv.count("lv"))
      P.ventricle = kv["lv"] != "0";
    if (kv.count("output"))
      P.output = kv["output"];
    return P;
  }
}

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;

  if (argc < 2)
  {
    std::cerr << "Usage: KokkosAtriumFlow mesh.mesh [key=value ...]" << std::endl;
    return EXIT_FAILURE;
  }
  ::Kokkos::ScopeGuard runtime(argc, argv);
  const Parameters P = parse(argc, argv);

  // Labels of the Atrium meshes.
  constexpr Attribute wall = 4;
  constexpr Attribute mitral = 6;
  constexpr std::array<Attribute, 4> veins{1, 2, 3, 5};

  Mesh host;
  host.load(argv[1], IO::FileFormat::MEDIT);
  if (host.getDimension() != 3 || host.getCellCount() == 0 ||
    host.getPolytopeCount(Polytope::Type::Tetrahedron) != host.getCellCount())
  {
    std::cerr << argv[1] << " is not a tetrahedral mesh of the fluid (a Git LFS pointer, "
              << "or a combined fluid-solid mesh?)." << std::endl;
    return EXIT_FAILURE;
  }
  host.scale(P.scale);
  // Face to cell incidence for the boundary; the remaining incidences for the
  // host space of the XDMF output.
  for (size_t d = 1; d <= 3; ++d)
  {
    for (size_t e = 0; e < d; ++e)
    {
      host.getConnectivity().compute(d, e);
      host.getConnectivity().compute(e, d);
    }
  }
  const Mesh<DeviceContext> mesh(DeviceContext{}, host);

  H1 vh(std::integral_constant<size_t, 1>{}, mesh, 3);
  H1 qh(std::integral_constant<size_t, 1>{}, mesh);
  TrialFunction u(vh);
  TrialFunction p(qh);
  TestFunction v(vh);
  TestFunction q(qh);
  GridFunction velocity(vh);
  GridFunction pressure(qh);
  velocity = Zero(3);
  pressure = 0.0;
  std::cout << "[mesh] " << argv[1] << "  cells=" << host.getCellCount()
            << " vertices=" << host.getVertexCount()
            << " unknowns=" << (vh.getSize() + qh.getSize()) << std::endl;

  // Areas of the openings, and the functionals that give the flow rates and
  // the mean pressure of the mitral plane: int_G phi n and int_G psi, assembled
  // once on the device and kept on the host.
  const BoundaryNormal normal(mesh);
  LinearForm veinFlux(v);
  LinearForm mitralFlux(v);
  LinearForm veinMean(q);
  LinearForm mitralMean(q);
  veinFlux = BoundaryIntegral(normal, v).over(veins[0], veins[1], veins[2], veins[3]);
  mitralFlux = BoundaryIntegral(normal, v).over(mitral);
  veinMean =
    BoundaryIntegral(RealFunction(1), q).over(veins[0], veins[1], veins[2], veins[3]);
  mitralMean = BoundaryIntegral(RealFunction(1), q).over(mitral);
  veinFlux.assemble();
  mitralFlux.assemble();
  veinMean.assemble();
  mitralMean.assemble();
  const auto veinFluxHost = veinFlux.getVector().toHost();
  const auto mitralFluxHost = mitralFlux.getVector().toHost();
  const auto veinMeanHost = veinMean.getVector().toHost();
  const auto mitralMeanHost = mitralMean.getVector().toHost();
  const Real veinArea = veinMeanHost.sum();
  const Real mitralArea = mitralMeanHost.sum();
  const Real volume = host.getVolume();
  const Real h = std::cbrt(volume / static_cast<Real>(host.getCellCount()));
  std::cout << "[boundary] venous area = " << veinArea
            << " m^2  mitral area = " << mitralArea << " m^2  volume = " << 1e6 * volume
            << " mL  mean cell size = " << 1e3 * h << " mm" << std::endl;
  if (!(veinArea > 0) || !(mitralArea > 0))
  {
    std::cerr << "No venous or mitral faces: labels 1, 2, 3, 5 and 6 are expected."
              << std::endl;
    return EXIT_FAILURE;
  }

  // Lumped ventricle, on the host. It reads the mean pressure of the mitral
  // plane at the last step as its atrial pressure.
  Real atrialPressure = P.rampCycles > 0 ? Real(0) : P.inletPressureMean;
  Ventricle::Input input = ventricleInput(P.period, &atrialPressure);
  Ventricle ventricle(input);
  {
    Ventricle::State s0;
    s0.t = 0.0;
    s0.y = 0.0;
    s0.v = 0.0;
    s0.pv = atrialPressure - 100.0;
    s0.par = 11000.0;
    s0.pd = 10000.0;
    s0.w = input.m0(s0.ec);
    ventricle.setMaxIterations(200)
      .setAbsoluteTolerance(1.0e-8)
      .setRelativeTolerance(1.0e-8)
      .setStepTolerance(1.0e-10)
      .setDampingFactor(1.0);
    ventricle.initialize(s0);
  }

  IO::XDMF xdmf(P.output);
  xdmf.grid().setMesh(mesh).add("velocity", velocity).add("pressure", pressure);
  xdmf.write(0).flush();
  std::ofstream csv(P.output + ".csv");
  csv
    << "t,cycle,pVenous,pOut,valveOpen,qIn,qOut,pInletMean,pOutletMean,maxSpeedComponent,"
       "assembly,solve\n";

  // The problem and the solver outlive the time loop: the sparsity pattern is
  // the same at every step, so a direct device solver keeps its analysis.
#if defined(RODIN_USE_CUDSS) || defined(RODIN_USE_ROCSOLVER)
#define DEVICE_SOLVER
  using System = Math::LinearSystem<Rodin::Kokkos::Math::SparseMatrix<Real, Backend>,
    Rodin::Kokkos::Math::Vector<Real, Backend>>;
  Problem<System, decltype(u), decltype(p), decltype(v), decltype(q)> flow(u, p, v, q);
#ifdef RODIN_USE_CUDSS
  Solver::CuDSS solver(flow);
#else
  Solver::RocSolverRF solver(flow);
#endif
#else
  Problem flow(u, p, v, q);
  Solver::GMRES solver(flow);
  solver.setTolerance(P.tolerance).setRestart(P.restart);
#endif

  const Real dt = P.period / static_cast<Real>(P.steps);
  const Real rho = P.rho;
  const Real mu = P.mu;
  const Real nu = mu / rho;
  const Real pi = Math::Constants::pi();
  const size_t total = P.cycles * P.steps;
  std::chrono::duration<double> assembly{0}, solve{0};
  for (size_t step = 1; step <= total; ++step)
  {
    const Real time = step * dt;
    const Real rampTime = P.rampCycles * P.period;
    const Real ramp = rampTime > 0 && time < rampTime
      ? 0.5 * (1 - std::cos(pi * time / rampTime))
      : Real(1);
    const Real pIn = ramp *
      (P.inletPressureMean +
        P.inletPressureAmplitude * std::sin(2 * pi * time / P.period));

    // Mitral plane: ventricular pressure and valve law, or a prescribed pressure.
    Real pOut = ramp * P.outletPressure;
    Real mitralImpedance = 0;
    bool valveOpen = true;
    if (P.ventricle)
    {
      if (!ventricle.step(dt).converged)
      {
        std::cerr << "The lumped ventricle did not converge at t = " << time << std::endl;
        return EXIT_FAILURE;
      }
      pOut = ventricle.getState().pv;
      valveOpen = atrialPressure > pOut;
      mitralImpedance = mitralArea / (valveOpen ? input.Kat : input.Kp);
    }

    // tau by cell, from w = u^n and the cell volume |det J|/6 of a tetrahedron.
    const ::Rodin::Kokkos::Variational::View<std::decay_t<decltype(velocity)>> w(
      velocity);
    RealFunction tau(
      KOKKOS_LAMBDA(const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>& point)->Real {
        const auto wp = w.getValue(point);
        Real speed2 = 0;
        for (Index i = 0; i < 3; ++i)
          speed2 += wp(i) * wp(i);
        const Real hk = ::Kokkos::cbrt(point.getDistortion() / 6);
        const Real a = 2 / dt;
        const Real b2 = 4 * speed2 / (hk * hk);
        const Real c = 4 * nu / (hk * hk);
        return 1 / ::Kokkos::sqrt(a * a + b2 + c * c);
      });
    RealFunction pressureTau(
      KOKKOS_LAMBDA(const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>& point)->Real {
        const auto wp = w.getValue(point);
        Real speed2 = 0;
        for (Index i = 0; i < 3; ++i)
          speed2 += wp(i) * wp(i);
        const Real hk = ::Kokkos::cbrt(point.getDistortion() / 6);
        const Real a = 2 / dt;
        const Real b2 = 4 * speed2 / (hk * hk);
        const Real c = 4 * nu / (hk * hk);
        return 1 / (rho * ::Kokkos::sqrt(a * a + b2 + c * c));
      });

    flow = (rho / dt) * Integral(u, v) +
      (0.5 * rho) * Integral(Dot(Mult(Jacobian(u), velocity), v)) -
      (0.5 * rho) * Integral(Dot(u, Mult(Jacobian(v), velocity))) +
      mu * Integral(Jacobian(u), Jacobian(v)) - Integral(p, Div(v)) +
      Integral(Div(u), q) +
      rho *
        Integral(
          Dot(Mult(tau, Mult(Jacobian(u), velocity)), Mult(Jacobian(v), velocity))) +
      Integral(Dot(Mult(pressureTau, Grad(p)), Grad(q))) +
      pIn * BoundaryIntegral(normal, v).over(veins[0], veins[1], veins[2], veins[3]) +
      pOut * BoundaryIntegral(normal, v).over(mitral) +
      P.inletImpedance *
        BoundaryIntegral(u, v).over(veins[0], veins[1], veins[2], veins[3]) +
      mitralImpedance * BoundaryIntegral(Dot(normal, u), Dot(normal, v)).over(mitral) -
      (rho / dt) * Integral(velocity, v) + DirichletBC(u, Zero(3)).on(wall);

    auto t0 = std::chrono::steady_clock::now();
    flow.assemble();
    auto t1 = std::chrono::steady_clock::now();
#ifdef DEVICE_SOLVER
    solver.factorize(flow.getLinearSystem());
    if (!solver.success())
    {
      std::cerr << "Factorisation failed at step " << step << std::endl;
      return EXIT_FAILURE;
    }
#endif
    solver.solve();
    auto t2 = std::chrono::steady_clock::now();
    if (!solver.success())
    {
      std::cerr << "The linear solver failed at step " << step << std::endl;
      return EXIT_FAILURE;
    }
    assembly += t1 - t0;
    solve += t2 - t1;
    velocity.setData(u.getSolution().getData());
    pressure.setData(p.getSolution().getData());

    // Diagnostics on the host copies of the solution. n is outward: the
    // venous flux is negated to read as the inflow.
    const auto uh = velocity.getData().toHost();
    const auto ph = pressure.getData().toHost();
    const Real qIn = -veinFluxHost.dot(uh);
    const Real qOut = mitralFluxHost.dot(uh);
    const Real pVein = veinMeanHost.dot(ph) / veinArea;
    atrialPressure = mitralMeanHost.dot(ph) / mitralArea;
    const Real speed = uh.cwiseAbs().maxCoeff();
    if (!std::isfinite(speed) || speed > 50)
    {
      std::cerr << "Diverged at step " << step << ": max |u_i| = " << speed << " m/s"
                << std::endl;
      return EXIT_FAILURE;
    }
    const size_t cycle = (step - 1) / P.steps + 1;
    csv << time << ',' << cycle << ',' << pIn << ',' << pOut << ',' << (valveOpen ? 1 : 0)
        << ',' << qIn << ',' << qOut << ',' << pVein << ',' << atrialPressure << ','
        << speed << ',' << std::chrono::duration<double>(t1 - t0).count() << ','
        << std::chrono::duration<double>(t2 - t1).count() << '\n';
    if (step % P.outputEvery == 0 || step == total)
    {
      xdmf.write(time).flush();
      csv.flush();
      std::cout << "Step " << step << " / " << total << "  t=" << time
                << " s  p_pv=" << pIn << "  p_out=" << pOut << " Pa  valve "
                << (valveOpen ? "open" : "closed") << "  qIn=" << qIn << "  qOut=" << qOut
                << " m^3/s  max|u_i|=" << speed << " m/s  assembly "
                << std::chrono::duration<double>(t1 - t0).count() << " s  solve "
                << std::chrono::duration<double>(t2 - t1).count() << " s" << std::endl;
    }
  }
  std::cout << "Assembly " << assembly.count() << " s, solve " << solve.count()
            << " s in " << total << " steps" << std::endl;
  return 0;
}
