#include <cmath>
#include <fstream>
#include <iostream>
#include <numbers>
#include <type_traits>
#include <algorithm>

#include "Rodin/Heart/CCMLC2014.h"

using Real = Rodin::Real;
using Model = Rodin::Heart::CCMLC2014T<>;

static Real periodic_activation(Real t)
{
  const Real T = 0.85;
  const Real tau = t - T * std::floor(t / T);

  if (tau < 0.13)  return 0.0;
  if (tau < 0.141) return 35.0 * ((tau - 0.13) / 0.011);
  if (tau < 0.281) return 35.0;
  if (tau < 0.361) return 35.0 - 55.0 * ((tau - 0.281) / 0.08);
  if (tau < 0.45)  return -20.0;
  return 0.0;
}

static Real load_dependent_relaxation_m0(Real ec)
{
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

static Real load_dependent_relaxation_dm0(Real ec)
{
  const Real low_ec = 0.0;
  const Real high_ec = 2.0;
  const Real low_value = 1.6;
  const Real high_value = 1.0;

  if (ec <= low_ec || ec >= high_ec)
    return 0.0;
  return (high_value - low_value) / (high_ec - low_ec);
}

static Real atrial_pressure(Real t)
{
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

  if (tau < t1)
  {
    alpha = -(tau - t1) / t1;
    value = alpha * min_value + (1.0 - alpha) * max_value;
  }
  else if (tau < t2)
  {
    value = max_value;
  }
  else if (tau < t3)
  {
    alpha = -(tau - t3) / (t3 - t2);
    value = alpha * max_value + (1.0 - alpha) * min_value;
  }
  else if (tau < t4)
  {
    alpha = -(tau - t4) / (t4 - t3);
    value = alpha * min_value + (1.0 - alpha) * second_threshold;
  }
  else if (tau < t5)
  {
    value = second_threshold;
  }
  else if (tau < t6)
  {
    alpha = -(tau - t6) / (t6 - t5);
    value = alpha * second_threshold + (1.0 - alpha) * min_value;
  }
  else
  {
    value = min_value;
  }

  return value;
}

int main()
{
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

  in.mu_0 = 0.04868;
  in.mu_Inf = 0.003605;
  in.lambda = 3.39;
  in.n = 0.198;
  in.yasuda = 1.235;
  in.proximalRadius = 0.015;
  in.proximalLength = 0.4;
  in.distalRadius = 0.0007;
  in.distalLength = 0.004;
  in.windkesselRheology =
    Rodin::Heart::CCMLC2014::Model::WindkesselRheology::CarreauYasuda;

  // Valve parameters
  in.Kat = 9.0e-6;
  in.Kp  = 5.0e-10;
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

  s0.y   = 0.0;
  s0.v   = 0.0;
  s0.pv  = in.pAt(0.0) - 100.0;
  s0.par = 11000.0;
  s0.pd  = 10000.0;

  // Initial local active state
  s0.ec = in.initFibDef;
  s0.gamma = std::sqrt(std::max<Real>(in.initActiveStiffness, 0.0));
  s0.beta = (s0.gamma > 0.0) ? (in.initActiveStress / s0.gamma) : 0.0;
  s0.kc = s0.gamma * s0.gamma;
  s0.tauc = s0.gamma * s0.beta;
  s0.w = in.m0(s0.ec);

  model.initialize(s0);

  {
    const auto& s = model.getState();
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

  const Real dt = 1e-3;
  const int nsteps = 3 * static_cast<int>(0.85 / dt);

  std::ofstream out("ccmlc2014_0d_cycle.csv");
  out << "t,y,v,pv,par,pd,ec,gamma,beta,w,kc,tauc,V,Q,pat\n";

  for (int i = 0; i < nsteps; ++i)
  {
    std::cout << "Step " << i << ": t = " << model.getState().t << "\n";

    const auto rep = model.step(dt);
    std::cout << "  Newton step: "
              << (rep.converged ? "converged" : "not converged")
              << ", iterations = " << rep.iterations
              << ", final residual = " << rep.finalResidual
              << ", final step norm = " << rep.finalStepNorm
              << '\n';

    if (!rep.converged)
    {
      std::cerr << "Solver failed to converge at step "
                << i << ", t = " << model.getState().t << "\n";
      break;
    }

    const auto& s = model.getState();

    const Real R = in.R0 + s.y;
    const Real V = (4.0 / 3.0) * std::numbers::pi_v<Real> * R * R * R;
    const Real Q = 4.0 * std::numbers::pi_v<Real> * R * R * s.v;
    const Real pat = in.pAt(s.t);

    out << s.t << ","
        << s.y << ","
        << s.v << ","
        << s.pv << ","
        << s.par << ","
        << s.pd << ","
        << s.ec << ","
        << s.gamma << ","
        << s.beta << ","
        << s.w << ","
        << s.kc << ","
        << s.tauc << ","
        << V << ","
        << Q << ","
        << pat << "\n";
  }

  return 0;
}
