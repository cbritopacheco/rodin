#include <cmath>
#include <fstream>
#include <vector>

#include "Rodin/Heart/CCMLC2014.h"

using Real = Rodin::Real;
using Model = Rodin::Heart::CCMLC2014;

static Real periodic_activation(Real t)
{
  const Real T = 0.85; // match Fig. 10 time scale roughly
  const Real tau = t - T * std::floor(t / T);

  // crude first approximation of Fig. 7(b)
  if (tau < 0.03) return 30.0 * (tau / 0.03);
  if (tau < 0.12) return 30.0;
  if (tau < 0.20) return 30.0 * (1.0 - (tau - 0.12) / 0.08);
  if (tau < 0.28) return -20.0 * ((tau - 0.20) / 0.08);
  if (tau < 0.45) return -20.0;
  if (tau < 0.55) return -20.0 * (1.0 - (tau - 0.45) / 0.10);
  return 0.0;
}

static Real n0_piecewise(Real ec)
{
  // start with the same piecewise-linear law used in later code versions
  const Real x1 = -0.4, y1 = 0.0;
  const Real x2 =  0.3, y2 = 0.38;
  const Real x3 =  0.73, y3 = 0.74;
  const Real x4 =  1.0, y4 = 1.0;
  const Real x5 =  1.3, y5 = 1.0;
  const Real x6 =  2.4, y6 = 0.0;

  Real n = 0.0;
  if (ec < x2) n = ((y2-y1)/(x2-x1)) * (ec-x2) + y2;
  else if (ec < x3) n = ((y3-y2)/(x3-x2)) * (ec-x3) + y3;
  else if (ec < x4) n = ((y4-y3)/(x4-x3)) * (ec-x4) + y4;
  else if (ec < x5) n = y4;
  else if (ec < x6) n = ((y6-y5)/(x6-x5)) * (ec-x6) + y6;

  return std::max<Real>(n, 0.0);
}

static Real m0_piecewise(Real ec)
{
  // only a first guess; must be tuned from Fig. 7(a)
  if (ec <= 0.8) return 1.6;
  if (ec <= 1.1) return 1.6 - (ec - 0.8) * (0.6 / 0.3);
  return 1.0;
}

static Real m0_prime_piecewise(Real ec)
{
  if (ec < 0.8) return 0.0;
  if (ec < 1.1) return -0.6 / 0.3;
  return 0.0;
}

int main()
{
  Model::Input in;

  // 1D-calibrated values reused in 0D
  in.rho = 1.0e3;
  in.Es = 3.0e7;
  in.mu = 70.0;
  in.eta = 70.0;
  in.alpha = 1.5;
  in.alphaR = 0.12;

  // 0D-specific values from Table 1
  in.R0 = 2.36e-2;
  in.d0 = 1.42e-2;
  in.k0 = 1.0e5;
  in.sigma0 = 1.24e5;
  in.Rp = 8.0e6;
  in.Cp = 2.5e-9;
  in.Rd = 1.0e8;
  in.Cd = 1.0e-8;
  in.Kat = 9.0e-6;
  in.Kp = 5.0e-10;
  in.Kar = 1.3e-5;

  in.pSv = [](Real) { return 1.0e3; };
  in.pAt = [](Real t) {
    const Real T = 0.85;
    const Real tau = t - T * std::floor(t / T);
    if (tau < 0.10) return 4.0e3;
    if (tau < 0.20) return 8.0e3;
    if (tau < 0.35) return 4.0e3;
    return 2.0e3;
  };

  in.u = periodic_activation;
  in.n0 = n0_piecewise;
  in.m0 = m0_piecewise;
  in.m0Prime = m0_prime_piecewise;

  {
    using PassiveEnergy = std::decay_t<decltype(in.passiveEnergy)>;
    PassiveEnergy::Parameters hp;
    hp.mu1 = 0.0;
    hp.mu2 = 0.0;
    hp.C0 = 1.9e3;
    hp.C1 = 1.1e-1;
    hp.C2 = 1.9e3;
    hp.C3 = 1.1e-1;
    in.passiveEnergy = PassiveEnergy(hp);
  }

  Model model(in);
  model.setMaxIterations(80)
       .setAbsoluteTolerance(1e-10)
       .setRelativeTolerance(1e-10)
       .setStepTolerance(1e-12)
       .setDampingFactor(1.0);

  Model::State s0;
  s0.y = 5.0e-3;
  s0.v = 0.0;
  s0.pv = in.pAt(0.0) + 1;
  s0.par = 1.0e4;
  s0.pd  = 8.0e3;
  const Real e1D0 = 0.5 * (std::pow(1.0 + s0.y / in.R0, 2) - 1.0);
  s0.ec = e1D0;
  s0.kc = 0.0;
  s0.tauc = 0.0;
  s0.w = in.m0(s0.ec);
  s0.t = 0.0;

  model.initialize(s0);

  const Real dt = 1e-4;
  const int nsteps = 15000;

  std::ofstream out("ccmlc2014_0d_cycle.csv");
  out << "t,y,v,pv,par,pd,ec,kc,tauc,w,V,Q,pat\n";

  for (int i = 0; i < nsteps; ++i)
  {
    std::cout << "Step " << i << ": t = " << model.getState().t << "\n";
    const auto rep = model.step(dt);
    if (!rep.converged)
    {
      std::cerr << "Solver failed to converge at step " << i << ", t = " << model.getState().t << "\n";
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
        << s.kc << ","
        << s.tauc << ","
        << s.w << ","
        << V << ","
        << Q << ","
        << pat << "\n";
  }
}
