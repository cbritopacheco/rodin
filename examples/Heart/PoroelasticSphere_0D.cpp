/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file PoroelasticSphere_0D.cpp
 * @brief Standalone driver for the 0D poroelastic sphere heart model.
 *
 * Integrates the thick-walled, solid-incompressible poroelastic sphere of
 * @ref Rodin::Heart::PoroelasticSphereT, which shares the passive, active,
 * valve and Windkessel laws of CCMLC2014 and adds a porosity balance fed by a
 * two-resistor coronary source, driven by a periodic activation. The parameters
 * are the healthy-adult calibration of CoronaryPoroelastic/HealthyLV.h
 * (Murray-sized Windkessel branches, porosity-dependent perfusion
 * conductances). It writes the resulting pressure, volume, flow and porosity
 * history. No mesh and no finite element space are involved.
 */
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>

#include "CoronaryPoroelastic/HealthyLV.h"

using Real = Rodin::Real;
using namespace Rodin::Examples::Heart::CoronaryPoroelastic;

int main()
{
  HealthyCalibration cal;
  const Model::Input in = makeHealthyInput(HealthyTargets(), &cal);

  std::cout << "Healthy calibration:\n"
            << "  V_w0      = " << cal.wallVolume * 1e6 << " mL\n"
            << "  Q_cor     = " << cal.coronaryFlow * 6e7 << " mL/min (target)\n"
            << "  R_p, R_d  = " << cal.Rp << ", " << cal.Rd << " Pa s/m^3\n"
            << "  Murray r0 = " << cal.proximal.radius * 1e3 << " mm, L_p = "
            << cal.proximal.length << " m (" << cal.proximal.generations
            << " generations), L_d = " << cal.distal.length << " m\n"
            << "  gamma_ar  = " << in.gammaAr << ", gamma_ven = " << in.gammaVen
            << " m^3/(Pa s)\n";

  Model model(in);
  model.setMaxIterations(200)
    .setAbsoluteTolerance(1e-8)
    .setRelativeTolerance(1e-8)
    .setStepTolerance(1e-10)
    .setDampingFactor(1.0);
  model.initialize(makeHealthyInitialState(in));

  {
    const auto& s = model.getState();
    std::cout << "Initial state:\n"
              << "  y     = " << s.y << '\n'
              << "  phi   = " << s.phi << '\n'
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
  const int nsteps = 10 * static_cast<int>(0.85 / dt);
  Real Vprev = (4.0 / 3.0) * std::numbers::pi_v<Real> * in.R0 * in.R0 * in.R0;

  std::ofstream out("poroelastic_sphere_0d_cycle.csv");
  // Full precision: the porosity varies by O(1e-2) about phi0, so the default
  // six significant digits quantize it too coarsely for its discrete rate --
  // and hence for the fluid mass balance -- to be recovered from the file.
  out << std::setprecision(12);
  out << "t,y,phi,pv,par,pd,ec,gamma,beta,w,kc,tauc,V,Q,pat,lambdaBar,pf,Qcor,Qven\n";

  for (int i = 0; i < nsteps; ++i)
  {
    std::cout << "Step " << i << ": t = " << model.getState().t << "\n";

    const auto rep = model.step(dt);
    std::cout << "  Newton step: " << (rep.converged ? "converged" : "not converged")
              << ", iterations = " << rep.iterations
              << ", final residual = " << rep.finalResidual
              << ", final step norm = " << rep.finalStepNorm << '\n';

    if (!rep.converged)
    {
      std::cerr << "Solver failed to converge at step " << i
                << ", t = " << model.getState().t << "\n";
      break;
    }

    const auto& s = model.getState();

    const Real R = in.R0 + s.y;
    const Real V = (4.0 / 3.0) * std::numbers::pi_v<Real> * R * R * R;
    const Real Q = (V - Vprev) / dt;
    Vprev = V;
    const Real pat = in.pAt(s.t);

    out << s.t << "," << s.y << "," << s.phi << "," << s.pv << "," << s.par << "," << s.pd
        << "," << s.ec << "," << s.gamma << "," << s.beta << "," << s.w << "," << s.kc
        << "," << s.tauc << "," << V << "," << Q << "," << pat << "," << s.lambdaBar
        << "," << s.pf << "," << s.qPerfusionIn << "," << s.qPerfusionOut << "\n";
  }

  return 0;
}
