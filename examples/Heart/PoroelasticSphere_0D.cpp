/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file PoroelasticSphere_0D.cpp
 * @brief Standalone driver for the reduced poroelastic left ventricle.
 *
 * Integrates @ref Rodin::Heart::PoroelasticSphereT (shell reduction,
 * Usyk-type passive energy, CCMLC2014 activation, valves and Windkessel,
 * exponential pore storage, shear-dependent venous drainage) under the
 * standalone circulatory closure: the myocardial supply is the lumped
 * conductance gammaAr drawn from the proximal Windkessel node.
 *
 * Usage: PoroelasticSphere_0D [scenario] [cycles] [csv]
 *   scenario: healthy (default), aortic_stenosis, hyperviscosity, lvh.
 */
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <string>

#include "CoronaryArtery/PoroelasticLV.h"

using Real = Rodin::Real;
using namespace Rodin::Examples::Heart;
using Model = PoroelasticLV;

int main(int argc, char** argv)
{
  const std::string name = (argc > 1) ? argv[1] : "healthy";
  const int cycles = (argc > 2) ? std::max(1, std::atoi(argv[2])) : 10;
  const std::string csv =
    (argc > 3) ? argv[3] : ("poroelastic_sphere_0d_" + name + ".csv");

  const Scenario scenario = Scenario::make(name);
  // Standalone closure: lumped supply of the order of the resolved LCA.
  const Model::Input in = makePoroelasticLVInput(scenario, 2.0e-10);

  Model model(in);
  model.setMaxIterations(200)
    .setAbsoluteTolerance(1e-8)
    .setRelativeTolerance(1e-8)
    .setStepTolerance(1e-10)
    .setDampingFactor(1.0);
  model.initialize(makePoroelasticLVInitialState(in));

  const Real dt = 1e-3;
  const int stepsPerCycle = static_cast<int>(std::lround(0.85 / dt));
  const int nsteps = cycles * stepsPerCycle;
  const Real pi = std::numbers::pi_v<Real>;
  Real Vprev = (4.0 / 3.0) * pi * in.R0 * in.R0 * in.R0;

  std::ofstream out(csv);
  out << std::setprecision(12);
  out << "t,y,y_epi,h,phi,pv,par,pd,ec,kc,tauc,w,V,Q,pat,lambdaBar,pP,pf,Qar,Qven,"
         "shear_rate,nu_eff\n";

  for (int i = 0; i < nsteps; ++i)
  {
    const auto rep = model.step(dt);
    if (!rep.converged)
    {
      std::cerr << "0D step " << i << " did not converge (t = " << model.getState().t
                << ", residual " << rep.finalResidual << ")\n";
      return 1;
    }
    const auto& s = model.getState();
    const Real R = in.R0 + s.y;
    const Real V = (4.0 / 3.0) * pi * R * R * R;
    const Real Q = (V - Vprev) / dt;
    Vprev = V;
    out << s.t << ',' << s.y << ',' << s.yEpi << ',' << s.thickness << ',' << s.phi << ','
        << s.pv << ',' << s.par << ',' << s.pd << ',' << s.ec << ',' << s.kc << ','
        << s.tauc << ',' << s.w << ',' << V << ',' << Q << ',' << in.pAt(s.t) << ','
        << s.lambdaBar << ',' << s.pP << ',' << s.pf << ',' << s.qArterial << ','
        << s.qVenous << ',' << s.shearRate << ',' << s.poreViscosity << '\n';
    if ((i + 1) % stepsPerCycle == 0)
      std::cout << "[" << scenario.name << "] cycle " << (i + 1) / stepsPerCycle
                << "  p_v=" << s.pv << "  p_ar=" << s.par << "  p_f=" << s.pf
                << "  phi=" << s.phi << '\n';
  }
  return 0;
}
