/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numbers>
#include <string>

#include "Rodin/Heart/CCMLC2014.h"
#include "Rodin/Heart/CCMLC2014/Physics/Windkessel.h"

using namespace Rodin;

namespace
{
  struct BranchGeometry
  {
      const char* name;
      Real radius;
      Real length;
  };

  using Model = Heart::CCMLC2014T<>;

  Model::Input makeInput()
  {
    Model::Input input;
    input.mu_0 = 0.353;
    input.mu_Inf = 0.004181;
    input.lambda = 15.6821;
    input.n = 0.2050;
    input.yasuda = 0.6497;
    return input;
  }

  Real newtonianResistance(Real mu, Real length, Real radius)
  {
    return 8.0 * mu * length / (std::numbers::pi_v<Real> * std::pow(radius, 4.0));
  }

  void writeSample(std::ostream& out, const char* family, const BranchGeometry& geometry,
    Real pressureDrop, const Model::Input& input)
  {
    using Law = Heart::CCMLC2014::Physics::Rheology::CarreauYasuda;
    const Real r0 = newtonianResistance(input.mu_0, geometry.length, geometry.radius);
    const Real rInf = newtonianResistance(input.mu_Inf, geometry.length, geometry.radius);
    const auto d =
      Law::diagnostics(pressureDrop, geometry.length, geometry.radius, r0, input);

    out << family << ',' << geometry.name << ',' << geometry.radius << ','
        << geometry.length << ',' << pressureDrop << ',' << d.tauW << ',' << d.gammaW
        << ',' << d.carreauNumber << ',' << d.wallViscosity << ',' << d.eta << ','
        << (d.wallViscosity / input.mu_Inf) << ',' << d.flow << ','
        << d.dFlow_dPressureDrop << ',' << d.effectiveResistance << ','
        << d.lowShearNewtonianFlow << ',' << d.highShearNewtonianFlow << ','
        << (d.effectiveResistance / rInf) << ',' << d.usedPoiseuilleFallback << '\n';
  }
}

int main(int argc, char** argv)
{
  const std::string path = (argc > 1) ? argv[1] : "CoronaryOutletFlowDiagnostic.csv";
  std::ofstream out(path);
  if (!out)
  {
    std::cerr << "Failed to open " << path << " for writing.\n";
    return 1;
  }

  const auto input = makeInput();

  out << "family,geometry,radius_m,length_m,pressure_drop_pa,tau_w_pa,gamma_w_per_s,"
      << "carreau_number,mu_w_pa_s,eta,mu_over_mu_inf,q_m3_per_s,dq_dpressure,"
      << "effective_resistance_pa_s_per_m3,q_newtonian_mu0,q_newtonian_mu_inf,"
      << "effective_resistance_over_r_inf,used_poiseuille_fallback\n";

  const std::array<BranchGeometry, 4> defaultGeometries{{
    {"proximal_0p4mm", 4.0e-4, 1.0e-2},
    {"proximal_0p6mm", 6.0e-4, 1.25e-2},
    {"distal_0p2mm", 2.0e-4, 2.5e-3},
    {"distal_0p3mm", 3.0e-4, 2.5e-3},
  }};

  const std::array<Real, 6> radii{{1.0e-4, 2.0e-4, 3.0e-4, 4.0e-4, 6.0e-4, 1.0e-3}};

  const Real minDp = 1.0e-3;
  const Real maxDp = 1.0e6;
  const int samples = 120;
  for (int i = 0; i < samples; ++i)
  {
    const Real s = static_cast<Real>(i) / static_cast<Real>(samples - 1);
    const Real dp = minDp * std::pow(maxDp / minDp, s);

    for (const auto& geometry : defaultGeometries)
      writeSample(out, "default_geometry", geometry, dp, input);

    for (const Real radius : radii)
      writeSample(
        out, "fixed_length_radius_sweep", {"radius_sweep", radius, 2.5e-3}, dp, input);

    writeSample(
      out, "constant_newtonian_resistance", {"base", 2.0e-4, 2.5e-3}, dp, input);
    writeSample(
      out, "constant_newtonian_resistance", {"double_radius", 4.0e-4, 4.0e-2}, dp, input);
    writeSample(out, "constant_newtonian_resistance", {"half_radius", 1.0e-4, 1.5625e-4},
      dp, input);
  }

  std::cout << "Wrote " << path << '\n';
  return 0;
}
