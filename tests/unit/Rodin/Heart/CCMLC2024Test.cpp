/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <cmath>

#include "Rodin/Heart/CCMLC2024.h"

using namespace Rodin;

namespace
{
  Math::Matrix<Real> finiteDifferenceJacobian(
      const Heart::CCMLC2024::Input& input,
      const Math::Vector<Real>& x,
      const Math::Vector<Real>& xPrev,
      Real t,
      Real dt)
  {
    Math::Vector<Real> R0;
    Heart::CCMLC2024::evaluateResidual(input, x, xPrev, t, dt, R0);

    Math::Matrix<Real> J(R0.size(), x.size());
    const Real eps = 1e-7;
    for (Eigen::Index j = 0; j < x.size(); ++j)
    {
      const Real h = eps * std::max<Real>(1.0, std::abs(x[j]));
      auto xp = x;
      auto xm = x;
      xp[j] += h;
      xm[j] -= h;

      Math::Vector<Real> Rp, Rm;
      Heart::CCMLC2024::evaluateResidual(input, xp, xPrev, t, dt, Rp);
      Heart::CCMLC2024::evaluateResidual(input, xm, xPrev, t, dt, Rm);

      J.col(j) = (Rp - Rm) / (2.0 * h);
    }
    return J;
  }
}

TEST(CCMLC2024Test, DerivedJacobianMatchesFiniteDifference)
{
  Heart::CCMLC2024::Input input;
  input.rho = 1.2;
  input.d0 = 0.8;
  input.R0 = 1.1;
  input.eta = 0.05;
  input.Es = 2.0;
  input.mu = 0.9;
  input.alpha = 0.3;
  input.n0 = 0.2;
  input.k0 = 0.4;
  input.sigma0 = 0.35;
  input.w = 0.7;
  input.Cp = 1.4;
  input.Cd = 1.8;
  input.Rp = 1.1;
  input.Rd = 1.6;
  input.Psv = 0.3;
  input.e1D = [](Real C) { return 0.5 * (C * C - 1.0); };
  input.e1DPrime = [](Real C) { return C; };
  input.ubar = [](Real t) { return 0.2 * std::sin(t); };
  input.pAt = [](Real t) { return 0.1 + 0.05 * std::cos(t); };
  input.valveFlow = [](Real pv, Real par, Real pat) { return 0.8 * (pv - par) + 0.1 * pat; };
  input.dValveFlow_dPv = [](Real, Real, Real) { return 0.8; };
  input.dValveFlow_dPar = [](Real, Real, Real) { return -0.8; };
  input.dWe_dJ1 = [](Real j1, Real, Real) { return 0.02 * j1; };
  input.dWe_dJ2 = [](Real, Real j2, Real) { return 0.01 * j2; };
  input.dWe_dJ4 = [](Real, Real, Real j4) { return 0.015 * j4; };
  input.dSigmaPassive_dC = [](Real C, Real d1, Real d2, Real d4)
  {
    (void)d4;
    return 12.0 * std::pow(C, -4) * (d1 + C * d2)
         + 4.0 * (1.0 - std::pow(C, -3)) * d2;
  };

  Math::Vector<Real> x(Heart::CCMLC2024::NVAR), xPrev(Heart::CCMLC2024::NVAR);
  x << 0.08, -0.03, 0.12, 0.07, 0.03, 0.015, 0.02, 0.04;
  xPrev << 0.05, -0.01, 0.10, 0.08, 0.04, 0.012, 0.018, 0.035;

  const Real t = 0.2;
  const Real dt = 1e-2;

  Math::Matrix<Real> Janalytic;
  Heart::CCMLC2024::evaluateJacobian(input, x, xPrev, t, dt, Janalytic);

  const auto Jfd = finiteDifferenceJacobian(input, x, xPrev, t, dt);

  const Real rel = (Janalytic - Jfd).norm() / std::max<Real>(Jfd.norm(), 1e-14);
  EXPECT_LT(rel, 5e-4);
}

TEST(CCMLC2024Test, StepConvergesForQuasiStaticCase)
{
  Heart::CCMLC2024::Input input;
  input.rho = 0.0;
  input.d0 = 1.0;
  input.R0 = 1.0;
  input.eta = 0.0;
  input.Es = 0.0;
  input.mu = 1.0;
  input.alpha = 0.0;
  input.n0 = 0.0;
  input.k0 = 0.0;
  input.sigma0 = 0.0;
  input.w = 1.0;
  input.Cp = 1.0;
  input.Cd = 1.0;
  input.Rp = 1.0;
  input.Rd = 1.0;
  input.Psv = 0.0;
  input.e1D = [](Real) { return 0.0; };
  input.e1DPrime = [](Real) { return 0.0; };
  input.pAt = [](Real) { return 0.0; };
  input.ubar = [](Real) { return 0.0; };
  input.valveFlow = [](Real pv, Real par, Real) { return pv - par; };
  input.dValveFlow_dPv = [](Real, Real, Real) { return 1.0; };
  input.dValveFlow_dPar = [](Real, Real, Real) { return -1.0; };

  Heart::CCMLC2024 model(input);
  model.setMaxIterations(50)
    .setAbsoluteTolerance(1e-10)
    .setRelativeTolerance(1e-10)
    .setStepTolerance(1e-12);

  Heart::CCMLC2024::State initial;
  model.initialize(initial);

  const auto report = model.step(1e-2);
  EXPECT_TRUE(report.converged);
}
