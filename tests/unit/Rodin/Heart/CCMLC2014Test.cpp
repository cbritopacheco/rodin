/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <cmath>
#include <type_traits>

#include "Rodin/Heart/CCMLC2014.h"

using namespace Rodin;

namespace
{
  using Model = Heart::CCMLC2014;

  Math::Matrix<Real> finiteDifferenceJacobian(
      const Model::Input& input,
      const Math::Vector<Real>& x,
      const Model::State& prev,
      Real t,
      Real dt)
  {
    Math::Vector<Real> R0;
    Model::evaluateResidual(input, x, prev, t, dt, R0);

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
      Model::evaluateResidual(input, xp, prev, t, dt, Rp);
      Model::evaluateResidual(input, xm, prev, t, dt, Rm);

      J.col(j) = (Rp - Rm) / (2.0 * h);
    }

    return J;
  }
}

TEST(CCMLC2014Test, DerivedJacobianMatchesFiniteDifference)
{
  Model::Input input;
  input.rho = 1.2;
  input.d0 = 0.8;
  input.R0 = 1.1;
  input.eta = 0.05;
  input.Es = 2.0;
  input.mu = 0.9;
  input.alpha = 0.3;
  input.k0 = 0.4;
  input.sigma0 = 0.35;
  input.Cp = 1.4;
  input.Cd = 1.8;
  input.Rp = 1.1;
  input.Rd = 1.6;
  input.Kat = 0.6;
  input.Kp = 0.4;
  input.Kar = 0.9;

  input.u = [](Real t) { return 0.2 * std::sin(t); };
  input.pAt = [](Real t) { return 0.1 + 0.05 * std::cos(t); };
  input.pSv = [](Real) { return 0.3; };
  input.n0 = [](Real) { return 0.9; };
  input.m0 = [](Real ec) { return 1.0 + 0.2 * ec; };
  input.m0Prime = [](Real) { return 0.2; };

  {
    using PassiveEnergy = std::decay_t<decltype(input.passiveEnergy)>;
    typename PassiveEnergy::Parameters hp;
    hp.mu1 = 0.0;
    hp.mu2 = 0.01;
    hp.C0 = 0.05;
    hp.C1 = 0.5;
    hp.C2 = 0.03;
    hp.C3 = 0.6;
    input.passiveEnergy = PassiveEnergy(hp);
  }

  Math::Vector<Real> x(Model::NVAR);
  x << 0.08, -0.03, 0.12, 0.07, 0.03, 0.015, 0.02, 0.04, 1.1;

  Model::State prev;
  prev.y = 0.05;
  prev.v = -0.01;
  prev.pv = 0.10;
  prev.par = 0.08;
  prev.pd = 0.04;
  prev.ec = 0.012;
  prev.kc = 0.018;
  prev.tauc = 0.035;
  prev.w = 1.0;
  prev.t = 0.2;

  const Real t = 0.2;
  const Real dt = 1e-2;

  Math::Matrix<Real> Janalytic;
  Model::evaluateJacobian(input, x, prev, t, dt, Janalytic);

  const auto Jfd = finiteDifferenceJacobian(input, x, prev, t, dt);

  const Real rel =
    (Janalytic - Jfd).norm() / std::max<Real>(Jfd.norm(), 1e-14);

  EXPECT_LT(rel, 5e-4);
}

TEST(CCMLC2014Test, StepConvergesForQuasiStaticCase)
{
  Model::Input input;
  input.rho = 0.0;
  input.d0 = 1.0;
  input.R0 = 1.0;
  input.eta = 0.0;
  input.Es = 0.0;
  input.mu = 1.0;
  input.alpha = 0.0;
  input.k0 = 0.0;
  input.sigma0 = 0.0;
  input.Cp = 1.0;
  input.Cd = 1.0;
  input.Rp = 1.0;
  input.Rd = 1.0;
  input.Kat = 1.0;
  input.Kp = 1.0;
  input.Kar = 1.0;

  input.pAt = [](Real) { return 0.0; };
  input.pSv = [](Real) { return 0.0; };
  input.u = [](Real) { return 0.0; };
  input.n0 = [](Real) { return 0.0; };
  input.m0 = [](Real) { return 1.0; };
  input.m0Prime = [](Real) { return 0.0; };

  Model model(input);
  model.setMaxIterations(50)
    .setAbsoluteTolerance(1e-10)
    .setRelativeTolerance(1e-10)
    .setStepTolerance(1e-12);

  Model::State initial;
  initial.w = 1.0;
  model.initialize(initial);

  const auto report = model.step(1e-2);
  EXPECT_TRUE(report.converged);
}
