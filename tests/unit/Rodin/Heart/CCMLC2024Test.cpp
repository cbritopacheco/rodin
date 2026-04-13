/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Heart/CCMLC2024.h"

using namespace Rodin;

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
  input.pAt = [](Real) { return 0.0; };
  input.ubar = [](Real) { return 0.0; };
  input.valveFlow = [](Real pv, Real par, Real) { return pv - par; };

  Heart::CCMLC2024 model(input);
  model.setMaxIterations(50)
    .setAbsoluteTolerance(1e-10)
    .setRelativeTolerance(1e-10)
    .setStepTolerance(1e-12);

  Heart::CCMLC2024::State initial;
  initial.y = 0.0;
  initial.v = 0.0;
  initial.pv = 0.0;
  initial.par = 0.0;
  initial.pd = 0.0;
  initial.ec = 0.0;
  initial.kc = 0.0;
  initial.tauc = 0.0;
  initial.t = 0.0;
  model.initialize(initial);

  const auto report = model.step(1e-2);
  EXPECT_TRUE(report.converged);

  const auto& s = model.getState();
  EXPECT_NEAR(s.y, 0.0, 1e-12);
  EXPECT_NEAR(s.v, 0.0, 1e-12);
  EXPECT_NEAR(s.pv, 0.0, 1e-12);
  EXPECT_NEAR(s.par, 0.0, 1e-12);
  EXPECT_NEAR(s.pd, 0.0, 1e-12);
  EXPECT_NEAR(s.ec, 0.0, 1e-12);
  EXPECT_NEAR(s.kc, 0.0, 1e-12);
  EXPECT_NEAR(s.tauc, 0.0, 1e-12);
}
