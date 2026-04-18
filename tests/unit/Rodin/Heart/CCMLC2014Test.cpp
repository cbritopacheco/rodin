/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>

#include "Rodin/Heart/CCMLC2014.h"
#include "Rodin/Heart/CCMLC2014/Numerics/Jacobian.h"
#include "Rodin/Heart/CCMLC2014/Numerics/Residual.h"

using namespace Rodin;

namespace
{
  using Model = Heart::CCMLC2014T<>;
  using PassiveLaw = Heart::CCMLC2014PassiveLaw<Real>;
  namespace ModelVars = Heart::CCMLC2014::Model;

  Model::Input makeInput()
  {
    Model::Input in;

    in.rho = 1.0e3;
    in.R0 = 2.36e-2;
    in.d0 = 1.42e-2;

    in.Es = 3.0e5;
    in.mu = 70.0;
    in.eta = 70.0;
    in.alpha = 3.0;
    in.k0 = 1.0e5;
    in.sigma0 = 5.0e5;

    in.Rp = 8.0e6;
    in.Cp = 5.0e-9;
    in.Rd = 1.0e8;
    in.Cd = 1.0e-8;

    in.Kat = 8.0e-7;
    in.Kp  = 5.0e-10;
    in.Kar = 1.3e-5;

    in.cavityCapacity = 5.0e-12;
    in.localTolerance = 1e-12;
    in.localMaxIterations = 50;
    in.localDamping = 1.0;
    in.absRegularization = 1e-14;

    in.u = [](Real) { return 25.0; };
    in.pAt = [](Real) { return 900.0; };
    in.pSv = [](Real) { return 1000.0; };

    typename std::decay_t<decltype(in.passiveEnergy)>::Parameters hp;
    hp.mu1 = 0.0;
    hp.mu2 = 0.0;
    hp.C0 = 1.9e3;
    hp.C1 = 1.1e-1;
    hp.C2 = 1.9e3;
    hp.C3 = 1.1e-1;
    in.passiveEnergy = std::decay_t<decltype(in.passiveEnergy)>(hp);

    return in;
  }
}

TEST(CCMLC2014Test, InitializeUsesInputActiveDefaultsWhenInitialActiveStateIsZero)
{
  auto in = makeInput();
  in.initFibDef = 0.12;
  in.initActiveStiffness = 0.25;
  in.initActiveStress = 0.5;

  Model model(in);

  Model::State initial;
  initial.t = 0.0;
  model.initialize(initial);

  const auto& s = model.getState();
  EXPECT_NEAR(s.ec, 0.12, 1e-14);
  EXPECT_NEAR(s.gamma, 0.5, 1e-14);
  EXPECT_NEAR(s.beta, 1.0, 1e-14);
  EXPECT_NEAR(s.kc, 0.25, 1e-14);
  EXPECT_NEAR(s.tauc, 0.5, 1e-14);
}

TEST(CCMLC2014Test, InitializeUsesProvidedGammaAndBeta)
{
  auto in = makeInput();

  Model model(in);

  Model::State initial;
  initial.t = 0.0;
  initial.ec = 0.03;
  initial.gamma = 0.4;
  initial.beta = 0.2;
  initial.kc = 1.0;
  initial.tauc = 1.0;
  model.initialize(initial);

  const auto& s = model.getState();
  EXPECT_NEAR(s.ec, 0.03, 1e-14);
  EXPECT_NEAR(s.gamma, 0.4, 1e-14);
  EXPECT_NEAR(s.beta, 0.2, 1e-14);
  EXPECT_NEAR(s.kc, 0.16, 1e-14);
  EXPECT_NEAR(s.tauc, 0.08, 1e-14);
}

TEST(CCMLC2014Test, StepConvergesAndAdvancesTime)
{
  auto in = makeInput();
  in.initFibDef = 0.0;
  in.initActiveStiffness = 0.0;
  in.initActiveStress = 0.0;

  Model model(in);
  model.setMaxIterations(200)
       .setAbsoluteTolerance(1e-8)
       .setRelativeTolerance(1e-8)
       .setStepTolerance(1e-10)
       .setDampingFactor(1.0);

  Model::State initial;
  initial.t = 0.0;
  initial.y = 0.0;
  initial.v = 0.0;
  initial.pv = in.pAt(0.0) - 100.0;
  initial.par = 11000.0;
  initial.pd = 10000.0;
  model.initialize(initial);

  const Real dt = 1e-3;
  const auto rep = model.step(dt);
  EXPECT_TRUE(rep.converged);
  EXPECT_NEAR(model.getState().t, dt, 1e-14);
}

TEST(CCMLC2014Test, DynamicJacobianMatchesFiniteDifference)
{
  auto in = makeInput();

  Model::DenseVector x(ModelVars::NVAR);
  x[ModelVars::DISP] = 8e-5;
  x[ModelVars::PV] = 1.0e4;
  x[ModelVars::PAR] = 9.0e3;
  x[ModelVars::PD] = 8.0e3;

  Model::State sn;
  sn.t = 0.1;
  sn.y = 7e-5;
  sn.pv = 9.8e3;
  sn.par = 8.8e3;
  sn.pd = 7.9e3;
  sn.ec = 0.01;
  sn.gamma = 0.2;
  sn.beta = 0.3;
  sn.kc = sn.gamma * sn.gamma;
  sn.tauc = sn.gamma * sn.beta;

  Model::State snm1 = sn;
  snm1.t = 0.099;
  snm1.y = 6e-5;

  const Real dt = 1e-3;
  const Real tnp1 = sn.t + dt;

  Model::EvalData d;
  Heart::CCMLC2014::Numerics::buildEvalData<PassiveLaw>(in, x, sn, snm1, tnp1, dt, d);

  Model::DenseMatrix J;
  Heart::CCMLC2014::Numerics::evaluateDynamicJacobian(in, J, d, dt);

  Model::DenseMatrix Jfd(ModelVars::NVAR, ModelVars::NVAR);
  Jfd.setZero();

  const Real eps = 1e-7;
  for (Index j = 0; j < ModelVars::NVAR; ++j)
  {
    const Real h = eps * std::max<Real>(1.0, std::abs(x[j]));
    auto xp = x;
    auto xm = x;
    xp[j] += h;
    xm[j] -= h;

    Model::EvalData dp;
    Model::EvalData dm;
    Heart::CCMLC2014::Numerics::buildEvalData<PassiveLaw>(in, xp, sn, snm1, tnp1, dt, dp);
    Heart::CCMLC2014::Numerics::buildEvalData<PassiveLaw>(in, xm, sn, snm1, tnp1, dt, dm);

    Model::DenseVector Rp;
    Model::DenseVector Rm;
    Heart::CCMLC2014::Numerics::evaluateDynamicResidual(in, dp, Rp);
    Heart::CCMLC2014::Numerics::evaluateDynamicResidual(in, dm, Rm);
    Jfd.col(j) = (Rp - Rm) / (2.0 * h);
  }

  const Real rel = (J - Jfd).norm() / std::max<Real>(Jfd.norm(), 1e-14);
  EXPECT_LT(rel, 1e-3);
}
