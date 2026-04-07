/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

#include "Rodin/Math/RungeKutta/RK2.h"
#include "Rodin/Math/RungeKutta/RK4.h"
#include "Rodin/Types.h"

using namespace Rodin;
using namespace Rodin::Math::RungeKutta;

class RungeKuttaTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// dy/dt = y, y(0) = 1 → y(t) = e^t
// Single step from t=0 to t=dt
TEST_F(RungeKuttaTest, RK2SingleStepExponential)
{
  RK2 rk2;
  Real q = 0.0;
  Real p = 1.0;
  Real dt = 0.1;
  auto f = [](const Real& y) -> Real { return y; };
  rk2.step(q, dt, p, f);
  // RK2 midpoint: k1=1, k2=f(1+0.05*1)=1.05, q = 1 + 0.1*1.05 = 1.105
  Real exact = std::exp(dt);
  EXPECT_NEAR(q, 1.105, 1e-12);
  EXPECT_NEAR(q, exact, 0.01); // 2nd order: error ~ O(dt^3)
}

TEST_F(RungeKuttaTest, RK4SingleStepExponential)
{
  RK4 rk4;
  Real q = 0.0;
  Real p = 1.0;
  Real dt = 0.1;
  auto f = [](const Real& y) -> Real { return y; };
  rk4.step(q, dt, p, f);
  Real exact = std::exp(dt);
  EXPECT_NEAR(q, exact, 1e-6); // 4th order: very accurate
}

TEST_F(RungeKuttaTest, RK2MultiStepExponential)
{
  RK2 rk2;
  Real y = 1.0;
  Real dt = 0.01;
  int steps = 100; // integrate to t = 1.0
  auto f = [](const Real& y) -> Real { return y; };
  for (int i = 0; i < steps; ++i)
  {
    Real q = 0.0;
    rk2.step(q, dt, y, f);
    y = q;
  }
  Real exact = std::exp(1.0);
  EXPECT_NEAR(y, exact, 0.01);
}

TEST_F(RungeKuttaTest, RK4MultiStepExponential)
{
  RK4 rk4;
  Real y = 1.0;
  Real dt = 0.01;
  int steps = 100; // integrate to t = 1.0
  auto f = [](const Real& y) -> Real { return y; };
  for (int i = 0; i < steps; ++i)
  {
    Real q = 0.0;
    rk4.step(q, dt, y, f);
    y = q;
  }
  Real exact = std::exp(1.0);
  EXPECT_NEAR(y, exact, 1e-8);
}

// dy/dt = -y, y(0) = 1 → y(t) = e^{-t}
TEST_F(RungeKuttaTest, RK4ExponentialDecay)
{
  RK4 rk4;
  Real y = 1.0;
  Real dt = 0.01;
  int steps = 100;
  auto f = [](const Real& y) -> Real { return -y; };
  for (int i = 0; i < steps; ++i)
  {
    Real q = 0.0;
    rk4.step(q, dt, y, f);
    y = q;
  }
  EXPECT_NEAR(y, std::exp(-1.0), 1e-8);
}

TEST_F(RungeKuttaTest, RK2QAndPSameObject)
{
  // Test aliasing: q and p are the same object
  RK2 rk2;
  Real y = 1.0;
  Real dt = 0.1;
  auto f = [](const Real& y) -> Real { return y; };
  rk2.step(y, dt, y, f);
  // Should produce correct result even with aliasing
  EXPECT_NEAR(y, 1.105, 1e-10);
}

TEST_F(RungeKuttaTest, RK4QAndPSameObject)
{
  // Test aliasing: q and p are the same object
  RK4 rk4;
  Real y = 1.0;
  Real dt = 0.1;
  auto f = [](const Real& y) -> Real { return y; };
  rk4.step(y, dt, y, f);
  Real exact = std::exp(0.1);
  EXPECT_NEAR(y, exact, 1e-6);
}

TEST_F(RungeKuttaTest, RK2QAndPDifferentObjects)
{
  RK2 rk2;
  Real p = 1.0;
  Real q = 0.0;
  Real dt = 0.1;
  auto f = [](const Real& y) -> Real { return y; };
  rk2.step(q, dt, p, f);
  EXPECT_NEAR(q, 1.105, 1e-10);
  EXPECT_DOUBLE_EQ(p, 1.0); // p unchanged
}

TEST_F(RungeKuttaTest, RK4QAndPDifferentObjects)
{
  RK4 rk4;
  Real p = 1.0;
  Real q = 0.0;
  Real dt = 0.1;
  auto f = [](const Real& y) -> Real { return y; };
  rk4.step(q, dt, p, f);
  Real exact = std::exp(0.1);
  EXPECT_NEAR(q, exact, 1e-6);
  EXPECT_DOUBLE_EQ(p, 1.0); // p unchanged
}

TEST_F(RungeKuttaTest, RK4ConvergenceOrder)
{
  // Verify 4th order convergence: halving dt should reduce error by ~16x
  RK4 rk4;
  auto f = [](const Real& y) -> Real { return y; };
  Real exact = std::exp(1.0);

  auto integrate = [&](Real dt) -> Real
  {
    Real y = 1.0;
    int steps = static_cast<int>(1.0 / dt + 0.5);
    for (int i = 0; i < steps; ++i)
    {
      Real q = 0.0;
      rk4.step(q, dt, y, f);
      y = q;
    }
    return std::abs(y - exact);
  };

  Real err1 = integrate(0.1);
  Real err2 = integrate(0.05);
  // Ratio should be approximately 2^4 = 16 for 4th order
  Real ratio = err1 / err2;
  EXPECT_GT(ratio, 10.0); // Should be ~16, give some margin
}
