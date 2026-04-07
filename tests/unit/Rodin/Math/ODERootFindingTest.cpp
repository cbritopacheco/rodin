#include <gtest/gtest.h>

#include <cmath>

#include "Rodin/Math/RungeKutta/RK2.h"
#include "Rodin/Math/RungeKutta/RK4.h"
#include "Rodin/Math/RootFinding/NewtonRaphson.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Types.h"

using namespace Rodin;
using namespace Rodin::Math;

class ODERootFindingTest : public ::testing::Test
{};

TEST_F(ODERootFindingTest, RK2ExponentialGrowth)
{
  RungeKutta::RK2 rk2;
  Vector<Real> p(1);
  p(0) = 1.0;
  Vector<Real> q = p;
  const Real dt = 0.1;

  auto f = [](const auto& state)
  {
    return state;
  };

  rk2.step(q, dt, p, f);

  EXPECT_NEAR(q(0), std::exp(dt), 2e-2);
}

TEST_F(ODERootFindingTest, RK4ExponentialGrowth)
{
  RungeKutta::RK4 rk4;
  Vector<Real> p(1);
  p(0) = 1.0;
  Vector<Real> q(1);
  q(0) = 0.0;
  const Real dt = 0.1;

  auto f = [](const auto& state)
  {
    return state;
  };

  rk4.step(q, dt, p, f);

  EXPECT_NEAR(q(0), std::exp(dt), 1e-6);
}

TEST_F(ODERootFindingTest, NewtonRaphsonFindsSqrt2)
{
  RootFinding::NewtonRaphson<Real> solver(1e-12, 1e-9, 1e-12, 64);
  auto f = [](Real x)
  {
    return std::make_pair(x * x - 2.0, 2.0 * x);
  };

  const auto root = solver.solve(f, 1.5, 1.0, 2.0);
  ASSERT_TRUE(root.has_value());
  EXPECT_NEAR(*root, std::sqrt(2.0), 1e-10);
}

TEST_F(ODERootFindingTest, NewtonRaphsonNaNDerivativeFallback)
{
  RootFinding::NewtonRaphson<Real> solver(1e-12, 1e-9, 1e-12, 128);
  auto f = [](Real x)
  {
    return std::make_pair(x * x - 2.0, std::numeric_limits<Real>::quiet_NaN());
  };

  const auto root = solver.solve(f, 1.2, 1.0, 2.0);
  ASSERT_TRUE(root.has_value());
  EXPECT_NEAR(*root, std::sqrt(2.0), 1e-8);
}

TEST_F(ODERootFindingTest, NewtonRaphsonNearZeroDerivative)
{
  RootFinding::NewtonRaphson<Real> solver(1e-12, 1e-9, 1e-12, 256);
  auto f = [](Real x)
  {
    const Real d = x - 1.0;
    return std::make_pair(d * d * d, 3.0 * d * d);
  };

  const auto root = solver.solve(f, 0.2, 0.0, 2.0);
  ASSERT_TRUE(root.has_value());
  EXPECT_NEAR(*root, 1.0, 1e-8);
}

TEST_F(ODERootFindingTest, NewtonRaphsonMaxIterationReturnsEmpty)
{
  RootFinding::NewtonRaphson<Real> solver(1e-20, 1e-20, 1e-20, 0);
  auto f = [](Real x)
  {
    return std::make_pair(x - 0.123456789, 1.0);
  };

  const auto root = solver.solve(f, 0.9, 0.0, 1.0);
  EXPECT_FALSE(root.has_value());
}
