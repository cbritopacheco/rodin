/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <cmath>

#include "Rodin/Math/RootFinding/NewtonRaphson.h"
#include "Rodin/Types.h"

using namespace Rodin;
using namespace Rodin::Math::RootFinding;

class NewtonRaphsonTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(NewtonRaphsonTest, DefaultConstruction)
{
  NewtonRaphson<Real> solver;
  // Should compile and construct with default tolerances
  (void)solver;
}

TEST_F(NewtonRaphsonTest, CustomTolerances)
{
  NewtonRaphson<Real> solver(1e-15, 1e-12, 1e-15, 50);
  (void)solver;
}

TEST_F(NewtonRaphsonTest, SqrtTwo)
{
  // f(x) = x^2 - 2, f'(x) = 2x, root = sqrt(2)
  NewtonRaphson<Real> solver(1e-12, 1e-9, 1e-12, 25);
  auto f = [](Real x) -> std::pair<Real, Real>
  {
    return {x * x - 2.0, 2.0 * x};
  };
  auto root = solver.solve(f, 1.5, 1.0, 2.0);
  ASSERT_TRUE(root.has_value());
  EXPECT_NEAR(*root, std::sqrt(2.0), 1e-10);
}

TEST_F(NewtonRaphsonTest, CubicRoot)
{
  // f(x) = x^3 - 8, f'(x) = 3x^2, root = 2
  NewtonRaphson<Real> solver(1e-12, 1e-9, 1e-12, 25);
  auto f = [](Real x) -> std::pair<Real, Real>
  {
    return {x * x * x - 8.0, 3.0 * x * x};
  };
  auto root = solver.solve(f, 1.5, 1.0, 3.0);
  ASSERT_TRUE(root.has_value());
  EXPECT_NEAR(*root, 2.0, 1e-10);
}

TEST_F(NewtonRaphsonTest, SinRoot)
{
  // f(x) = sin(x), f'(x) = cos(x), root at x = pi in [2, 4]
  NewtonRaphson<Real> solver(1e-12, 1e-9, 1e-12, 25);
  auto f = [](Real x) -> std::pair<Real, Real>
  {
    return {std::sin(x), std::cos(x)};
  };
  auto root = solver.solve(f, 3.0, 2.0, 4.0);
  ASSERT_TRUE(root.has_value());
  EXPECT_NEAR(*root, M_PI, 1e-10);
}

TEST_F(NewtonRaphsonTest, RootNearBoundary)
{
  // f(x) = x - 1.001, f'(x) = 1, root at x = 1.001 (near boundary of [1, 2])
  NewtonRaphson<Real> solver(1e-12, 1e-9, 1e-12, 25);
  auto f = [](Real x) -> std::pair<Real, Real>
  {
    return {x - 1.001, 1.0};
  };
  auto root = solver.solve(f, 1.5, 1.0, 2.0);
  ASSERT_TRUE(root.has_value());
  EXPECT_NEAR(*root, 1.001, 1e-10);
}

TEST_F(NewtonRaphsonTest, MaxIterationLimit)
{
  // Very tight tolerance with only 1 iteration allowed
  NewtonRaphson<Real> solver(1e-30, 1e-30, 1e-30, 1);
  auto f = [](Real x) -> std::pair<Real, Real>
  {
    return {x * x - 2.0, 2.0 * x};
  };
  auto root = solver.solve(f, 1.5, 1.0, 2.0);
  // With only 1 iteration and extremely tight tolerance, may or may not converge
  // but should not crash
  (void)root;
}

TEST_F(NewtonRaphsonTest, LinearFunction)
{
  // f(x) = 2x - 6, f'(x) = 2, root = 3
  NewtonRaphson<Real> solver(1e-12, 1e-9, 1e-12, 25);
  auto f = [](Real x) -> std::pair<Real, Real>
  {
    return {2.0 * x - 6.0, 2.0};
  };
  auto root = solver.solve(f, 2.0, 0.0, 5.0);
  ASSERT_TRUE(root.has_value());
  EXPECT_NEAR(*root, 3.0, 1e-10);
}
