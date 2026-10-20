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

/// @brief Verifies default construction for newton raphson test by checking solver behavior.
TEST_F(NewtonRaphsonTest, DefaultConstruction)
{
  NewtonRaphson<Real> solver;
  // Should compile and construct with default tolerances
  (void)solver;
}

/// @brief Verifies custom tolerances for newton raphson test by checking solver behavior.
TEST_F(NewtonRaphsonTest, CustomTolerances)
{
  NewtonRaphson<Real> solver(1e-15, 1e-12, 1e-15, 50);
  (void)solver;
}

/// @brief F(x) = x^2 - 2, f'(x) = 2x, root = sqrt(2).
TEST_F(NewtonRaphsonTest, SqrtTwo)
{
  // f(x) = x^2 - 2, f'(x) = 2x, root = sqrt(2)
  NewtonRaphson<Real> solver(1e-12, 1e-9, 1e-12, 25);
  auto f = [](Real x) -> std::pair<Real, Real> { return {x * x - 2.0, 2.0 * x}; };
  auto root = solver.solve(f, 1.5, 1.0, 2.0);
  ASSERT_TRUE(root.has_value());
  EXPECT_NEAR(*root, std::sqrt(2.0), 1e-10);
}

/// @brief F(x) = x^3 - 8, f'(x) = 3x^2, root = 2.
TEST_F(NewtonRaphsonTest, CubicRoot)
{
  // f(x) = x^3 - 8, f'(x) = 3x^2, root = 2
  NewtonRaphson<Real> solver(1e-12, 1e-9, 1e-12, 25);
  auto f = [](Real x) -> std::pair<Real, Real> { return {x * x * x - 8.0, 3.0 * x * x}; };
  auto root = solver.solve(f, 1.5, 1.0, 3.0);
  ASSERT_TRUE(root.has_value());
  EXPECT_NEAR(*root, 2.0, 1e-10);
}

/// @brief F(x) = sin(x), f'(x) = cos(x), root at x = pi in [2, 4].
TEST_F(NewtonRaphsonTest, SinRoot)
{
  // f(x) = sin(x), f'(x) = cos(x), root at x = pi in [2, 4]
  NewtonRaphson<Real> solver(1e-12, 1e-9, 1e-12, 25);
  auto f = [](Real x) -> std::pair<Real, Real> { return {std::sin(x), std::cos(x)}; };
  auto root = solver.solve(f, 3.0, 2.0, 4.0);
  ASSERT_TRUE(root.has_value());
  EXPECT_NEAR(*root, M_PI, 1e-10);
}

/// @brief F(x) = x - 1.001, f'(x) = 1, root at x = 1.001 (near boundary of [1, 2]).
TEST_F(NewtonRaphsonTest, RootNearBoundary)
{
  // f(x) = x - 1.001, f'(x) = 1, root at x = 1.001 (near boundary of [1, 2])
  NewtonRaphson<Real> solver(1e-12, 1e-9, 1e-12, 25);
  auto f = [](Real x) -> std::pair<Real, Real> { return {x - 1.001, 1.0}; };
  auto root = solver.solve(f, 1.5, 1.0, 2.0);
  ASSERT_TRUE(root.has_value());
  EXPECT_NEAR(*root, 1.001, 1e-10);
}

/// @brief Very tight tolerance with only 1 iteration allowed.
TEST_F(NewtonRaphsonTest, MaxIterationLimit)
{
  // Very tight tolerance with only 1 iteration allowed
  NewtonRaphson<Real> solver(1e-30, 1e-30, 1e-30, 1);
  auto f = [](Real x) -> std::pair<Real, Real> { return {x * x - 2.0, 2.0 * x}; };
  auto root = solver.solve(f, 1.5, 1.0, 2.0);
  // With only 1 iteration and extremely tight tolerance, may or may not converge
  // but should not crash
  (void)root;
}

/// @brief F(x) = 2x - 6, f'(x) = 2, root = 3.
TEST_F(NewtonRaphsonTest, LinearFunction)
{
  // f(x) = 2x - 6, f'(x) = 2, root = 3
  NewtonRaphson<Real> solver(1e-12, 1e-9, 1e-12, 25);
  auto f = [](Real x) -> std::pair<Real, Real> { return {2.0 * x - 6.0, 2.0}; };
  auto root = solver.solve(f, 2.0, 0.0, 5.0);
  ASSERT_TRUE(root.has_value());
  EXPECT_NEAR(*root, 3.0, 1e-10);
}

namespace
{
  template <class Scalar>
  class NewtonSafeguards : public ::testing::Test
  {};
  using RootScalars = ::testing::Types<Float, Real>;
  TYPED_TEST_SUITE(NewtonSafeguards, RootScalars);

  /// @brief Exact bracket endpoints are returned without evaluating an interior point.
  TYPED_TEST(NewtonSafeguards, ExactEndpoints)
  {
    using Scalar = TypeParam;
    const Scalar lo = Scalar(-2), hi = Scalar(3);
    for (const Scalar expected : {lo, hi})
    {
      for (const Scalar slope : {Scalar(-3), Scalar(3)})
      {
        size_t calls = 0;
        const auto function = [&](Scalar x) {
          ++calls;
          EXPECT_TRUE(x == lo || x == hi);
          return std::pair{(x - expected) * slope, slope};
        };
        const auto root = NewtonRaphson<Scalar>().solve(function, Scalar(0), lo, hi);
        ASSERT_TRUE(root);
        EXPECT_EQ(*root, expected);
        EXPECT_EQ(calls, 2u);
      }
    }
  }

  /// @brief Nonfinite and out-of-bracket initial guesses are replaced by a valid midpoint.
  TYPED_TEST(NewtonSafeguards, InvalidInitialGuess)
  {
    using Scalar = TypeParam;
    const auto function = [](Scalar x) { return std::pair{x - Scalar(0.25), Scalar(1)}; };
    for (const Scalar guess : {Scalar(-10), Scalar(0), Scalar(1), Scalar(10),
           std::numeric_limits<Scalar>::quiet_NaN()})
    {
      const auto root =
        NewtonRaphson<Scalar>().solve(function, guess, Scalar(0), Scalar(1));
      ASSERT_TRUE(root);
      EXPECT_EQ(*root, Scalar(0.25));
    }
  }

  /// @brief Missing or zero derivatives still converge with a valid bracket.
  TYPED_TEST(NewtonSafeguards, DerivativeSafeguards)
  {
    using Scalar = TypeParam;
    constexpr Scalar Tolerance = Scalar(64) * std::numeric_limits<Scalar>::epsilon();
    const NewtonRaphson<Scalar> solver(Tolerance, Tolerance, Tolerance, 100);
    for (const Scalar derivative : {Scalar(0), std::numeric_limits<Scalar>::quiet_NaN()})
    {
      const auto function = [derivative](Scalar x) {
        return std::pair{x * x - Scalar(2), derivative};
      };
      const auto root = solver.solve(function, Scalar(1.5), Scalar(1), Scalar(2));
      ASSERT_TRUE(root);
      EXPECT_GE(*root, Scalar(1));
      EXPECT_LE(*root, Scalar(2));
      EXPECT_NEAR(*root, std::sqrt(Scalar(2)), Scalar(4) * Tolerance);
      EXPECT_LE(std::abs(function(*root).first), Scalar(8) * Tolerance);
    }
  }

  /// @brief Rescaling a finite linear function preserves its bracket and analytic root.
  TYPED_TEST(NewtonSafeguards, FiniteFunctionScalingPreservesBracket)
  {
    using Scalar = TypeParam;
    constexpr Scalar CoordinateTolerance =
      Scalar(64) * std::numeric_limits<Scalar>::epsilon();
    // Covers all binary64 bisection levels plus safeguarded startup steps.
    constexpr size_t WorkBudget = 100;
    constexpr Scalar Exact = Scalar(0.75);
    // The small scale keeps evaluations finite while their products underflow.
    for (const Scalar scale : {Scalar(1), Scalar(4) * std::numeric_limits<Scalar>::min(),
           std::numeric_limits<Scalar>::max() / Scalar(4)})
    {
      for (const Scalar sign : {Scalar(-1), Scalar(1)})
      {
        for (const Scalar initial : {Scalar(0.125), Scalar(0.875)})
        {
          SCOPED_TRACE(::testing::Message()
            << "scale=" << scale << " sign=" << sign << " initial=" << initial);
          const Scalar functionTolerance = scale * CoordinateTolerance;
          const auto function = [&](Scalar x) {
            EXPECT_GE(x, Scalar(0));
            EXPECT_LE(x, Scalar(1));
            return std::pair{sign * scale * (x - Exact), sign * scale};
          };
          const NewtonRaphson<Scalar> solver(
            CoordinateTolerance, Scalar(0), functionTolerance, WorkBudget);
          const auto root = solver.solve(function, initial, Scalar(0), Scalar(1));
          ASSERT_TRUE(root);
          EXPECT_NEAR(*root, Exact, CoordinateTolerance);
          EXPECT_LE(std::abs(function(*root).first), functionTolerance);
        }
      }
    }
  }

  /// @brief Exact interior roots terminate even when the function tolerance is zero.
  TYPED_TEST(NewtonSafeguards, ExactInteriorRootWithZeroFunctionTolerance)
  {
    using Scalar = TypeParam;
    constexpr Scalar CoordinateTolerance =
      Scalar(64) * std::numeric_limits<Scalar>::epsilon();
    // Covers all binary64 bisection levels plus safeguarded startup steps.
    constexpr size_t WorkBudget = 100;
    constexpr Scalar Exact = Scalar(0.75);
    const auto function = [](Scalar x) { return std::pair{x - Exact, Scalar(1)}; };
    const NewtonRaphson<Scalar> solver(
      CoordinateTolerance, Scalar(0), Scalar(0), WorkBudget);
    for (const Scalar initial : {Exact, Scalar(0.875)})
    {
      const auto root = solver.solve(function, initial, Scalar(0), Scalar(1));
      ASSERT_TRUE(root);
      EXPECT_EQ(*root, Exact);
    }
  }

  /// @brief Finite large brackets keep every midpoint finite and inside the supplied bounds.
  TYPED_TEST(NewtonSafeguards, LargeFiniteBoundsHaveFiniteMidpoints)
  {
    using Scalar = TypeParam;
    constexpr Scalar Magnitude = std::numeric_limits<Scalar>::max();
    constexpr Scalar Tolerance =
      Magnitude * (Scalar(64) * std::numeric_limits<Scalar>::epsilon());
    // Covers all binary64 bisection levels plus safeguarded startup steps.
    constexpr size_t WorkBudget = 100;
    for (const Scalar sign : {Scalar(-1), Scalar(1)})
    {
      const Scalar lower = sign > 0 ? Magnitude / Scalar(2) : -Magnitude;
      const Scalar upper = sign > 0 ? Magnitude : -Magnitude / Scalar(2);
      const Scalar exact = sign * (Magnitude / Scalar(4)) * Scalar(3);
      const Scalar interior = lower + (upper - lower) / Scalar(4);
      for (const bool coarse : {false, true})
      {
        for (const Scalar initial : {Scalar(0), interior})
        {
          SCOPED_TRACE(::testing::Message()
            << "sign=" << sign << " coarse=" << coarse << " initial=" << initial);
          const auto function = [&](Scalar x) {
            EXPECT_TRUE(std::isfinite(x));
            EXPECT_GE(x, lower);
            EXPECT_LE(x, upper);
            return std::pair{x - exact, std::numeric_limits<Scalar>::quiet_NaN()};
          };
          const NewtonRaphson<Scalar> solver(
            coarse ? Magnitude : Tolerance, Scalar(0), Tolerance, WorkBudget);
          const auto root = solver.solve(function, initial, lower, upper);
          ASSERT_TRUE(root);
          EXPECT_TRUE(std::isfinite(*root));
          EXPECT_GE(*root, lower);
          EXPECT_LE(*root, upper);
          EXPECT_NEAR(*root, exact, Tolerance);
          EXPECT_LE(std::abs(function(*root).first), Tolerance);
        }
      }
    }
  }

  /// @brief A bracket narrower than the coordinate tolerance returns its midpoint and error bound.
  TYPED_TEST(NewtonSafeguards, IntervalToleranceReturnsMidpoint)
  {
    using Scalar = TypeParam;
    constexpr Scalar Width = Scalar(1) / Scalar(128);
    constexpr Scalar Lower = Scalar(1), Upper = Lower + Width;
    constexpr Scalar Exact = Lower + Scalar(3) * Width / Scalar(4);
    constexpr Scalar Initial = Lower + Width / Scalar(4);
    constexpr Scalar AbsoluteIntervalTolerance = Scalar(2) * Width;
    constexpr Scalar FunctionTolerance = std::numeric_limits<Scalar>::epsilon();
    size_t calls = 0;
    const auto function = [&](Scalar x) {
      ++calls;
      EXPECT_GE(x, Lower);
      EXPECT_LE(x, Upper);
      return std::pair{x - Exact, Scalar(1)};
    };
    const NewtonRaphson<Scalar> solver(
      AbsoluteIntervalTolerance, Scalar(0), FunctionTolerance);
    const auto root = solver.solve(function, Initial, Lower, Upper);
    ASSERT_TRUE(root);
    EXPECT_EQ(*root, (Lower + Upper) / Scalar(2));
    EXPECT_LE(std::abs(*root - Exact), Width / Scalar(2));
    EXPECT_EQ(calls, 3u);
  }

  /// @brief A depleted work budget returns no unconverged interior result.
  TYPED_TEST(NewtonSafeguards, ZeroWorkBudget)
  {
    using Scalar = TypeParam;
    const auto function = [](Scalar x) {
      return std::pair{x * x - Scalar(2), Scalar(2) * x};
    };
    const NewtonRaphson<Scalar> solver(Scalar(0), Scalar(0), Scalar(0), 0);
    EXPECT_FALSE(solver.solve(function, Scalar(1.5), Scalar(1), Scalar(2)));
  }
}
