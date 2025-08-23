/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <cmath>

#include "Rodin/Math/Rad.h"
#include "Rodin/Types.h"

using namespace Rodin;
using namespace Rodin::Math;

class RadTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// Test Rad construction
TEST_F(RadTest, Construction)
{
  // Test value constructor
  Rad r1(M_PI);
  EXPECT_DOUBLE_EQ(static_cast<Real>(r1), M_PI);

  // Test copy constructor
  Rad r2(r1);
  EXPECT_DOUBLE_EQ(static_cast<Real>(r2), M_PI);

  // Test move constructor
  Rad tmp(M_PI_2);
  Rad r3(std::move(tmp));
  EXPECT_DOUBLE_EQ(static_cast<Real>(r3), M_PI_2);
}

// Test Rad static factory methods inherited from Unit
TEST_F(RadTest, StaticFactoryMethods)
{
  // Test One() - should represent 1 radian
  auto one = Rad::One();
  EXPECT_DOUBLE_EQ(static_cast<Real>(one), 1.0);

  // Test Zero() - should represent 0 radians
  auto zero = Rad::Zero();
  EXPECT_DOUBLE_EQ(static_cast<Real>(zero), 0.0);
}

// Test Rad with common angular values
TEST_F(RadTest, CommonAngularValues)
{
  // Test pi radians (180 degrees)
  Rad pi(M_PI);
  EXPECT_DOUBLE_EQ(static_cast<Real>(pi), M_PI);

  // Test pi/2 radians (90 degrees)
  Rad half_pi(M_PI_2);
  EXPECT_DOUBLE_EQ(static_cast<Real>(half_pi), M_PI_2);

  // Test pi/4 radians (45 degrees)
  Rad quarter_pi(M_PI_4);
  EXPECT_DOUBLE_EQ(static_cast<Real>(quarter_pi), M_PI_4);

  // Test 2*pi radians (360 degrees)
  Rad two_pi(2.0 * M_PI);
  EXPECT_DOUBLE_EQ(static_cast<Real>(two_pi), 2.0 * M_PI);
}

// Test Rad arithmetic operations
TEST_F(RadTest, ArithmeticOperations)
{
  Rad r1(M_PI_2);  // 90 degrees
  Rad r2(M_PI_4);  // 45 degrees

  // Test addition (90° + 45° = 135°)
  auto sum = r1 + r2;
  EXPECT_DOUBLE_EQ(static_cast<Real>(sum), M_PI_2 + M_PI_4);

  // Test subtraction (90° - 45° = 45°)
  auto diff = r1 - r2;
  EXPECT_DOUBLE_EQ(static_cast<Real>(diff), M_PI_2 - M_PI_4);

  // Test multiplication (scaling angle)
  auto doubled = r2 * Rad(2.0);
  EXPECT_DOUBLE_EQ(static_cast<Real>(doubled), M_PI_2);

  // Test division (halving angle)
  auto halved = r1 / Rad(2.0);
  EXPECT_DOUBLE_EQ(static_cast<Real>(halved), M_PI_4);
}

// Test Rad comparison operations
TEST_F(RadTest, ComparisonOperations)
{
  Rad r1(M_PI_4);     // 45 degrees
  Rad r2(M_PI_2);     // 90 degrees
  Rad r3(M_PI_4);     // 45 degrees

  // Test equality
  EXPECT_TRUE(r1 == r3);
  EXPECT_FALSE(r1 == r2);

  // Test inequality
  EXPECT_FALSE(r1 != r3);
  EXPECT_TRUE(r1 != r2);

  // Test less than
  EXPECT_TRUE(r1 < r2);
  EXPECT_FALSE(r2 < r1);
  EXPECT_FALSE(r1 < r3);

  // Test greater than
  EXPECT_TRUE(r2 > r1);
  EXPECT_FALSE(r1 > r2);
  EXPECT_FALSE(r1 > r3);

  // Test less than or equal
  EXPECT_TRUE(r1 <= r2);
  EXPECT_TRUE(r1 <= r3);
  EXPECT_FALSE(r2 <= r1);

  // Test greater than or equal
  EXPECT_TRUE(r2 >= r1);
  EXPECT_TRUE(r1 >= r3);
  EXPECT_FALSE(r1 >= r2);
}

// Test Rad unary operations
TEST_F(RadTest, UnaryOperations)
{
  Rad r1(M_PI_4);

  // Test unary plus
  auto plus_r1 = +r1;
  EXPECT_DOUBLE_EQ(static_cast<Real>(plus_r1), M_PI_4);

  // Test unary minus (negative angle)
  auto minus_r1 = -r1;
  EXPECT_DOUBLE_EQ(static_cast<Real>(minus_r1), -M_PI_4);

  // Test double negative
  auto double_minus = -(-r1);
  EXPECT_DOUBLE_EQ(static_cast<Real>(double_minus), M_PI_4);
}

// Test Rad compound assignment operations
TEST_F(RadTest, CompoundAssignmentOperations)
{
  Rad r1(M_PI_4);  // 45 degrees
  Rad r2(M_PI_4);  // 45 degrees

  // Test +=
  r1 += r2;
  EXPECT_DOUBLE_EQ(static_cast<Real>(r1), M_PI_2);  // Should be 90 degrees

  // Test -=
  r1 -= r2;
  EXPECT_DOUBLE_EQ(static_cast<Real>(r1), M_PI_4);  // Back to 45 degrees

  // Test *=
  r1 *= Rad(2.0);
  EXPECT_DOUBLE_EQ(static_cast<Real>(r1), M_PI_2);  // Should be 90 degrees

  // Test /=
  r1 /= Rad(2.0);
  EXPECT_DOUBLE_EQ(static_cast<Real>(r1), M_PI_4);  // Back to 45 degrees
}

// Test Rad with trigonometric functions
TEST_F(RadTest, TrigonometricConsistency)
{
  // Test that Rad works correctly with trigonometric functions
  Rad r_0(0.0);
  Rad r_90(M_PI_2);
  Rad r_180(M_PI);
  Rad r_270(3.0 * M_PI_2);

  // Test sin values
  EXPECT_NEAR(std::sin(static_cast<Real>(r_0)), 0.0, 1e-10);
  EXPECT_NEAR(std::sin(static_cast<Real>(r_90)), 1.0, 1e-10);
  EXPECT_NEAR(std::sin(static_cast<Real>(r_180)), 0.0, 1e-10);
  EXPECT_NEAR(std::sin(static_cast<Real>(r_270)), -1.0, 1e-10);

  // Test cos values
  EXPECT_NEAR(std::cos(static_cast<Real>(r_0)), 1.0, 1e-10);
  EXPECT_NEAR(std::cos(static_cast<Real>(r_90)), 0.0, 1e-10);
  EXPECT_NEAR(std::cos(static_cast<Real>(r_180)), -1.0, 1e-10);
  EXPECT_NEAR(std::cos(static_cast<Real>(r_270)), 0.0, 1e-10);
}

// Test Rad assignment operations
TEST_F(RadTest, AssignmentOperations)
{
  Rad r1(M_PI_4);
  Rad r2(M_PI_2);

  // Test copy assignment
  r1 = r2;
  EXPECT_DOUBLE_EQ(static_cast<Real>(r1), M_PI_2);

  // Test move assignment
  r1 = Rad(M_PI);
  EXPECT_DOUBLE_EQ(static_cast<Real>(r1), M_PI);
}

// Test Rad edge cases
TEST_F(RadTest, EdgeCases)
{
  // Test with zero
  Rad zero = Rad::Zero();
  EXPECT_DOUBLE_EQ(static_cast<Real>(zero), 0.0);

  // Test with negative angles
  Rad negative(-M_PI_2);
  EXPECT_DOUBLE_EQ(static_cast<Real>(negative), -M_PI_2);

  // Test with very large angles (multiple rotations)
  Rad large(10.0 * M_PI);
  EXPECT_DOUBLE_EQ(static_cast<Real>(large), 10.0 * M_PI);

  // Test with very small angles
  Rad small(1e-10);
  EXPECT_DOUBLE_EQ(static_cast<Real>(small), 1e-10);
}

// Test Rad type properties
TEST_F(RadTest, TypeProperties)
{
  // Test that Rad inherits from Unit<Rad, Real>
  static_assert(std::is_base_of_v<Unit<Rad, Real>, Rad>);

  // Test that the underlying type is Real
  using Type = Rad::Type;
  static_assert(std::is_same_v<Type, Real>);

  // Test that Rad can be constructed and copied
  static_assert(std::is_constructible_v<Rad, Real>);
  static_assert(std::is_copy_constructible_v<Rad>);
  static_assert(std::is_move_constructible_v<Rad>);
  static_assert(std::is_copy_assignable_v<Rad>);
  static_assert(std::is_move_assignable_v<Rad>);
}
