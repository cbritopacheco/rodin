/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Math/Unit.h"
#include "Rodin/Types.h"

using namespace Rodin;
using namespace Rodin::Math;

// Test unit class for testing Unit template
class TestUnit : public Unit<TestUnit, Real>
{
  public:
    using Parent = Unit<TestUnit, Real>;
    using Parent::Parent;
};

class UnitTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// Test Unit construction
TEST_F(UnitTest, Construction)
{
  // Test value constructor
  TestUnit u1(5.0);
  EXPECT_DOUBLE_EQ(static_cast<Real>(u1), 5.0);

  // Test copy constructor
  TestUnit u2(u1);
  EXPECT_DOUBLE_EQ(static_cast<Real>(u2), 5.0);

  // Test move constructor
  TestUnit tmp(3.14);
  TestUnit u3(std::move(tmp));
  EXPECT_DOUBLE_EQ(static_cast<Real>(u3), 3.14);
}

// Test Unit static factory methods
TEST_F(UnitTest, StaticFactoryMethods)
{
  // Test One()
  auto one = TestUnit::One();
  EXPECT_DOUBLE_EQ(static_cast<Real>(one), 1.0);

  // Test Zero()
  auto zero = TestUnit::Zero();
  EXPECT_DOUBLE_EQ(static_cast<Real>(zero), 0.0);
}

// Test Unit assignment operators
TEST_F(UnitTest, AssignmentOperators)
{
  TestUnit u1(5.0);
  TestUnit u2(10.0);

  // Test copy assignment
  u1 = u2;
  EXPECT_DOUBLE_EQ(static_cast<Real>(u1), 10.0);

  // Test move assignment
  TestUnit tmp(7.5);
  u1 = std::move(tmp);
  EXPECT_DOUBLE_EQ(static_cast<Real>(u1), 7.5);
}

// Test Unit comparison operators
TEST_F(UnitTest, ComparisonOperators)
{
  TestUnit u1(5.0);
  TestUnit u2(5.0);
  TestUnit u3(10.0);
  TestUnit u4(2.0);

  // Test equality
  EXPECT_TRUE(u1 == u2);
  EXPECT_FALSE(u1 == u3);

  // Test inequality
  EXPECT_FALSE(u1 != u2);
  EXPECT_TRUE(u1 != u3);

  // Test less than
  EXPECT_TRUE(u4 < u1);
  EXPECT_FALSE(u1 < u4);
  EXPECT_FALSE(u1 < u2);

  // Test greater than
  EXPECT_TRUE(u1 > u4);
  EXPECT_FALSE(u4 > u1);
  EXPECT_FALSE(u1 > u2);

  // Test less than or equal
  EXPECT_TRUE(u4 <= u1);
  EXPECT_TRUE(u1 <= u2);
  EXPECT_FALSE(u1 <= u4);

  // Test greater than or equal
  EXPECT_TRUE(u1 >= u4);
  EXPECT_TRUE(u1 >= u2);
  EXPECT_FALSE(u4 >= u1);
}

// Test Unit arithmetic operators
TEST_F(UnitTest, ArithmeticOperators)
{
  TestUnit u1(5.0);
  TestUnit u2(3.0);

  // Test addition
  auto sum = u1 + u2;
  EXPECT_DOUBLE_EQ(static_cast<Real>(sum), 8.0);

  // Test subtraction
  auto diff = u1 - u2;
  EXPECT_DOUBLE_EQ(static_cast<Real>(diff), 2.0);

  // Test multiplication
  auto prod = u1 * u2;
  EXPECT_DOUBLE_EQ(static_cast<Real>(prod), 15.0);

  // Test division
  auto quot = u1 / u2;
  EXPECT_DOUBLE_EQ(static_cast<Real>(quot), 5.0 / 3.0);
}

// Test Unit unary operators
TEST_F(UnitTest, UnaryOperators)
{
  TestUnit u1(5.0);
  TestUnit u2(-3.0);

  // Test unary plus
  auto plus_u1 = +u1;
  EXPECT_DOUBLE_EQ(static_cast<Real>(plus_u1), 5.0);

  // Test unary minus
  auto minus_u1 = -u1;
  EXPECT_DOUBLE_EQ(static_cast<Real>(minus_u1), -5.0);

  auto minus_u2 = -u2;
  EXPECT_DOUBLE_EQ(static_cast<Real>(minus_u2), 3.0);
}

// Test Unit compound assignment operators
TEST_F(UnitTest, CompoundAssignmentOperators)
{
  TestUnit u1(5.0);
  TestUnit u2(3.0);

  // Test +=
  u1 += u2;
  EXPECT_DOUBLE_EQ(static_cast<Real>(u1), 8.0);

  // Test -=
  u1 -= u2;
  EXPECT_DOUBLE_EQ(static_cast<Real>(u1), 5.0);

  // Test *=
  u1 *= u2;
  EXPECT_DOUBLE_EQ(static_cast<Real>(u1), 15.0);

  // Test /=
  u1 /= u2;
  EXPECT_DOUBLE_EQ(static_cast<Real>(u1), 5.0);
}

// Test Unit explicit type conversion
TEST_F(UnitTest, TypeConversion)
{
  TestUnit u1(3.14159);

  // Test explicit conversion to underlying type
  Real value = static_cast<Real>(u1);
  EXPECT_DOUBLE_EQ(value, 3.14159);
}

// Test Unit with different underlying types
TEST_F(UnitTest, DifferentTypes)
{
  // Test with integer type
  class IntUnit : public Unit<IntUnit, int> 
  {
    public:
      using Parent = Unit<IntUnit, int>;
      using Parent::Parent;
  };

  IntUnit iu1(42);
  IntUnit iu2(10);

  EXPECT_EQ(static_cast<int>(iu1), 42);

  auto sum = iu1 + iu2;
  EXPECT_EQ(static_cast<int>(sum), 52);

  // Test with float type
  class FloatUnit : public Unit<FloatUnit, float>
  {
    public:
      using Parent = Unit<FloatUnit, float>;
      using Parent::Parent;
  };

  FloatUnit fu1(2.5f);
  FloatUnit fu2(1.5f);

  EXPECT_FLOAT_EQ(static_cast<float>(fu1), 2.5f);

  auto prod = fu1 * fu2;
  EXPECT_FLOAT_EQ(static_cast<float>(prod), 3.75f);
}

// Test Unit edge cases
TEST_F(UnitTest, EdgeCases)
{
  // Test with zero values
  TestUnit zero1 = TestUnit::Zero();
  TestUnit zero2(0.0);
  EXPECT_TRUE(zero1 == zero2);

  // Test division by zero behavior (undefined, but shouldn't crash)
  TestUnit u1(5.0);
  TestUnit zero(0.0);
  // Note: Division by zero is undefined behavior, but the operation should compile
  auto result = u1 / zero;
  EXPECT_TRUE(std::isinf(static_cast<Real>(result)) || std::isnan(static_cast<Real>(result)));

  // Test with very small values
  TestUnit small1(1e-10);
  TestUnit small2(1e-10);
  EXPECT_TRUE(small1 == small2);

  // Test with very large values
  TestUnit large1(1e10);
  TestUnit large2(1e10);
  EXPECT_TRUE(large1 == large2);
}
