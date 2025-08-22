/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Math/Constants.h"
#include "Rodin/Types.h"

using namespace Rodin;
using namespace Rodin::Math::Constants;

class ConstantsTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// Test pi() function
TEST_F(ConstantsTest, PiFunction)
{
  // Test that pi() returns the correct value
  EXPECT_DOUBLE_EQ(pi(), M_PI);
  
  // Test basic mathematical properties of pi
  EXPECT_GT(pi(), 3.14);
  EXPECT_LT(pi(), 3.15);
  
  // Test that pi is greater than 3 and less than 4
  EXPECT_GT(pi(), 3.0);
  EXPECT_LT(pi(), 4.0);
  
  // Test constexpr nature
  constexpr Real pi_val = pi();
  EXPECT_DOUBLE_EQ(pi_val, M_PI);
}

// Test epsilon() function
TEST_F(ConstantsTest, EpsilonFunction)
{
  // Test that epsilon returns machine epsilon for Real type
  EXPECT_DOUBLE_EQ(epsilon(), std::numeric_limits<Real>::epsilon());
  
  // Test that epsilon is positive and very small
  EXPECT_GT(epsilon(), 0.0);
  EXPECT_LT(epsilon(), 1e-10);
  
  // Test constexpr nature
  constexpr Real eps_val = epsilon();
  EXPECT_DOUBLE_EQ(eps_val, std::numeric_limits<Real>::epsilon());
  
  // Test mathematical property: 1.0 + epsilon() > 1.0
  EXPECT_GT(1.0 + epsilon(), 1.0);
}

// Test zero() function
TEST_F(ConstantsTest, ZeroFunction)
{
  // Test that zero() returns exactly 0.0
  EXPECT_DOUBLE_EQ(zero(), 0.0);
  EXPECT_DOUBLE_EQ(zero(), Real(0));
  
  // Test constexpr nature
  constexpr Real zero_val = zero();
  EXPECT_DOUBLE_EQ(zero_val, 0.0);
  
  // Test mathematical properties
  EXPECT_EQ(zero() + 1.0, 1.0);
  EXPECT_EQ(zero() * 5.0, 0.0);
  EXPECT_EQ(zero() - zero(), 0.0);
}

// Test one() function
TEST_F(ConstantsTest, OneFunction)
{
  // Test that one() returns exactly 1.0
  EXPECT_DOUBLE_EQ(one(), 1.0);
  EXPECT_DOUBLE_EQ(one(), Real(1));
  
  // Test constexpr nature
  constexpr Real one_val = one();
  EXPECT_DOUBLE_EQ(one_val, 1.0);
  
  // Test mathematical properties
  EXPECT_EQ(one() + zero(), 1.0);
  EXPECT_EQ(one() * 5.0, 5.0);
  EXPECT_EQ(one() - one(), 0.0);
  EXPECT_EQ(one() / one(), 1.0);
}

// Test isZero() function
TEST_F(ConstantsTest, IsZeroFunction)
{
  // Test with actual zero
  EXPECT_TRUE(isZero(0.0));
  EXPECT_TRUE(isZero(zero()));
  EXPECT_TRUE(isZero(Real(0)));
  
  // Test with non-zero values
  EXPECT_FALSE(isZero(1.0));
  EXPECT_FALSE(isZero(-1.0));
  EXPECT_FALSE(isZero(epsilon()));
  EXPECT_FALSE(isZero(1e-100));
  EXPECT_FALSE(isZero(-1e-100));
  
  // Test constexpr nature
  constexpr Boolean is_zero_test = isZero(0.0);
  EXPECT_TRUE(is_zero_test);
  
  constexpr Boolean is_not_zero_test = isZero(1.0);
  EXPECT_FALSE(is_not_zero_test);
}

// Test isOne() function
TEST_F(ConstantsTest, IsOneFunction)
{
  // Test with actual one
  EXPECT_TRUE(isOne(1.0));
  EXPECT_TRUE(isOne(one()));
  EXPECT_TRUE(isOne(Real(1)));
  
  // Test with non-one values
  EXPECT_FALSE(isOne(0.0));
  EXPECT_FALSE(isOne(2.0));
  EXPECT_FALSE(isOne(-1.0));
  EXPECT_FALSE(isOne(1.0 + epsilon()));
  EXPECT_FALSE(isOne(1.0 - epsilon()));
  
  // Test constexpr nature
  constexpr Boolean is_one_test = isOne(1.0);
  EXPECT_TRUE(is_one_test);
  
  constexpr Boolean is_not_one_test = isOne(2.0);
  EXPECT_FALSE(is_not_one_test);
}

// Test mathematical relationships between constants
TEST_F(ConstantsTest, MathematicalRelationships)
{
  // Test zero and one relationship
  EXPECT_EQ(one() - one(), zero());
  EXPECT_EQ(zero() + one(), one());
  EXPECT_EQ(one() * zero(), zero());
  
  // Test pi relationships
  EXPECT_GT(pi(), one());
  EXPECT_GT(pi(), 2.0 * one());
  EXPECT_LT(pi(), 4.0 * one());
  
  // Test epsilon relationships
  EXPECT_GT(epsilon(), zero());
  EXPECT_LT(epsilon(), one());
  EXPECT_LT(epsilon(), pi());
  
  // Test that 1 + epsilon is not exactly 1 (but very close)
  EXPECT_NE(one() + epsilon(), one());
  EXPECT_NEAR(one() + epsilon(), one(), epsilon() * 10);
}

// Test usage in practical computations
TEST_F(ConstantsTest, PracticalUsage)
{
  // Test circle area computation: A = π * r²
  Real radius = 2.0;
  Real area = pi() * radius * radius;
  EXPECT_NEAR(area, 4.0 * M_PI, 1e-10);
  
  // Test circle circumference: C = 2 * π * r
  Real circumference = 2.0 * pi() * radius;
  EXPECT_NEAR(circumference, 4.0 * M_PI, 1e-10);
  
  // Test that we can use zero() and one() in computations
  Real computation = (one() + one()) * pi() / (2.0 * one());
  EXPECT_DOUBLE_EQ(computation, pi());
  
  // Test epsilon in numerical comparisons
  Real a = 1.0;
  Real b = 1.0 + epsilon() / 2.0;
  EXPECT_TRUE(abs(a - b) < epsilon());
}

// Test type consistency
TEST_F(ConstantsTest, TypeConsistency)
{
  // Test that all functions return Real type
  static_assert(std::is_same_v<decltype(pi()), Real>);
  static_assert(std::is_same_v<decltype(epsilon()), Real>);
  static_assert(std::is_same_v<decltype(zero()), Real>);
  static_assert(std::is_same_v<decltype(one()), Real>);
  
  // Test that Boolean functions return Boolean type
  static_assert(std::is_same_v<decltype(isZero(1.0)), Boolean>);
  static_assert(std::is_same_v<decltype(isOne(1.0)), Boolean>);
}

// Test constexpr evaluation
TEST_F(ConstantsTest, ConstexprEvaluation)
{
  // Test that all constant functions can be evaluated at compile time
  constexpr Real pi_val = pi();
  constexpr Real eps_val = epsilon();
  constexpr Real zero_val = zero();
  constexpr Real one_val = one();
  
  // Test Boolean functions at compile time
  constexpr Boolean zero_check = isZero(0.0);
  constexpr Boolean one_check = isOne(1.0);
  constexpr Boolean not_zero_check = isZero(1.0);
  constexpr Boolean not_one_check = isOne(0.0);
  
  // Verify values
  EXPECT_DOUBLE_EQ(pi_val, M_PI);
  EXPECT_DOUBLE_EQ(eps_val, std::numeric_limits<Real>::epsilon());
  EXPECT_DOUBLE_EQ(zero_val, 0.0);
  EXPECT_DOUBLE_EQ(one_val, 1.0);
  
  EXPECT_TRUE(zero_check);
  EXPECT_TRUE(one_check);
  EXPECT_FALSE(not_zero_check);
  EXPECT_FALSE(not_one_check);
}

// Test precision and numerical stability
TEST_F(ConstantsTest, PrecisionAndStability)
{
  // Test that pi() is stable across multiple calls
  Real pi1 = pi();
  Real pi2 = pi();
  EXPECT_EQ(pi1, pi2);
  
  // Test that epsilon is the smallest representable difference
  Real one_plus_eps = 1.0 + epsilon();
  Real one_plus_half_eps = 1.0 + epsilon() / 2.0;
  EXPECT_GT(one_plus_eps, 1.0);
  EXPECT_EQ(one_plus_half_eps, 1.0);  // Should round to 1.0
  
  // Test numerical limits consistency
  EXPECT_EQ(epsilon(), std::numeric_limits<Real>::epsilon());
  EXPECT_LE(zero(), std::numeric_limits<Real>::min());
  EXPECT_LE(one(), std::numeric_limits<Real>::max());
}