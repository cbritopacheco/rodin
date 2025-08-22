/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Test/Random/RandomFloat.h"
#include "Rodin/Types.h"

using namespace Rodin;
using namespace Rodin::Test::Random;

class RandomFloatTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// Test basic construction and type constraints
TEST_F(RandomFloatTest, Construction)
{
  // Test default construction
  RandomFloat<Real> rf1;
  EXPECT_NO_THROW(rf1());
  
  // Test construction with range
  RandomFloat<Real> rf2(0.0, 1.0);
  EXPECT_NO_THROW(rf2());
  
  // Test construction with seed
  RandomFloat<Real> rf3(0.0, 1.0, 42);
  EXPECT_NO_THROW(rf3());
  
  // Test that it works with different floating point types
  RandomFloat<float> rf_float;
  RandomFloat<double> rf_double;
  EXPECT_NO_THROW(rf_float());
  EXPECT_NO_THROW(rf_double());
}

// Test range constraints
TEST_F(RandomFloatTest, RangeConstraints)
{
  // Test that values are within specified range
  Real min_val = 2.0;
  Real max_val = 5.0;
  RandomFloat<Real> rf(min_val, max_val, 12345);
  
  // Generate multiple values and check they're in range
  for (int i = 0; i < 100; ++i)
  {
    Real value = rf();
    EXPECT_GE(value, min_val);
    EXPECT_LE(value, max_val);
  }
}

// Test seed functionality
TEST_F(RandomFloatTest, SeedFunctionality)
{
  // Test getSeed
  unsigned int seed = 42;
  RandomFloat<Real> rf(0.0, 1.0, seed);
  EXPECT_EQ(rf.getSeed(), seed);
  
  // Test setSeed and reproducibility
  RandomFloat<Real> rf1(0.0, 1.0, 123);
  RandomFloat<Real> rf2(0.0, 1.0, 456);
  
  // Set the same seed for both
  rf1.setSeed(789);
  rf2.setSeed(789);
  
  // They should generate the same sequence
  for (int i = 0; i < 10; ++i)
  {
    EXPECT_EQ(rf1(), rf2());
  }
}

// Test deterministic behavior with fixed seed
TEST_F(RandomFloatTest, DeterministicBehavior)
{
  unsigned int seed = 12345;
  
  // Create two generators with the same seed
  RandomFloat<Real> rf1(0.0, 1.0, seed);
  RandomFloat<Real> rf2(0.0, 1.0, seed);
  
  // They should generate identical sequences
  std::vector<Real> sequence1, sequence2;
  for (int i = 0; i < 50; ++i)
  {
    sequence1.push_back(rf1());
    sequence2.push_back(rf2());
  }
  
  EXPECT_EQ(sequence1, sequence2);
}

// Test different range scenarios
TEST_F(RandomFloatTest, DifferentRanges)
{
  // Test positive range
  RandomFloat<Real> rf_pos(1.0, 10.0, 111);
  for (int i = 0; i < 20; ++i)
  {
    Real val = rf_pos();
    EXPECT_GE(val, 1.0);
    EXPECT_LE(val, 10.0);
  }
  
  // Test negative range
  RandomFloat<Real> rf_neg(-10.0, -1.0, 222);
  for (int i = 0; i < 20; ++i)
  {
    Real val = rf_neg();
    EXPECT_GE(val, -10.0);
    EXPECT_LE(val, -1.0);
  }
  
  // Test range spanning zero
  RandomFloat<Real> rf_span(-5.0, 5.0, 333);
  for (int i = 0; i < 20; ++i)
  {
    Real val = rf_span();
    EXPECT_GE(val, -5.0);
    EXPECT_LE(val, 5.0);
  }
  
  // Test very small range
  RandomFloat<Real> rf_small(1.0, 1.1, 444);
  for (int i = 0; i < 20; ++i)
  {
    Real val = rf_small();
    EXPECT_GE(val, 1.0);
    EXPECT_LE(val, 1.1);
  }
}

// Test with different floating point types
TEST_F(RandomFloatTest, DifferentTypes)
{
  // Test with float
  RandomFloat<float> rf_float(0.0f, 1.0f, 555);
  for (int i = 0; i < 10; ++i)
  {
    float val = rf_float();
    EXPECT_GE(val, 0.0f);
    EXPECT_LE(val, 1.0f);
  }
  
  // Test with double
  RandomFloat<double> rf_double(0.0, 1.0, 666);
  for (int i = 0; i < 10; ++i)
  {
    double val = rf_double();
    EXPECT_GE(val, 0.0);
    EXPECT_LE(val, 1.0);
  }
  
  // Test with long double
  RandomFloat<long double> rf_ldouble(0.0L, 1.0L, 777);
  for (int i = 0; i < 10; ++i)
  {
    long double val = rf_ldouble();
    EXPECT_GE(val, 0.0L);
    EXPECT_LE(val, 1.0L);
  }
}

// Test statistical properties (basic distribution test)
TEST_F(RandomFloatTest, StatisticalProperties)
{
  RandomFloat<Real> rf(0.0, 1.0, 888);
  
  std::vector<Real> samples;
  const int numSamples = 1000;
  
  // Generate samples
  for (int i = 0; i < numSamples; ++i)
  {
    samples.push_back(rf());
  }
  
  // Check that we have values across the range
  Real minSample = *std::min_element(samples.begin(), samples.end());
  Real maxSample = *std::max_element(samples.begin(), samples.end());
  
  EXPECT_GT(minSample, 0.0);   // Should be greater than lower bound
  EXPECT_LT(maxSample, 1.0);   // Should be less than upper bound
  EXPECT_GT(maxSample - minSample, 0.5);  // Should span significant portion of range
  
  // Basic mean test (should be around 0.5 for uniform distribution)
  Real sum = std::accumulate(samples.begin(), samples.end(), 0.0);
  Real mean = sum / numSamples;
  EXPECT_NEAR(mean, 0.5, 0.1);  // Within 10% of expected mean
}

// Test edge cases
TEST_F(RandomFloatTest, EdgeCases)
{
  // Test with very large range
  RandomFloat<Real> rf_large(-1e6, 1e6, 999);
  for (int i = 0; i < 10; ++i)
  {
    Real val = rf_large();
    EXPECT_GE(val, -1e6);
    EXPECT_LE(val, 1e6);
  }
  
  // Test with very small range
  RandomFloat<Real> rf_tiny(1.0, 1.0 + 1e-10, 101112);
  for (int i = 0; i < 10; ++i)
  {
    Real val = rf_tiny();
    EXPECT_GE(val, 1.0);
    EXPECT_LE(val, 1.0 + 1e-10);
  }
  
  // Test with equal bounds (should always return the same value)
  RandomFloat<Real> rf_equal(5.0, 5.0, 131415);
  for (int i = 0; i < 10; ++i)
  {
    Real val = rf_equal();
    EXPECT_DOUBLE_EQ(val, 5.0);
  }
}

// Test operator() multiple calls
TEST_F(RandomFloatTest, OperatorCall)
{
  RandomFloat<Real> rf(-2.0, 3.0, 161718);
  
  // Test that multiple calls work
  std::vector<Real> values;
  for (int i = 0; i < 100; ++i)
  {
    values.push_back(rf());
  }
  
  // Check all values are valid
  for (Real val : values)
  {
    EXPECT_GE(val, -2.0);
    EXPECT_LE(val, 3.0);
  }
  
  // Check that we get different values (with high probability)
  std::set<Real> uniqueValues(values.begin(), values.end());
  EXPECT_GT(uniqueValues.size(), 50);  // Should have many unique values
}

// Test copy and assignment behavior
TEST_F(RandomFloatTest, CopyAndAssignment)
{
  RandomFloat<Real> rf1(0.0, 1.0, 192021);
  
  // Test copy constructor
  RandomFloat<Real> rf2 = rf1;
  EXPECT_EQ(rf1.getSeed(), rf2.getSeed());
  
  // Test assignment
  RandomFloat<Real> rf3(2.0, 3.0, 999);
  rf3 = rf1;
  EXPECT_EQ(rf1.getSeed(), rf3.getSeed());
}

// Test type safety (compile-time test)
TEST_F(RandomFloatTest, TypeSafety)
{
  // These should compile fine
  RandomFloat<float> rf_float;
  RandomFloat<double> rf_double;
  RandomFloat<long double> rf_ldouble;
  
  // Test that type deduction works correctly
  auto val_float = rf_float();
  auto val_double = rf_double();
  auto val_ldouble = rf_ldouble();
  
  static_assert(std::is_same_v<decltype(val_float), float>);
  static_assert(std::is_same_v<decltype(val_double), double>);
  static_assert(std::is_same_v<decltype(val_ldouble), long double>);
}