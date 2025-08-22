/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Test/Random/RandomInteger.h"
#include "Rodin/Types.h"

using namespace Rodin;
using namespace Rodin::Test::Random;

class RandomIntegerTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// Test basic construction and type constraints
TEST_F(RandomIntegerTest, Construction)
{
  // Test default construction
  RandomInteger<int> ri1;
  EXPECT_NO_THROW(ri1());
  
  // Test construction with range
  RandomInteger<int> ri2(0, 100);
  EXPECT_NO_THROW(ri2());
  
  // Test construction with seed
  RandomInteger<int> ri3(0, 100, 42);
  EXPECT_NO_THROW(ri3());
  
  // Test that it works with different integer types
  RandomInteger<short> ri_short;
  RandomInteger<long> ri_long;
  RandomInteger<unsigned int> ri_uint;
  EXPECT_NO_THROW(ri_short());
  EXPECT_NO_THROW(ri_long());
  EXPECT_NO_THROW(ri_uint());
}

// Test range constraints
TEST_F(RandomIntegerTest, RangeConstraints)
{
  // Test that values are within specified range
  int min_val = 10;
  int max_val = 50;
  RandomInteger<int> ri(min_val, max_val, 12345);
  
  // Generate multiple values and check they're in range
  for (int i = 0; i < 100; ++i)
  {
    int value = ri();
    EXPECT_GE(value, min_val);
    EXPECT_LE(value, max_val);
  }
}

// Test seed functionality
TEST_F(RandomIntegerTest, SeedFunctionality)
{
  // Test getSeed
  unsigned int seed = 42;
  RandomInteger<int> ri(0, 100, seed);
  EXPECT_EQ(ri.getSeed(), seed);
  
  // Test setSeed and reproducibility
  RandomInteger<int> ri1(0, 100, 123);
  RandomInteger<int> ri2(0, 100, 456);
  
  // Set the same seed for both
  ri1.setSeed(789);
  ri2.setSeed(789);
  
  // They should generate the same sequence
  for (int i = 0; i < 10; ++i)
  {
    EXPECT_EQ(ri1(), ri2());
  }
}

// Test deterministic behavior with fixed seed
TEST_F(RandomIntegerTest, DeterministicBehavior)
{
  unsigned int seed = 12345;
  
  // Create two generators with the same seed
  RandomInteger<int> ri1(0, 1000, seed);
  RandomInteger<int> ri2(0, 1000, seed);
  
  // They should generate identical sequences
  std::vector<int> sequence1, sequence2;
  for (int i = 0; i < 50; ++i)
  {
    sequence1.push_back(ri1());
    sequence2.push_back(ri2());
  }
  
  EXPECT_EQ(sequence1, sequence2);
}

// Test different range scenarios
TEST_F(RandomIntegerTest, DifferentRanges)
{
  // Test positive range
  RandomInteger<int> ri_pos(1, 100, 111);
  for (int i = 0; i < 20; ++i)
  {
    int val = ri_pos();
    EXPECT_GE(val, 1);
    EXPECT_LE(val, 100);
  }
  
  // Test negative range
  RandomInteger<int> ri_neg(-100, -1, 222);
  for (int i = 0; i < 20; ++i)
  {
    int val = ri_neg();
    EXPECT_GE(val, -100);
    EXPECT_LE(val, -1);
  }
  
  // Test range spanning zero
  RandomInteger<int> ri_span(-50, 50, 333);
  for (int i = 0; i < 20; ++i)
  {
    int val = ri_span();
    EXPECT_GE(val, -50);
    EXPECT_LE(val, 50);
  }
  
  // Test single value range
  RandomInteger<int> ri_single(42, 42, 444);
  for (int i = 0; i < 20; ++i)
  {
    int val = ri_single();
    EXPECT_EQ(val, 42);
  }
}

// Test with different integer types
TEST_F(RandomIntegerTest, DifferentTypes)
{
  // Test with short
  RandomInteger<short> ri_short(0, 100, 555);
  for (int i = 0; i < 10; ++i)
  {
    short val = ri_short();
    EXPECT_GE(val, 0);
    EXPECT_LE(val, 100);
  }
  
  // Test with long
  RandomInteger<long> ri_long(0, 1000000L, 666);
  for (int i = 0; i < 10; ++i)
  {
    long val = ri_long();
    EXPECT_GE(val, 0L);
    EXPECT_LE(val, 1000000L);
  }
  
  // Test with unsigned int
  RandomInteger<unsigned int> ri_uint(100u, 200u, 777);
  for (int i = 0; i < 10; ++i)
  {
    unsigned int val = ri_uint();
    EXPECT_GE(val, 100u);
    EXPECT_LE(val, 200u);
  }
  
  // Test with char
  RandomInteger<char> ri_char('A', 'Z', 888);
  for (int i = 0; i < 10; ++i)
  {
    char val = ri_char();
    EXPECT_GE(val, 'A');
    EXPECT_LE(val, 'Z');
  }
}

// Test statistical properties (basic distribution test)
TEST_F(RandomIntegerTest, StatisticalProperties)
{
  RandomInteger<int> ri(0, 9, 999);  // 10 possible values (0-9)
  
  std::map<int, int> counts;
  const int numSamples = 10000;
  
  // Generate samples
  for (int i = 0; i < numSamples; ++i)
  {
    int val = ri();
    counts[val]++;
  }
  
  // Check that all values in range appeared
  for (int i = 0; i <= 9; ++i)
  {
    EXPECT_GT(counts[i], 0) << "Value " << i << " never appeared";
  }
  
  // Check roughly uniform distribution (each value should appear ~1000 times)
  for (int i = 0; i <= 9; ++i)
  {
    EXPECT_NEAR(counts[i], 1000, 200) << "Value " << i << " appeared " << counts[i] << " times";
  }
}

// Test edge cases
TEST_F(RandomIntegerTest, EdgeCases)
{
  // Test with very large range
  RandomInteger<long> ri_large(-1000000L, 1000000L, 101);
  for (int i = 0; i < 10; ++i)
  {
    long val = ri_large();
    EXPECT_GE(val, -1000000L);
    EXPECT_LE(val, 1000000L);
  }
  
  // Test with minimum possible values for int
  RandomInteger<int> ri_min(std::numeric_limits<int>::min(), 
                           std::numeric_limits<int>::min() + 10, 102);
  for (int i = 0; i < 10; ++i)
  {
    int val = ri_min();
    EXPECT_GE(val, std::numeric_limits<int>::min());
    EXPECT_LE(val, std::numeric_limits<int>::min() + 10);
  }
  
  // Test with maximum possible values for int
  RandomInteger<int> ri_max(std::numeric_limits<int>::max() - 10,
                           std::numeric_limits<int>::max(), 103);
  for (int i = 0; i < 10; ++i)
  {
    int val = ri_max();
    EXPECT_GE(val, std::numeric_limits<int>::max() - 10);
    EXPECT_LE(val, std::numeric_limits<int>::max());
  }
}

// Test operator() multiple calls
TEST_F(RandomIntegerTest, OperatorCall)
{
  RandomInteger<int> ri(-50, 50, 104);
  
  // Test that multiple calls work
  std::vector<int> values;
  for (int i = 0; i < 100; ++i)
  {
    values.push_back(ri());
  }
  
  // Check all values are valid
  for (int val : values)
  {
    EXPECT_GE(val, -50);
    EXPECT_LE(val, 50);
  }
  
  // Check that we get different values (with high probability)
  std::set<int> uniqueValues(values.begin(), values.end());
  EXPECT_GT(uniqueValues.size(), 50);  // Should have many unique values
}

// Test copy and assignment behavior
TEST_F(RandomIntegerTest, CopyAndAssignment)
{
  RandomInteger<int> ri1(0, 100, 105);
  
  // Test copy constructor
  RandomInteger<int> ri2 = ri1;
  EXPECT_EQ(ri1.getSeed(), ri2.getSeed());
  
  // Test assignment
  RandomInteger<int> ri3(200, 300, 999);
  ri3 = ri1;
  EXPECT_EQ(ri1.getSeed(), ri3.getSeed());
}

// Test boundary values
TEST_F(RandomIntegerTest, BoundaryValues)
{
  // Test that boundary values can be generated
  RandomInteger<int> ri(0, 2, 106);  // Only 3 possible values: 0, 1, 2
  
  std::set<int> generatedValues;
  for (int i = 0; i < 1000; ++i)
  {
    generatedValues.insert(ri());
  }
  
  // Should have generated all possible values
  EXPECT_TRUE(generatedValues.count(0) > 0);
  EXPECT_TRUE(generatedValues.count(1) > 0);
  EXPECT_TRUE(generatedValues.count(2) > 0);
  EXPECT_EQ(generatedValues.size(), 3);
}

// Test with bool-like range
TEST_F(RandomIntegerTest, BooleanLikeRange)
{
  RandomInteger<int> ri(0, 1, 107);  // Binary values
  
  std::map<int, int> counts;
  for (int i = 0; i < 1000; ++i)
  {
    int val = ri();
    EXPECT_TRUE(val == 0 || val == 1);
    counts[val]++;
  }
  
  // Both values should appear with roughly equal frequency
  EXPECT_GT(counts[0], 300);  // At least 30% 
  EXPECT_GT(counts[1], 300);  // At least 30%
  EXPECT_EQ(counts[0] + counts[1], 1000);
}

// Test unsigned types specifically
TEST_F(RandomIntegerTest, UnsignedTypes)
{
  // Test unsigned char
  RandomInteger<unsigned char> ri_uchar(100, 200, 108);
  for (int i = 0; i < 20; ++i)
  {
    unsigned char val = ri_uchar();
    EXPECT_GE(val, 100);
    EXPECT_LE(val, 200);
  }
  
  // Test unsigned long
  RandomInteger<unsigned long> ri_ulong(1000000UL, 2000000UL, 109);
  for (int i = 0; i < 20; ++i)
  {
    unsigned long val = ri_ulong();
    EXPECT_GE(val, 1000000UL);
    EXPECT_LE(val, 2000000UL);
  }
}

// Test type safety (compile-time test)
TEST_F(RandomIntegerTest, TypeSafety)
{
  // These should compile fine
  RandomInteger<int> ri_int;
  RandomInteger<short> ri_short;
  RandomInteger<long> ri_long;
  RandomInteger<unsigned int> ri_uint;
  RandomInteger<char> ri_char;
  
  // Test that type deduction works correctly
  auto val_int = ri_int();
  auto val_short = ri_short();
  auto val_long = ri_long();
  auto val_uint = ri_uint();
  auto val_char = ri_char();
  
  static_assert(std::is_same_v<decltype(val_int), int>);
  static_assert(std::is_same_v<decltype(val_short), short>);
  static_assert(std::is_same_v<decltype(val_long), long>);
  static_assert(std::is_same_v<decltype(val_uint), unsigned int>);
  static_assert(std::is_same_v<decltype(val_char), char>);
}

// Test large range efficiency
TEST_F(RandomIntegerTest, LargeRange)
{
  // Test with a very large range to ensure no overflow issues
  RandomInteger<long long> ri(std::numeric_limits<long long>::min() / 2,
                             std::numeric_limits<long long>::max() / 2, 110);
  
  for (int i = 0; i < 10; ++i)
  {
    long long val = ri();
    EXPECT_GE(val, std::numeric_limits<long long>::min() / 2);
    EXPECT_LE(val, std::numeric_limits<long long>::max() / 2);
  }
}