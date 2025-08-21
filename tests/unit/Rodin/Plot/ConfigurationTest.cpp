/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Plot/Configuration.h"

using namespace Rodin;
using namespace Rodin::Plot;

class PlotConfigurationTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// Test basic construction
TEST_F(PlotConfigurationTest, Construction)
{
  // Test default construction
  Configuration config;
  EXPECT_NO_THROW(config.getSampleCount());
  
  // Default sample count should be reasonable (> 0)
  EXPECT_GT(config.getSampleCount(), 0);
}

// Test sample count getter and setter
TEST_F(PlotConfigurationTest, SampleCount)
{
  Configuration config;
  
  // Test setting sample count
  int testSampleCount = 16;
  Configuration& result = config.setSampleCount(testSampleCount);
  
  // Check method chaining (returns reference to self)
  EXPECT_EQ(&result, &config);
  
  // Check that value was set correctly
  EXPECT_EQ(config.getSampleCount(), testSampleCount);
}

// Test method chaining
TEST_F(PlotConfigurationTest, MethodChaining)
{
  Configuration config;
  
  // Test that multiple calls can be chained
  Configuration& result = config.setSampleCount(8);
  
  // Verify it returns a reference to the same object
  EXPECT_EQ(&result, &config);
  
  // Test that we can chain another call
  Configuration& result2 = result.setSampleCount(32);
  EXPECT_EQ(&result2, &config);
  EXPECT_EQ(config.getSampleCount(), 32);
}

// Test various sample count values
TEST_F(PlotConfigurationTest, VariousSampleCounts)
{
  Configuration config;
  
  // Test small values
  config.setSampleCount(1);
  EXPECT_EQ(config.getSampleCount(), 1);
  
  config.setSampleCount(2);
  EXPECT_EQ(config.getSampleCount(), 2);
  
  config.setSampleCount(4);
  EXPECT_EQ(config.getSampleCount(), 4);
  
  // Test typical anti-aliasing values
  config.setSampleCount(8);
  EXPECT_EQ(config.getSampleCount(), 8);
  
  config.setSampleCount(16);
  EXPECT_EQ(config.getSampleCount(), 16);
  
  config.setSampleCount(32);
  EXPECT_EQ(config.getSampleCount(), 32);
  
  // Test large values
  config.setSampleCount(64);
  EXPECT_EQ(config.getSampleCount(), 64);
  
  config.setSampleCount(128);
  EXPECT_EQ(config.getSampleCount(), 128);
}

// Test edge cases
TEST_F(PlotConfigurationTest, EdgeCases)
{
  Configuration config;
  
  // Test zero (might be valid depending on implementation)
  config.setSampleCount(0);
  EXPECT_EQ(config.getSampleCount(), 0);
  
  // Test negative values (behavior depends on implementation)
  config.setSampleCount(-1);
  EXPECT_EQ(config.getSampleCount(), -1);
  
  // Test very large values
  config.setSampleCount(1000000);
  EXPECT_EQ(config.getSampleCount(), 1000000);
}

// Test copy semantics
TEST_F(PlotConfigurationTest, CopySemantics)
{
  Configuration config1;
  config1.setSampleCount(42);
  
  // Test copy constructor
  Configuration config2 = config1;
  EXPECT_EQ(config2.getSampleCount(), 42);
  
  // Verify they are independent
  config1.setSampleCount(24);
  EXPECT_EQ(config1.getSampleCount(), 24);
  EXPECT_EQ(config2.getSampleCount(), 42);  // Should remain unchanged
  
  // Test copy assignment
  Configuration config3;
  config3.setSampleCount(100);
  config3 = config1;
  EXPECT_EQ(config3.getSampleCount(), 24);
}

// Test assignment semantics
TEST_F(PlotConfigurationTest, AssignmentSemantics)
{
  Configuration config1;
  Configuration config2;
  
  config1.setSampleCount(77);
  config2.setSampleCount(88);
  
  // Test assignment
  config2 = config1;
  EXPECT_EQ(config2.getSampleCount(), 77);
  
  // Verify they remain independent after assignment
  config1.setSampleCount(99);
  EXPECT_EQ(config1.getSampleCount(), 99);
  EXPECT_EQ(config2.getSampleCount(), 77);  // Should remain unchanged
}

// Test default values
TEST_F(PlotConfigurationTest, DefaultValues)
{
  Configuration config;
  
  // Default sample count should be positive and reasonable for anti-aliasing
  int defaultSampleCount = config.getSampleCount();
  EXPECT_GT(defaultSampleCount, 0);
  EXPECT_LE(defaultSampleCount, 128);  // Reasonable upper bound for defaults
  
  // Common anti-aliasing values are powers of 2
  bool isPowerOfTwo = (defaultSampleCount & (defaultSampleCount - 1)) == 0;
  EXPECT_TRUE(isPowerOfTwo) << "Default sample count " << defaultSampleCount 
                           << " should typically be a power of 2";
}

// Test multiple instances independence
TEST_F(PlotConfigurationTest, MultipleInstancesIndependence)
{
  Configuration config1;
  Configuration config2;
  Configuration config3;
  
  // Set different values for each
  config1.setSampleCount(4);
  config2.setSampleCount(8);
  config3.setSampleCount(16);
  
  // Verify they maintain their values independently
  EXPECT_EQ(config1.getSampleCount(), 4);
  EXPECT_EQ(config2.getSampleCount(), 8);
  EXPECT_EQ(config3.getSampleCount(), 16);
  
  // Modify one and verify others are unaffected
  config2.setSampleCount(32);
  EXPECT_EQ(config1.getSampleCount(), 4);
  EXPECT_EQ(config2.getSampleCount(), 32);
  EXPECT_EQ(config3.getSampleCount(), 16);
}

// Test const correctness
TEST_F(PlotConfigurationTest, ConstCorrectness)
{
  Configuration config;
  config.setSampleCount(42);
  
  // Test that getSampleCount works with const objects
  const Configuration& constConfig = config;
  EXPECT_EQ(constConfig.getSampleCount(), 42);
  
  // Verify const method doesn't modify the object
  int value1 = constConfig.getSampleCount();
  int value2 = constConfig.getSampleCount();
  EXPECT_EQ(value1, value2);
}