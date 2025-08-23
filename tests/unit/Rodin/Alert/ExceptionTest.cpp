/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <sstream>
#include <string>

#include <Rodin/Alert.h>

using namespace Rodin::Alert;

// Test class for Exception functionality
class ExceptionTest : public ::testing::Test
{
  protected:
    void SetUp() override 
    {
      // Reset string streams for each test
      errorStream.str("");
      errorStream.clear();
    }

    void TearDown() override {}

    std::ostringstream errorStream;
};

// Test Exception Prefix
TEST_F(ExceptionTest, ExceptionPrefix)
{
  ExceptionPrefix prefix;
  std::ostringstream oss;

  // Test that prefix outputs correctly - we can't easily test the color
  // but we can test that it contains the prefix text
  oss << prefix;
  std::string output = oss.str();

  // Should contain "Error" somewhere in the output
  EXPECT_TRUE(output.find("Error") != std::string::npos);
}

// Test Exception class construction (but not raising since it aborts)
TEST_F(ExceptionTest, ExceptionConstruction)
{
  Exception ex(errorStream);
  EXPECT_NO_THROW(ex << "Test exception message");

  // Note: We cannot test raise() because it calls std::abort()
  // which would terminate the test program
}

// Test exception prefix type
TEST_F(ExceptionTest, ExceptionPrefixType)
{
  ExceptionPrefix exPrefix;

  // Test that it can be constructed and used
  std::ostringstream oss;
  EXPECT_NO_THROW(oss << exPrefix);
}