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

// Test class for Success functionality
class SuccessTest : public ::testing::Test
{
  protected:
    void SetUp() override 
    {
      // Reset string streams for each test
      successStream.str("");
      successStream.clear();
    }

    void TearDown() override {}

    std::ostringstream successStream;
};

// Test Success class functionality
TEST_F(SuccessTest, SuccessConstruction)
{
  // Test default construction
  Success success1;

  // Test construction with stream
  Success success2(successStream);

  // Test copy construction
  Success success3(success2);

  // Test move construction
  Success success4(std::move(success3));

  // Basic construction should not throw
  EXPECT_NO_THROW(Success());
}

TEST_F(SuccessTest, SuccessMessage)
{
  Success success(successStream);
  success << "Operation completed successfully";
  success.raise();

  std::string output = successStream.str();

  // Should contain success message
  EXPECT_TRUE(output.find("Operation completed successfully") != std::string::npos);
}

TEST_F(SuccessTest, SuccessStreamOperator)
{
  Success success(successStream);

  // Test streaming various types
  success << "Count: " << 456 << ", Status: " << "OK" << ", Rate: " << 0.95;
  success.raise();

  std::string output = successStream.str();

  EXPECT_TRUE(output.find("Count: 456") != std::string::npos);
  EXPECT_TRUE(output.find("Status: OK") != std::string::npos);
  EXPECT_TRUE(output.find("Rate: 0.95") != std::string::npos);
}