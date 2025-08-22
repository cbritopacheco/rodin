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

// Test class for Alert integration functionality
class AlertIntegrationTest : public ::testing::Test
{
  protected:
    void SetUp() override 
    {
      // Reset string streams for each test
      errorStream.str("");
      errorStream.clear();
      warningStream.str("");
      warningStream.clear();
      infoStream.str("");
      infoStream.clear();
      successStream.str("");
      successStream.clear();
    }

    void TearDown() override {}

    std::ostringstream errorStream;
    std::ostringstream warningStream;
    std::ostringstream infoStream;
    std::ostringstream successStream;
};

// Integration test with multiple alert types
TEST_F(AlertIntegrationTest, MultipleAlertTypes)
{
  // Create multiple alerts with different streams
  Warning warn(warningStream);
  Info info(infoStream);
  Success success(successStream);

  warn << "Warning message";
  info << "Info message";
  success << "Success message";

  warn.raise();
  info.raise();
  success.raise();

  EXPECT_TRUE(warningStream.str().find("Warning message") != std::string::npos);
  EXPECT_TRUE(infoStream.str().find("Info message") != std::string::npos);
  EXPECT_TRUE(successStream.str().find("Success message") != std::string::npos);
}