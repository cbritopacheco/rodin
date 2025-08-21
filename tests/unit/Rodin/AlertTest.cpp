/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <sstream>
#include <string>

#include "Rodin/Alert.h"

using namespace Rodin::Alert;

// Test class for Alert functionality
class AlertTest : public ::testing::Test
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

// Test Exception Prefix
TEST_F(AlertTest, ExceptionPrefix)
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

// Test Warning class functionality
TEST_F(AlertTest, WarningConstruction)
{
  // Test default construction
  Warning warn1;

  // Test construction with stream
  Warning warn2(warningStream);

  // Test copy construction
  Warning warn3(warn2);

  // Test move construction
  Warning warn4(std::move(warn3));

  // Basic construction should not throw
  EXPECT_NO_THROW(Warning());
}

TEST_F(AlertTest, WarningMessage)
{
  Warning warn(warningStream);
  warn << "This is a test warning message";
  warn.raise();

  std::string output = warningStream.str();

  // Should contain warning prefix and message
  EXPECT_TRUE(output.find("Warning") != std::string::npos);
  EXPECT_TRUE(output.find("This is a test warning message") != std::string::npos);
}

TEST_F(AlertTest, WarningStreamOperator)
{
  Warning warn(warningStream);

  // Test streaming various types
  warn << "Number: " << 42 << ", Float: " << 3.14 << ", Bool: " << true;
  warn.raise();

  std::string output = warningStream.str();

  EXPECT_TRUE(output.find("Number: 42") != std::string::npos);
  EXPECT_TRUE(output.find("Float: 3.14") != std::string::npos);
  EXPECT_TRUE(output.find("Bool: 1") != std::string::npos);
}

// Test Info class functionality
TEST_F(AlertTest, InfoConstruction)
{
  // Test default construction
  Info info1;

  // Test construction with stream
  Info info2(infoStream);

  // Test copy construction
  Info info3(info2);

  // Test move construction
  Info info4(std::move(info3));

  // Basic construction should not throw
  EXPECT_NO_THROW(Info());
}

TEST_F(AlertTest, InfoMessage)
{
  Info info(infoStream);
  info << "This is an informational message";
  info.raise();

  std::string output = infoStream.str();

  // Should contain info prefix and message
  EXPECT_TRUE(output.find("Info") != std::string::npos);
  EXPECT_TRUE(output.find("This is an informational message") != std::string::npos);
}

TEST_F(AlertTest, InfoStreamOperator)
{
  Info info(infoStream);

  // Test streaming various types
  info << "Value: " << 123 << ", Text: " << "hello" << ", Double: " << 2.718;
  info.raise();

  std::string output = infoStream.str();

  EXPECT_TRUE(output.find("Value: 123") != std::string::npos);
  EXPECT_TRUE(output.find("Text: hello") != std::string::npos);
  EXPECT_TRUE(output.find("Double: 2.718") != std::string::npos);
}

// Test Success class functionality
TEST_F(AlertTest, SuccessConstruction)
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

TEST_F(AlertTest, SuccessMessage)
{
  Success success(successStream);
  success << "Operation completed successfully";
  success.raise();

  std::string output = successStream.str();

  // Should contain success message
  EXPECT_TRUE(output.find("Operation completed successfully") != std::string::npos);
}

TEST_F(AlertTest, SuccessStreamOperator)
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

// Test message chaining
TEST_F(AlertTest, MessageChaining)
{
  Warning warn(warningStream);

  // Test method chaining
  warn << "First part" << " - " << "Second part" << " - " << 42;
  warn.raise();

  std::string output = warningStream.str();

  EXPECT_TRUE(output.find("First part - Second part - 42") != std::string::npos);
}

// Test empty messages
TEST_F(AlertTest, EmptyMessages)
{
  Warning warn(warningStream);
  warn.raise(); // Empty message

  std::string output = warningStream.str();

  // Should still contain the prefix even with empty message
  EXPECT_TRUE(output.find("Warning") != std::string::npos);
}

// Test multiple raises
TEST_F(AlertTest, MultipleRaises)
{
  Info info(infoStream);

  info << "First message";
  info.raise();

  std::string output1 = infoStream.str();
  EXPECT_TRUE(output1.find("First message") != std::string::npos);

  // Clear and add more content
  infoStream.str("");
  infoStream.clear();

  info << "Second message";
  info.raise();

  std::string output2 = infoStream.str();
  EXPECT_TRUE(output2.find("Second message") != std::string::npos);
}

// Test Exception class construction (but not raising since it aborts)
TEST_F(AlertTest, ExceptionConstruction)
{
  Exception ex(errorStream);
  EXPECT_NO_THROW(ex << "Test exception message");

  // Note: We cannot test raise() because it calls std::abort()
  // which would terminate the test program
}

// Test prefixes have correct types
TEST_F(AlertTest, PrefixTypes)
{
  ExceptionPrefix exPrefix;
  WarningPrefix warnPrefix; 
  InfoPrefix infoPrefix;

  // Test that they can be constructed and used
  std::ostringstream oss;
  EXPECT_NO_THROW(oss << exPrefix);
  EXPECT_NO_THROW(oss << warnPrefix);
  EXPECT_NO_THROW(oss << infoPrefix);
}

// Test color functionality (basic test)
TEST_F(AlertTest, ColorOutputPresence)
{
  Warning warn(warningStream);
  warn << "Colored message";
  warn.raise();

  std::string output = warningStream.str();

  // The output should be non-empty and contain our message
  EXPECT_FALSE(output.empty());
  EXPECT_TRUE(output.find("Colored message") != std::string::npos);

  // May contain ANSI color codes (depending on implementation)
  // but we can't easily test for specific codes without knowing the implementation
}

// Integration test with multiple alert types
TEST_F(AlertTest, MultipleAlertTypes)
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
