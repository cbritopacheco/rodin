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

// Test class for Warning functionality
class WarningTest : public ::testing::Test
{
  protected:
    void SetUp() override 
    {
      // Reset string streams for each test
      warningStream.str("");
      warningStream.clear();
    }

    void TearDown() override {}

    std::ostringstream warningStream;
};

// Test Warning class functionality
TEST_F(WarningTest, WarningConstruction)
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

TEST_F(WarningTest, WarningMessage)
{
  Warning warn(warningStream);
  warn << "This is a test warning message";
  warn.raise();

  std::string output = warningStream.str();

  // Should contain warning prefix and message
  EXPECT_TRUE(output.find("Warning") != std::string::npos);
  EXPECT_TRUE(output.find("This is a test warning message") != std::string::npos);
}

TEST_F(WarningTest, WarningStreamOperator)
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

// Test message chaining
TEST_F(WarningTest, MessageChaining)
{
  Warning warn(warningStream);

  // Test method chaining
  warn << "First part" << " - " << "Second part" << " - " << 42;
  warn.raise();

  std::string output = warningStream.str();

  EXPECT_TRUE(output.find("First part - Second part - 42") != std::string::npos);
}

// Test empty messages
TEST_F(WarningTest, EmptyMessages)
{
  Warning warn(warningStream);
  warn.raise(); // Empty message

  std::string output = warningStream.str();

  // Should still contain the prefix even with empty message
  EXPECT_TRUE(output.find("Warning") != std::string::npos);
}

// Test color functionality (basic test)
TEST_F(WarningTest, ColorOutputPresence)
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

// Test warning prefix type
TEST_F(WarningTest, WarningPrefixType)
{
  WarningPrefix warnPrefix; 

  // Test that it can be constructed and used
  std::ostringstream oss;
  EXPECT_NO_THROW(oss << warnPrefix);
}