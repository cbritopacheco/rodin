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

// Test class for Info functionality
class InfoTest : public ::testing::Test
{
  protected:
    void SetUp() override 
    {
      // Reset string streams for each test
      infoStream.str("");
      infoStream.clear();
    }

    void TearDown() override {}

    std::ostringstream infoStream;
};

// Test Info class functionality
TEST_F(InfoTest, InfoConstruction)
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

TEST_F(InfoTest, InfoMessage)
{
  Info info(infoStream);
  info << "This is an informational message";
  info.raise();

  std::string output = infoStream.str();

  // Should contain info prefix and message
  EXPECT_TRUE(output.find("Info") != std::string::npos);
  EXPECT_TRUE(output.find("This is an informational message") != std::string::npos);
}

TEST_F(InfoTest, InfoStreamOperator)
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

// Test multiple raises
TEST_F(InfoTest, MultipleRaises)
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

// Test info prefix type
TEST_F(InfoTest, InfoPrefixType)
{
  InfoPrefix infoPrefix;

  // Test that it can be constructed and used
  std::ostringstream oss;
  EXPECT_NO_THROW(oss << infoPrefix);
}