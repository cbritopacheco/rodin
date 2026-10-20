/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <sstream>
#include <string>
#include <utility>

#include <Rodin/Alert.h>
#include <Rodin/Alert/ClassException.h>
#include <Rodin/Alert/MemberFunctionException.h>
#include <Rodin/Alert/NamespacedException.h>

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
/// @brief Verifies exception prefix for exception test by checking true predicates.
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

// Test Exception class construction without raising.
/// @brief Verifies exception construction for exception test by checking no-throw behavior.
TEST_F(ExceptionTest, ExceptionConstruction)
{
  Exception ex(errorStream);
  EXPECT_NO_THROW(ex << "Test exception message");

  // Raising and propagation are exercised by the contextual cases below.
}

// Test exception prefix type
/// @brief Verifies exception prefix type for exception test by checking no-throw behavior.
TEST_F(ExceptionTest, ExceptionPrefixType)
{
  ExceptionPrefix exPrefix;

  // Test that it can be constructed and used
  std::ostringstream oss;
  EXPECT_NO_THROW(oss << exPrefix);
}

namespace Rodin::Tests::Unit
{
  struct AlertExceptionContext
  {};

  /// @brief Class context survives streaming and propagation as a native exception.
  TEST(AlertContextException, ClassMessageAndRaise)
  {
    const AlertExceptionContext context;
    ClassException error(context);
    std::ostringstream output;
    error.setOutputStream(output);
    error << "invalid index " << 17;
    const std::string expected =
      "In class Rodin::Tests::Unit::AlertExceptionContext. invalid index 17";
    EXPECT_STREQ(error.what(), expected.c_str());
    try
    {
      error << Raise;
      FAIL() << "The Alert raise tag must propagate an exception.";
    }
    catch (const Exception& caught)
    {
      const std::exception& standardView = caught;
      EXPECT_STREQ(standardView.what(), expected.c_str());
    }
    EXPECT_NE(output.str().find("invalid index 17"), std::string::npos);
    EXPECT_NE(output.str().find("Error"), std::string::npos);
  }

  /// @brief Literal and owned function names produce the same member context.
  TEST(AlertContextException, MemberFunctionNameVariants)
  {
    AlertExceptionContext context;
    MemberFunctionException literal(context, "assemble");
    const std::string function = "assemble";
    MemberFunctionException owned(std::as_const(context), function);
    literal << "invalid shape";
    owned << "invalid shape";
    const std::string expected =
      "In member function assemble of class Rodin::Tests::Unit::AlertExceptionContext: "
      "invalid shape";
    EXPECT_STREQ(literal.what(), expected.c_str());
    EXPECT_STREQ(owned.what(), expected.c_str());
    std::ostringstream output;
    owned.setOutputStream(output);
    const Exception& base = owned;
    try
    {
      base.raise();
      FAIL() << "Virtual dispatch must retain the native exception contract.";
    }
    catch (const Exception& caught)
    {
      EXPECT_STREQ(caught.what(), expected.c_str());
    }
    EXPECT_NE(output.str().find("assemble"), std::string::npos);
    EXPECT_NE(output.str().find("invalid shape"), std::string::npos);
  }

  /// @brief Namespace context is retained in the payload and captured output.
  TEST(AlertContextException, NamespaceMessageAndRaise)
  {
    NamespacedException error("Rodin::Geometry");
    std::ostringstream output;
    error.setOutputStream(output);
    error << "invalid dimension " << 4;
    const std::string expected = "Rodin::Geometry. invalid dimension 4";
    EXPECT_STREQ(error.what(), expected.c_str());
    try
    {
      error << Raise;
      FAIL() << "Namespaced errors must propagate as native Alert exceptions.";
    }
    catch (const Exception& caught)
    {
      EXPECT_STREQ(caught.what(), expected.c_str());
    }
    EXPECT_NE(output.str().find("Rodin::Geometry"), std::string::npos);
    EXPECT_NE(output.str().find("invalid dimension 4"), std::string::npos);
  }
}
