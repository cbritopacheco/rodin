/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Test/Utility/ExpectedResultTable.h"
#include "Rodin/Types.h"

using namespace Rodin;
using namespace Rodin::Test::Utility;

class ExpectedResultTableTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// Test basic construction and single parameter function
TEST_F(ExpectedResultTableTest, BasicConstruction)
{
  // Test with single parameter function
  std::function<int(int)> squareFunc = [](int x) { return x * x; };
  ExpectedResultTable<int, int> table(squareFunc);

  // Should construct without throwing
  EXPECT_NO_THROW(table.emplace_back(4, 2));
  EXPECT_NO_THROW(table.emplace_back(9, 3));
  EXPECT_NO_THROW(table.emplace_back(16, 4));
}

// Test ExpectedResult class
TEST_F(ExpectedResultTableTest, ExpectedResult)
{
  using ExpectedResult = ExpectedResultTable<Real, Real, Real>::ExpectedResult;

  // Test construction
  ExpectedResult entry(5.0, 2.0, 3.0);

  // Test getResult
  EXPECT_DOUBLE_EQ(entry.getResult(), 5.0);

  // Test getParameters
  auto params = entry.getParameters();
  EXPECT_DOUBLE_EQ(std::get<0>(params), 2.0);
  EXPECT_DOUBLE_EQ(std::get<1>(params), 3.0);
}

// Test two-parameter function evaluation
TEST_F(ExpectedResultTableTest, TwoParameterFunction)
{
  std::function<Real(Real, Real)> addFunc = [](Real a, Real b) { return a + b; };
  ExpectedResultTable<Real, Real, Real> table(addFunc);

  // Add expected results
  table.emplace_back(5.0, 2.0, 3.0);  // 2 + 3 = 5
  table.emplace_back(7.0, 3.0, 4.0);  // 3 + 4 = 7
  table.emplace_back(0.0, -1.0, 1.0); // -1 + 1 = 0

  // All should pass
  EXPECT_TRUE(table.evaluate());
}

// Test failed evaluation
TEST_F(ExpectedResultTableTest, FailedEvaluation)
{
  std::function<Real(Real, Real)> addFunc = [](Real a, Real b) { return a + b; };
  ExpectedResultTable<Real, Real, Real> table(addFunc);

  // Add correct and incorrect expected results
  table.emplace_back(5.0, 2.0, 3.0);  // 2 + 3 = 5 (correct)
  table.emplace_back(8.0, 3.0, 4.0);  // 3 + 4 = 7, but expecting 8 (incorrect)

  // Should fail because one result is wrong
  EXPECT_FALSE(table.evaluate());
}

// Test with multiplication function
TEST_F(ExpectedResultTableTest, MultiplicationFunction)
{
  std::function<Real(Real, Real)> multiplyFunc = [](Real a, Real b) { return a * b; };
  ExpectedResultTable<Real, Real, Real> table(multiplyFunc);

  // Add expected results for multiplication
  table.emplace_back(6.0, 2.0, 3.0);   // 2 * 3 = 6
  table.emplace_back(12.0, 3.0, 4.0);  // 3 * 4 = 12
  table.emplace_back(0.0, 0.0, 5.0);   // 0 * 5 = 0
  table.emplace_back(-6.0, -2.0, 3.0); // -2 * 3 = -6

  // All should pass
  EXPECT_TRUE(table.evaluate());
}

// Test with custom comparison function
TEST_F(ExpectedResultTableTest, CustomComparison)
{
  std::function<Real(Real, Real)> addFunc = [](Real a, Real b) { return a + b; };

  // Custom comparison with tolerance
  auto toleranceCompare = [](const Real& modelResult, const Real& expectedResult) {
    return std::abs(modelResult - expectedResult) < 1e-10;
  };

  ExpectedResultTable<Real, Real, Real> table(addFunc, toleranceCompare);

  // Add results that might have small floating point errors
  table.emplace_back(0.3, 0.1, 0.2);  // 0.1 + 0.2 might not exactly equal 0.3

  // Should pass with tolerance comparison
  EXPECT_TRUE(table.evaluate());
}

// Test with integer types
TEST_F(ExpectedResultTableTest, IntegerTypes)
{
  std::function<int(int)> squareFunc = [](int x) { return x * x; };
  ExpectedResultTable<int, int> table(squareFunc);

  // Add expected results for squares
  table.emplace_back(0, 0);    // 0^2 = 0
  table.emplace_back(1, 1);    // 1^2 = 1
  table.emplace_back(4, 2);    // 2^2 = 4
  table.emplace_back(9, 3);    // 3^2 = 9
  table.emplace_back(16, 4);   // 4^2 = 16
  table.emplace_back(25, -5);  // (-5)^2 = 25

  // All should pass
  EXPECT_TRUE(table.evaluate());
}

// Test with string types
TEST_F(ExpectedResultTableTest, StringTypes)
{
  std::function<std::string(std::string, std::string)> concatFunc = 
    [](const std::string& a, const std::string& b) { return a + b; };
  ExpectedResultTable<std::string, std::string, std::string> table(concatFunc);

  // Add expected results for string concatenation
  table.emplace_back(std::string("hello world"), std::string("hello "), std::string("world"));
  table.emplace_back(std::string("ab"), std::string("a"), std::string("b"));
  table.emplace_back(std::string("test"), std::string("test"), std::string(""));
  table.emplace_back(std::string(""), std::string(""), std::string(""));

  // All should pass
  EXPECT_TRUE(table.evaluate());
}

// Test with three parameters
TEST_F(ExpectedResultTableTest, ThreeParameters)
{
  std::function<Real(Real, Real, Real)> sumThreeFunc = 
    [](Real a, Real b, Real c) { return a + b + c; };
  ExpectedResultTable<Real, Real, Real, Real> table(sumThreeFunc);

  // Add expected results for sum of three numbers
  table.emplace_back(6.0, 1.0, 2.0, 3.0);    // 1 + 2 + 3 = 6
  table.emplace_back(0.0, -1.0, 0.0, 1.0);   // -1 + 0 + 1 = 0
  table.emplace_back(15.0, 5.0, 5.0, 5.0);   // 5 + 5 + 5 = 15

  // All should pass
  EXPECT_TRUE(table.evaluate());
}

// Test push_back with ExpectedResult objects
TEST_F(ExpectedResultTableTest, PushBackExpectedResult)
{
  using ExpectedResult = ExpectedResultTable<Real, Real, Real>::ExpectedResult;

  std::function<Real(Real, Real)> addFunc = [](Real a, Real b) { return a + b; };
  ExpectedResultTable<Real, Real, Real> table(addFunc);

  // Create ExpectedResult objects and push them
  ExpectedResult entry1(5.0, 2.0, 3.0);
  ExpectedResult entry2(7.0, 3.0, 4.0);

  table.push_back(entry1);
  table.push_back(entry2);

  // Should evaluate correctly
  EXPECT_TRUE(table.evaluate());
}

// Test with lambda functions
TEST_F(ExpectedResultTableTest, LambdaFunctions)
{
  // Test with a more complex lambda
  auto polynomialFunc = [](Real x) { return x*x*x - 2*x*x + x - 1; };
  std::function<Real(Real)> polyFunc = polynomialFunc;
  ExpectedResultTable<Real, Real> table(polyFunc);

  // Add some known results
  table.emplace_back(-1.0, 0.0);  // 0^3 - 2*0^2 + 0 - 1 = -1
  table.emplace_back(-1.0, 1.0);  // 1^3 - 2*1^2 + 1 - 1 = -1
  table.emplace_back(1.0, 2.0);   // 2^3 - 2*2^2 + 2 - 1 = 8 - 8 + 2 - 1 = 1

  // All should pass
  EXPECT_TRUE(table.evaluate());
}

// Test error cases and edge conditions
TEST_F(ExpectedResultTableTest, EdgeCases)
{
  std::function<Real(Real, Real)> addFunc = [](Real a, Real b) { return a + b; };
  ExpectedResultTable<Real, Real, Real> table(addFunc);

  // Test with empty table
  EXPECT_TRUE(table.evaluate());  // Empty table should pass

  // Test with very small numbers
  table.emplace_back(1e-15, 1e-15, 0.0);
  EXPECT_TRUE(table.evaluate());

  // Test with very large numbers
  ExpectedResultTable<Real, Real, Real> largeTable(addFunc);
  largeTable.emplace_back(2e10, 1e10, 1e10);
  EXPECT_TRUE(largeTable.evaluate());
}

// Test with boolean return type
TEST_F(ExpectedResultTableTest, BooleanReturnType)
{
  std::function<bool(int, int)> greaterThanFunc = 
    [](int a, int b) { return a > b; };
  ExpectedResultTable<bool, int, int> table(greaterThanFunc);

  // Add expected boolean results
  table.emplace_back(true, 5, 3);   // 5 > 3 = true
  table.emplace_back(false, 2, 4);  // 2 > 4 = false
  table.emplace_back(false, 3, 3);  // 3 > 3 = false
  table.emplace_back(true, 10, -5); // 10 > -5 = true

  // All should pass
  EXPECT_TRUE(table.evaluate());
}

// Test mixed success and failure
TEST_F(ExpectedResultTableTest, MixedResults)
{
  std::function<int(int, int)> subtractFunc = [](int a, int b) { return a - b; };
  ExpectedResultTable<int, int, int> table(subtractFunc);

  // Mix of correct and incorrect results
  table.emplace_back(2, 5, 3);   // 5 - 3 = 2 (correct)
  table.emplace_back(1, 4, 2);   // 4 - 2 = 2, but expecting 1 (incorrect)
  table.emplace_back(-1, 3, 4);  // 3 - 4 = -1 (correct)

  // Should fail because of the incorrect middle result
  EXPECT_FALSE(table.evaluate());
}

// Test performance with many entries
TEST_F(ExpectedResultTableTest, PerformanceTest)
{
  std::function<Real(Real)> identityFunc = [](Real x) { return x; };
  ExpectedResultTable<Real, Real> table(identityFunc);

  // Add many entries
  for (int i = 0; i < 1000; ++i)
  {
    Real val = static_cast<Real>(i) / 10.0;
    table.emplace_back(val, val);
  }

  // Should evaluate quickly and correctly
  EXPECT_TRUE(table.evaluate());
}
