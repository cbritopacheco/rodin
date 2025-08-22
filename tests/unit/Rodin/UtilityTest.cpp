/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <variant>
#include <string>
#include <vector>
#include <type_traits>

#include "Rodin/Utility.h"

using namespace Rodin::Utility;

// Test class for Overloaded functionality
class UtilityOverloadedTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(UtilityOverloadedTest, VariantVisitation)
{
  std::variant<int, std::string, double> var1 = 42;
  std::variant<int, std::string, double> var2 = std::string("hello");
  std::variant<int, std::string, double> var3 = 3.14;

  auto visitor = Overloaded {
    [](int i) -> std::string { return "integer: " + std::to_string(i); },
    [](const std::string& s) -> std::string { return "string: " + s; },
    [](double d) -> std::string { return "double: " + std::to_string(d); }
  };

  EXPECT_EQ(std::visit(visitor, var1), "integer: 42");
  EXPECT_EQ(std::visit(visitor, var2), "string: hello");
  EXPECT_EQ(std::visit(visitor, var3), "double: 3.140000");
}

TEST_F(UtilityOverloadedTest, MultipleArgumentTypes)
{
  auto visitor = Overloaded {
    [](int i) { return i * 2; },
    [](double d) { return d * 2.0; },
    [](const std::string& s) { return s + s; }
  };

  EXPECT_EQ(visitor(5), 10);
  EXPECT_DOUBLE_EQ(visitor(2.5), 5.0);
  EXPECT_EQ(visitor(std::string("test")), "testtest");
}

TEST_F(UtilityOverloadedTest, NestedOverloads)
{
  using NestedVariant = std::variant<int, std::variant<std::string, double>>;
  
  NestedVariant nested1 = 42;
  NestedVariant nested2 = std::variant<std::string, double>(std::string("nested"));
  NestedVariant nested3 = std::variant<std::string, double>(1.23);

  auto innerVisitor = Overloaded {
    [](const std::string& s) -> std::string { return "inner_string: " + s; },
    [](double d) -> std::string { return "inner_double: " + std::to_string(d); }
  };

  auto outerVisitor = Overloaded {
    [](int i) -> std::string { return "outer_int: " + std::to_string(i); },
    [&innerVisitor](const std::variant<std::string, double>& inner) -> std::string {
      return std::visit(innerVisitor, inner);
    }
  };

  EXPECT_EQ(std::visit(outerVisitor, nested1), "outer_int: 42");
  EXPECT_EQ(std::visit(outerVisitor, nested2), "inner_string: nested");
  EXPECT_EQ(std::visit(outerVisitor, nested3), "inner_double: 1.230000");
}

// Test class for IsSpecialization functionality
class UtilityIsSpecializationTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(UtilityIsSpecializationTest, VectorSpecialization)
{
  EXPECT_TRUE((IsSpecialization<std::vector<int>, std::vector>::Value));
  EXPECT_TRUE((IsSpecialization<std::vector<double>, std::vector>::Value));
  EXPECT_TRUE((IsSpecialization<std::vector<std::string>, std::vector>::Value));
  
  EXPECT_FALSE((IsSpecialization<int, std::vector>::Value));
  EXPECT_FALSE((IsSpecialization<std::string, std::vector>::Value));
}

TEST_F(UtilityIsSpecializationTest, VariantSpecialization)
{
  EXPECT_TRUE((IsSpecialization<std::variant<int>, std::variant>::Value));
  EXPECT_TRUE((IsSpecialization<std::variant<int, double>, std::variant>::Value));
  EXPECT_TRUE((IsSpecialization<std::variant<int, std::string, double>, std::variant>::Value));
  
  EXPECT_FALSE((IsSpecialization<int, std::variant>::Value));
  EXPECT_FALSE((IsSpecialization<std::vector<int>, std::variant>::Value));
}

// Custom template for testing
template<typename T, typename U = int>
struct CustomTemplate {};

TEST_F(UtilityIsSpecializationTest, CustomTemplateSpecialization)
{
  EXPECT_TRUE((IsSpecialization<CustomTemplate<int>, CustomTemplate>::Value));
  EXPECT_TRUE((IsSpecialization<CustomTemplate<double>, CustomTemplate>::Value));
  EXPECT_TRUE((IsSpecialization<CustomTemplate<int, double>, CustomTemplate>::Value));
  
  EXPECT_FALSE((IsSpecialization<int, CustomTemplate>::Value));
  EXPECT_FALSE((IsSpecialization<std::vector<int>, CustomTemplate>::Value));
}

TEST_F(UtilityIsSpecializationTest, NestedTemplateSpecialization)
{
  EXPECT_TRUE((IsSpecialization<std::vector<std::vector<int>>, std::vector>::Value));
  EXPECT_TRUE((IsSpecialization<CustomTemplate<std::vector<int>>, CustomTemplate>::Value));
  EXPECT_TRUE((IsSpecialization<std::vector<CustomTemplate<int>>, std::vector>::Value));
  
  // The inner template doesn't match the outer check
  EXPECT_FALSE((IsSpecialization<std::vector<CustomTemplate<int>>, CustomTemplate>::Value));
  EXPECT_FALSE((IsSpecialization<CustomTemplate<std::vector<int>>, std::vector>::Value));
}

// Test class for False functionality
class UtilityFalseTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(UtilityFalseTest, BasicFalseBehavior)
{
  EXPECT_FALSE(False<>);
  EXPECT_FALSE((False<int>));
  EXPECT_FALSE((False<double>));
  EXPECT_FALSE((False<std::string>));
}

TEST_F(UtilityFalseTest, MultipleTypesFalse)
{
  EXPECT_FALSE((False<int, double>));
  EXPECT_FALSE((False<int, double, std::string>));
  EXPECT_FALSE((False<std::vector<int>, CustomTemplate<double>>));
}

TEST_F(UtilityFalseTest, CompileTimeConstant)
{
  static_assert(False<> == false);
  static_assert(False<int> == false);
  static_assert(False<int, double> == false);
  
  // Test that it's truly constexpr
  constexpr bool result1 = False<>;
  constexpr bool result2 = False<int>;
  constexpr bool result3 = False<int, double, std::string>;
  
  EXPECT_FALSE(result1);
  EXPECT_FALSE(result2);
  EXPECT_FALSE(result3);
}

// Test usage of False in template metaprogramming contexts
template<typename T>
constexpr bool always_false_helper()
{
  return False<T>;
}

TEST_F(UtilityFalseTest, TemplateMetaprogrammingUsage)
{
  EXPECT_FALSE(always_false_helper<int>());
  EXPECT_FALSE(always_false_helper<double>());
  EXPECT_FALSE(always_false_helper<std::string>());
}

// Integration test combining multiple utilities
class UtilityIntegrationTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(UtilityIntegrationTest, CombinedUtilities)
{
  // Use IsSpecialization to detect vector types and False for default case
  auto typeChecker = Overloaded {
    [](const auto& value) -> std::string {
      using ValueType = std::decay_t<decltype(value)>;
      if constexpr (IsSpecialization<ValueType, std::vector>::Value) {
        return "vector type";
      } else if constexpr (std::is_integral_v<ValueType>) {
        return "integral type";
      } else if constexpr (False<ValueType>) {
        // This branch will never be taken due to False<T> being false
        return "impossible";
      } else {
        return "other type";
      }
    }
  };

  std::vector<int> vec{1, 2, 3};
  int integer = 42;
  double floating = 3.14;
  std::string text = "hello";

  EXPECT_EQ(typeChecker(vec), "vector type");
  EXPECT_EQ(typeChecker(integer), "integral type");
  EXPECT_EQ(typeChecker(floating), "other type");
  EXPECT_EQ(typeChecker(text), "other type");
}

TEST_F(UtilityIntegrationTest, VariantWithSpecializationDetection)
{
  using TestVariant = std::variant<std::vector<int>, int, std::string>;
  
  TestVariant var1 = std::vector<int>{1, 2, 3};
  TestVariant var2 = 42;
  TestVariant var3 = std::string("test");

  auto visitor = Overloaded {
    [](const auto& value) -> std::string {
      using ValueType = std::decay_t<decltype(value)>;
      if constexpr (IsSpecialization<ValueType, std::vector>::Value) {
        return "Found vector with " + std::to_string(value.size()) + " elements";
      } else if constexpr (std::is_same_v<ValueType, int>) {
        return "Found integer: " + std::to_string(value);
      } else {
        return "Found other: " + std::string(typeid(ValueType).name());
      }
    }
  };

  EXPECT_EQ(std::visit(visitor, var1), "Found vector with 3 elements");
  EXPECT_EQ(std::visit(visitor, var2), "Found integer: 42");
  EXPECT_TRUE(std::visit(visitor, var3).find("Found other:") == 0);
}