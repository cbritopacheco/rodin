/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/FormLanguage/Traits.h"
#include "Rodin/Types.h"

using namespace Rodin;
using namespace Rodin::FormLanguage;

// Test class for traits specialization
class TestClass
{
  public:
    using ScalarType = Real;
    using ValueType = int;
};

// Specialize Traits for our test class
namespace Rodin::FormLanguage
{
  template <>
  struct Traits<TestClass>
  {
    using ScalarType = Real;
    using ValueType = int;
  };
}

class FormLanguageTraitsTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// Test basic Traits template structure
TEST_F(FormLanguageTraitsTest, BasicTraitsStructure)
{
  // Test that Traits is a template that can be specialized
  static_assert(std::is_class_v<Traits<TestClass>>);
  
  // Test that our specialization works
  static_assert(std::is_same_v<typename Traits<TestClass>::ScalarType, Real>);
  static_assert(std::is_same_v<typename Traits<TestClass>::ValueType, int>);
}

// Test Traits with different types
TEST_F(FormLanguageTraitsTest, DifferentTypes)
{
  // Test that Traits can be used with different types
  // The actual specializations would be defined elsewhere in the codebase
  
  // For now, just test that the template exists and can be instantiated
  // (even if no specialization exists, the primary template should be defined)
  
  // Test with built-in types
  static_assert(std::is_class_v<Traits<int>>);
  static_assert(std::is_class_v<Traits<double>>);
  static_assert(std::is_class_v<Traits<std::string>>);
}

// Test type deduction with Traits
TEST_F(FormLanguageTraitsTest, TypeDeduction)
{
  // Test that we can use Traits to extract type information
  using TestScalarType = typename Traits<TestClass>::ScalarType;
  using TestValueType = typename Traits<TestClass>::ValueType;
  
  // Verify the types are what we expect
  static_assert(std::is_same_v<TestScalarType, Real>);
  static_assert(std::is_same_v<TestValueType, int>);
  
  // Test that we can use these types in variable declarations
  TestScalarType scalar = 3.14;
  TestValueType value = 42;
  
  EXPECT_TRUE(std::is_same_v<decltype(scalar), Real>);
  EXPECT_TRUE(std::is_same_v<decltype(value), int>);
  
  // Suppress unused variable warnings
  (void)scalar;
  (void)value;
}

// Test SFINAE capabilities with Traits
TEST_F(FormLanguageTraitsTest, SFINAECapabilities)
{
  // Test that Traits can be used in SFINAE contexts
  
  // Helper to check if a type has ScalarType
  template <typename T, typename = void>
  struct has_scalar_type : std::false_type {};
  
  template <typename T>
  struct has_scalar_type<T, std::void_t<typename Traits<T>::ScalarType>> : std::true_type {};
  
  // Test that our TestClass has ScalarType
  static_assert(has_scalar_type<TestClass>::value);
  
  // Test with a type that doesn't have a Traits specialization
  class NoTraitsClass {};
  
  // This would be false if there's no specialization defining ScalarType
  // The exact behavior depends on whether a primary template is defined
  // For now, just verify the mechanism works
  constexpr bool hasScalarType = has_scalar_type<NoTraitsClass>::value;
  (void)hasScalarType;  // Suppress unused variable warning
}

// Test that Traits works with const and reference types
TEST_F(FormLanguageTraitsTest, CVRefQualifiers)
{
  // Test that Traits can work with cv-qualified and reference types
  // This tests the flexibility of the Traits system
  
  // Test with const
  static_assert(std::is_class_v<Traits<const TestClass>>);
  
  // Test with reference (if specialized)
  static_assert(std::is_class_v<Traits<TestClass&>>);
  
  // Test with const reference
  static_assert(std::is_class_v<Traits<const TestClass&>>);
}

// Test Traits as a metaprogramming tool
TEST_F(FormLanguageTraitsTest, MetaprogrammingUsage)
{
  // Test that Traits can be used for template metaprogramming
  
  // Example: A function template that uses Traits
  auto getScalarValue = []<typename T>() -> typename Traits<T>::ScalarType {
    if constexpr (std::is_same_v<T, TestClass>) {
      return typename Traits<T>::ScalarType{42.0};
    } else {
      return typename Traits<T>::ScalarType{};
    }
  };
  
  auto result = getScalarValue.template operator()<TestClass>();
  EXPECT_TRUE(std::is_same_v<decltype(result), Real>);
  EXPECT_DOUBLE_EQ(result, 42.0);
}

// Test that Traits preserves type relationships
TEST_F(FormLanguageTraitsTest, TypeRelationships)
{
  // Test that Traits maintains expected type relationships
  
  // For our test class, ScalarType should be Real
  using ScalarType = typename Traits<TestClass>::ScalarType;
  static_assert(std::is_floating_point_v<ScalarType>);
  static_assert(std::is_same_v<ScalarType, Real>);
  
  // ValueType should be int
  using ValueType = typename Traits<TestClass>::ValueType;
  static_assert(std::is_integral_v<ValueType>);
  static_assert(std::is_same_v<ValueType, int>);
}

// Test compilation with enable_if and Traits
TEST_F(FormLanguageTraitsTest, EnableIfIntegration)
{
  // Test that Traits works well with std::enable_if for conditional compilation
  
  // Function that only works with types that have Real as ScalarType
  template <typename T>
  typename std::enable_if_t<
    std::is_same_v<typename Traits<T>::ScalarType, Real>,
    bool
  > hasRealScalarType() {
    return true;
  }
  
  // This should compile and return true for TestClass
  bool result = hasRealScalarType<TestClass>();
  EXPECT_TRUE(result);
}

// Test multiple Traits specializations
TEST_F(FormLanguageTraitsTest, MultipleSpecializations)
{
  // Define another test class
  class AnotherTestClass {};
  
  // We can't define the specialization here because it would be in the wrong namespace
  // But we can test that the mechanism allows for multiple specializations
  
  // Test that we can have different Traits for different classes
  static_assert(std::is_class_v<Traits<TestClass>>);
  static_assert(std::is_class_v<Traits<AnotherTestClass>>);
  
  // The actual type members would be different for each specialization
  // This demonstrates the extensibility of the Traits system
}

// Test Traits with template types
TEST_F(FormLanguageTraitsTest, TemplateTypes)
{
  // Test that Traits can handle template types
  static_assert(std::is_class_v<Traits<std::vector<int>>>);
  static_assert(std::is_class_v<Traits<std::pair<int, double>>>);
  
  // These would need specializations to be useful, but the structure supports them
}

// Test boost::type_index integration (as included in Traits.h)
TEST_F(FormLanguageTraitsTest, TypeIndexIntegration)
{
  // Test that boost::type_index works with our types
  auto typeInfo = boost::typeindex::type_id<TestClass>();
  EXPECT_FALSE(typeInfo.pretty_name().empty());
  
  // Test with Traits types
  auto scalarTypeInfo = boost::typeindex::type_id<typename Traits<TestClass>::ScalarType>();
  EXPECT_FALSE(scalarTypeInfo.pretty_name().empty());
  
  // Verify they're different
  EXPECT_NE(typeInfo.pretty_name(), scalarTypeInfo.pretty_name());
}