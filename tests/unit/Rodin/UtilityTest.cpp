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
#include <optional>
#include <functional>

#include "Rodin/Utility.h"
#include "Rodin/Utility/Zip.h"
#include "Rodin/Utility/Product.h"
#include "Rodin/Utility/Repeat.h"
#include "Rodin/Utility/UnwrapReference.h"
#include "Rodin/Utility/HasTypeMember.h"
#include "Rodin/Utility/HasValueMember.h"
#include "Rodin/Utility/Extract.h"
#include "Rodin/Pair.h"
#include "Rodin/Tuple.h"

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

// -----------------------------------------------------------------------------
// Test class for Zip functionality
// -----------------------------------------------------------------------------
class UtilityZipTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(UtilityZipTest, ZipTwoSingleElementTuples)
{
  using T1 = Rodin::Tuple<int>;
  using T2 = Rodin::Tuple<double>;
  using ZippedType = Zip<T1, T2>::Type<std::pair>;
  
  // Should produce Tuple<std::pair<int, double>>
  EXPECT_TRUE((std::is_same_v<ZippedType, Rodin::Tuple<std::pair<int, double>>>));
}

TEST_F(UtilityZipTest, ZipTwoTuples)
{
  using T1 = Rodin::Tuple<int, float>;
  using T2 = Rodin::Tuple<double, char>;
  using ZippedType = Zip<T1, T2>::Type<std::pair>;
  
  // Should produce Tuple<std::pair<int, double>, std::pair<float, char>>
  EXPECT_TRUE((std::is_same_v<
    ZippedType,
    Rodin::Tuple<std::pair<int, double>, std::pair<float, char>>
  >));
}

TEST_F(UtilityZipTest, ZipThreeElementTuples)
{
  using T1 = Rodin::Tuple<int, float, char>;
  using T2 = Rodin::Tuple<double, long, bool>;
  using ZippedType = Zip<T1, T2>::Type<std::pair>;
  
  using ExpectedType = Rodin::Tuple<
    std::pair<int, double>,
    std::pair<float, long>,
    std::pair<char, bool>
  >;
  
  EXPECT_TRUE((std::is_same_v<ZippedType, ExpectedType>));
}

TEST_F(UtilityZipTest, ZipWithRodinPair)
{
  using T1 = Rodin::Tuple<int, float>;
  using T2 = Rodin::Tuple<double, char>;
  using ZippedType = Zip<T1, T2>::Type<Rodin::Pair>;
  
  // Should produce Tuple<Pair<int, double>, Pair<float, char>>
  EXPECT_TRUE((std::is_same_v<
    ZippedType,
    Rodin::Tuple<Rodin::Pair<int, double>, Rodin::Pair<float, char>>
  >));
}

// -----------------------------------------------------------------------------
// Test class for Product functionality
// -----------------------------------------------------------------------------
class UtilityProductTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(UtilityProductTest, ProductSingleElements)
{
  using T1 = Rodin::Tuple<int>;
  using T2 = Rodin::Tuple<double>;
  using ProductType = Product<T1, T2>::Type<std::pair>;
  
  // Should produce Tuple<std::pair<int, double>>
  EXPECT_TRUE((std::is_same_v<ProductType, Rodin::Tuple<std::pair<int, double>>>));
}

TEST_F(UtilityProductTest, ProductOneByTwo)
{
  using T1 = Rodin::Tuple<int>;
  using T2 = Rodin::Tuple<double, char>;
  using ProductType = Product<T1, T2>::Type<std::pair>;
  
  // Should produce Tuple<std::pair<int, double>, std::pair<int, char>>
  using ExpectedType = Rodin::Tuple<std::pair<int, double>, std::pair<int, char>>;
  EXPECT_TRUE((std::is_same_v<ProductType, ExpectedType>));
}

TEST_F(UtilityProductTest, ProductTwoByTwo)
{
  using T1 = Rodin::Tuple<int, float>;
  using T2 = Rodin::Tuple<double, char>;
  using ProductType = Product<T1, T2>::Type<std::pair>;
  
  // Should produce Tuple<
  //   std::pair<int, double>, std::pair<int, char>,
  //   std::pair<float, double>, std::pair<float, char>
  // >
  using ExpectedType = Rodin::Tuple<
    std::pair<int, double>, std::pair<int, char>,
    std::pair<float, double>, std::pair<float, char>
  >;
  EXPECT_TRUE((std::is_same_v<ProductType, ExpectedType>));
}

TEST_F(UtilityProductTest, ProductTwoByThree)
{
  using T1 = Rodin::Tuple<int, float>;
  using T2 = Rodin::Tuple<double, char, bool>;
  using ProductType = Product<T1, T2>::Type<std::pair>;
  
  using ExpectedType = Rodin::Tuple<
    std::pair<int, double>, std::pair<int, char>, std::pair<int, bool>,
    std::pair<float, double>, std::pair<float, char>, std::pair<float, bool>
  >;
  EXPECT_TRUE((std::is_same_v<ProductType, ExpectedType>));
}

TEST_F(UtilityProductTest, ProductThreeByTwo)
{
  using T1 = Rodin::Tuple<int, float, char>;
  using T2 = Rodin::Tuple<double, bool>;
  using ProductType = Product<T1, T2>::Type<std::pair>;
  
  using ExpectedType = Rodin::Tuple<
    std::pair<int, double>, std::pair<int, bool>,
    std::pair<float, double>, std::pair<float, bool>,
    std::pair<char, double>, std::pair<char, bool>
  >;
  EXPECT_TRUE((std::is_same_v<ProductType, ExpectedType>));
}


TEST_F(UtilityProductTest, ProductWithRodinPair)
{
  using T1 = Rodin::Tuple<int, float>;
  using T2 = Rodin::Tuple<double, char>;
  using ProductType = Product<T1, T2>::Type<Rodin::Pair>;
  
  using ExpectedType = Rodin::Tuple<
    Rodin::Pair<int, double>, Rodin::Pair<int, char>,
    Rodin::Pair<float, double>, Rodin::Pair<float, char>
  >;
  EXPECT_TRUE((std::is_same_v<ProductType, ExpectedType>));
}

// -----------------------------------------------------------------------------
// Test class for IsOneOf functionality
// -----------------------------------------------------------------------------
class UtilityRepeatTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(UtilityRepeatTest, RepeatZeroTimes)
{
  using ResultType = Repeat<0, int>::Type;
  EXPECT_TRUE((std::is_same_v<ResultType, Rodin::Tuple<>>));
}

TEST_F(UtilityRepeatTest, RepeatOnce)
{
  using ResultType = Repeat<1, int>::Type;
  EXPECT_TRUE((std::is_same_v<ResultType, Rodin::Tuple<int>>));
}

TEST_F(UtilityRepeatTest, RepeatTwice)
{
  using ResultType = Repeat<2, double>::Type;
  EXPECT_TRUE((std::is_same_v<ResultType, Rodin::Tuple<double, double>>));
}

TEST_F(UtilityRepeatTest, RepeatThreeTimes)
{
  using ResultType = Repeat<3, char>::Type;
  EXPECT_TRUE((std::is_same_v<ResultType, Rodin::Tuple<char, char, char>>));
}

TEST_F(UtilityRepeatTest, RepeatFiveTimes)
{
  using ResultType = Repeat<5, bool>::Type;
  using ExpectedType = Rodin::Tuple<bool, bool, bool, bool, bool>;
  EXPECT_TRUE((std::is_same_v<ResultType, ExpectedType>));
}

TEST_F(UtilityRepeatTest, RepeatWithComplexType)
{
  using ResultType = Repeat<3, std::string>::Type;
  using ExpectedType = Rodin::Tuple<std::string, std::string, std::string>;
  EXPECT_TRUE((std::is_same_v<ResultType, ExpectedType>));
}

// -----------------------------------------------------------------------------
// Test class for UnwrapReference functionality
// -----------------------------------------------------------------------------
class UtilityUnwrapReferenceTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(UtilityUnwrapReferenceTest, NonReferenceWrapper)
{
  using ResultType = UnwrapReference<int>::Type;
  EXPECT_TRUE((std::is_same_v<ResultType, int>));
}

TEST_F(UtilityUnwrapReferenceTest, UnwrapReferenceWrapper)
{
  using ResultType = UnwrapReference<std::reference_wrapper<int>>::Type;
  EXPECT_TRUE((std::is_same_v<ResultType, int&>));
}

TEST_F(UtilityUnwrapReferenceTest, UnwrapConstReferenceWrapper)
{
  using ResultType = UnwrapReference<std::reference_wrapper<const int>>::Type;
  EXPECT_TRUE((std::is_same_v<ResultType, const int&>));
}

TEST_F(UtilityUnwrapReferenceTest, ConstTypeNoWrapper)
{
  using ResultType = UnwrapReference<const double>::Type;
  EXPECT_TRUE((std::is_same_v<ResultType, const double>));
}

TEST_F(UtilityUnwrapReferenceTest, UnwrapRefDecay)
{
  using ResultType = UnwrapRefDecay<std::reference_wrapper<int>>::Type;
  EXPECT_TRUE((std::is_same_v<ResultType, int&>));
}

TEST_F(UtilityUnwrapReferenceTest, UnwrapRefDecayWithConst)
{
  using ResultType = UnwrapRefDecay<const std::reference_wrapper<int>>::Type;
  EXPECT_TRUE((std::is_same_v<ResultType, int&>));
}

TEST_F(UtilityUnwrapReferenceTest, UnwrapRefDecayPlainType)
{
  using ResultType = UnwrapRefDecay<const int&>::Type;
  EXPECT_TRUE((std::is_same_v<ResultType, int>));
}

// -----------------------------------------------------------------------------
// Test class for HasTypeMember functionality
// -----------------------------------------------------------------------------
class UtilityHasTypeMemberTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
    
    struct WithType { using Type = int; };
    struct WithoutType { int value; };
    struct WithNestedType { using Type = std::string; };
};

TEST_F(UtilityHasTypeMemberTest, TypeWithTypeMember)
{
  EXPECT_TRUE((HasTypeMember<WithType>::Value));
}

TEST_F(UtilityHasTypeMemberTest, TypeWithoutTypeMember)
{
  EXPECT_FALSE((HasTypeMember<WithoutType>::Value));
}

TEST_F(UtilityHasTypeMemberTest, PrimitiveType)
{
  EXPECT_FALSE((HasTypeMember<int>::Value));
  EXPECT_FALSE((HasTypeMember<double>::Value));
}

TEST_F(UtilityHasTypeMemberTest, TypeWithNestedType)
{
  EXPECT_TRUE((HasTypeMember<WithNestedType>::Value));
}

TEST_F(UtilityHasTypeMemberTest, StdTypes)
{
  EXPECT_FALSE((HasTypeMember<std::string>::Value));
  EXPECT_FALSE((HasTypeMember<std::vector<int>>::Value));
}

// -----------------------------------------------------------------------------
// Test class for HasValueMember functionality
// -----------------------------------------------------------------------------
class UtilityHasValueMemberTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
    
    struct WithValue { using Value = int; };
    struct WithoutValue { int data; };
    struct WithStaticValue { static constexpr int Value = 42; };
};

TEST_F(UtilityHasValueMemberTest, TypeWithValueMember)
{
  EXPECT_TRUE((HasValueMember<WithValue>::Value));
}

TEST_F(UtilityHasValueMemberTest, TypeWithoutValueMember)
{
  EXPECT_FALSE((HasValueMember<WithoutValue>::Value));
}

TEST_F(UtilityHasValueMemberTest, PrimitiveType)
{
  EXPECT_FALSE((HasValueMember<int>::Value));
  EXPECT_FALSE((HasValueMember<double>::Value));
}

TEST_F(UtilityHasValueMemberTest, TypeWithStaticValue)
{
  // HasValueMember checks for a type member, not a static value
  // So this should be FALSE since WithStaticValue has a static const, not a type
  EXPECT_FALSE((HasValueMember<WithStaticValue>::Value));
}

TEST_F(UtilityHasValueMemberTest, StdTypes)
{
  EXPECT_FALSE((HasValueMember<std::string>::Value));
  EXPECT_FALSE((HasValueMember<std::vector<int>>::Value));
}

// -----------------------------------------------------------------------------
// Test class for Extract functionality
// -----------------------------------------------------------------------------
class UtilityExtractTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
    
    template<typename T>
    struct AddConst { using Type = const T; };
    
    template<typename T>
    struct AddPointer { using Type = T*; };
    
    template<typename T>
    struct Identity { using Type = T; };
};

TEST_F(UtilityExtractTest, ExtractWithIdentity)
{
  using OriginalTuple = Rodin::Tuple<int, double, char>;
  using ResultType = Extract<OriginalTuple>::Type<Identity>;
  EXPECT_TRUE((std::is_same_v<ResultType, Rodin::Tuple<int, double, char>>));
}

TEST_F(UtilityExtractTest, ExtractWithAddConst)
{
  using OriginalTuple = Rodin::Tuple<int, double, char>;
  using ResultType = Extract<OriginalTuple>::Type<AddConst>;
  using ExpectedType = Rodin::Tuple<const int, const double, const char>;
  EXPECT_TRUE((std::is_same_v<ResultType, ExpectedType>));
}

TEST_F(UtilityExtractTest, ExtractWithAddPointer)
{
  using OriginalTuple = Rodin::Tuple<int, double, char>;
  using ResultType = Extract<OriginalTuple>::Type<AddPointer>;
  using ExpectedType = Rodin::Tuple<int*, double*, char*>;
  EXPECT_TRUE((std::is_same_v<ResultType, ExpectedType>));
}

TEST_F(UtilityExtractTest, ExtractSingleElement)
{
  using OriginalTuple = Rodin::Tuple<int>;
  using ResultType = Extract<OriginalTuple>::Type<AddConst>;
  EXPECT_TRUE((std::is_same_v<ResultType, Rodin::Tuple<const int>>));
}

TEST_F(UtilityExtractTest, ExtractMultipleElements)
{
  using OriginalTuple = Rodin::Tuple<int, float, double, char, bool>;
  using ResultType = Extract<OriginalTuple>::Type<AddPointer>;
  using ExpectedType = Rodin::Tuple<int*, float*, double*, char*, bool*>;
  EXPECT_TRUE((std::is_same_v<ResultType, ExpectedType>));
}
