#include <gtest/gtest.h>

#include <Rodin/Cast.h>

using namespace Rodin;

namespace Rodin::Tests::Unit
{
  // Test classes for casting functionality
  class BaseClass
  {
    public:
      virtual ~BaseClass() = default;
      virtual int getValue() const { return 42; }
  };

  class DerivedClass : public BaseClass
  {
    public:
      int getValue() const override { return 100; }
      int getSpecialValue() const { return 200; }
  };

  class UnrelatedClass
  {
    public:
      int getValue() const { return 999; }
  };
}

// Template specializations for Cast - these demonstrate how Cast might be specialized
template <>
class Rodin::Cast<int, double>
{
  public:
    double operator()(const int& value) const
    {
      return static_cast<double>(value);
    }
};

template <>
class Rodin::Cast<double, int>
{
  public:
    int operator()(const double& value) const
    {
      return static_cast<int>(value);
    }
};

template <>
class Rodin::Cast<Rodin::Tests::Unit::BaseClass*, Rodin::Tests::Unit::DerivedClass*>
{
  public:
    Rodin::Tests::Unit::DerivedClass* operator()(Rodin::Tests::Unit::BaseClass* ptr) const
    {
      return dynamic_cast<Rodin::Tests::Unit::DerivedClass*>(ptr);
    }
};

namespace Rodin::Tests::Unit
{



  //=== Template Declaration Tests ===========================================

  TEST(Rodin_Cast, TemplateDeclaration)
  {
    // Test that the Cast template class is declared and can be used
    // This is primarily a compilation test
    
    // The basic template should be declared (though not implemented)
    static_assert(std::is_class_v<Cast<int, double>>);
    static_assert(std::is_class_v<Cast<double, int>>);
    static_assert(std::is_class_v<Cast<BaseClass*, DerivedClass*>>);
  }

  //=== Specialization Tests =================================================

  TEST(Rodin_Cast, IntToDoubleCast)
  {
    Cast<int, double> caster;
    
    int intValue = 42;
    double result = caster(intValue);
    
    EXPECT_DOUBLE_EQ(result, 42.0);
  }

  TEST(Rodin_Cast, DoubleToIntCast)
  {
    Cast<double, int> caster;
    
    double doubleValue = 3.14;
    int result = caster(doubleValue);
    
    EXPECT_EQ(result, 3);  // Truncation expected
  }

  TEST(Rodin_Cast, PointerCastSuccess)
  {
    DerivedClass derived;
    BaseClass* basePtr = &derived;
    
    Cast<BaseClass*, DerivedClass*> caster;
    DerivedClass* result = caster(basePtr);
    
    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->getValue(), 100);
    EXPECT_EQ(result->getSpecialValue(), 200);
  }

  TEST(Rodin_Cast, PointerCastFailure)
  {
    BaseClass base;
    BaseClass* basePtr = &base;
    
    Cast<BaseClass*, DerivedClass*> caster;
    DerivedClass* result = caster(basePtr);
    
    EXPECT_EQ(result, nullptr);  // dynamic_cast should fail
  }

  //=== Edge Cases ============================================================

  TEST(Rodin_Cast, ZeroValues)
  {
    Cast<int, double> intToDouble;
    Cast<double, int> doubleToInt;
    
    EXPECT_DOUBLE_EQ(intToDouble(0), 0.0);
    EXPECT_EQ(doubleToInt(0.0), 0);
  }

  TEST(Rodin_Cast, NegativeValues)
  {
    Cast<int, double> intToDouble;
    Cast<double, int> doubleToInt;
    
    EXPECT_DOUBLE_EQ(intToDouble(-42), -42.0);
    EXPECT_EQ(doubleToInt(-3.14), -3);
  }

  TEST(Rodin_Cast, LargeValues)
  {
    Cast<int, double> intToDouble;
    
    int largeInt = 1000000;
    double result = intToDouble(largeInt);
    
    EXPECT_DOUBLE_EQ(result, 1000000.0);
  }

  TEST(Rodin_Cast, NullPointer)
  {
    Cast<BaseClass*, DerivedClass*> caster;
    
    BaseClass* nullPtr = nullptr;
    DerivedClass* result = caster(nullPtr);
    
    EXPECT_EQ(result, nullptr);
  }

  //=== Functional Tests =====================================================

  TEST(Rodin_Cast, FunctorBehavior)
  {
    // Test that Cast objects can be used as functors
    Cast<int, double> caster;
    
    // Should be callable with operator()
    auto result = caster(123);
    EXPECT_DOUBLE_EQ(result, 123.0);
    
    // Should be copyable
    auto casterCopy = caster;
    auto result2 = casterCopy(456);
    EXPECT_DOUBLE_EQ(result2, 456.0);
  }

  TEST(Rodin_Cast, ConstCorrectness)
  {
    const Cast<int, double> constCaster;
    
    int value = 789;
    double result = constCaster(value);
    
    EXPECT_DOUBLE_EQ(result, 789.0);
  }

  //=== Type Safety Tests ====================================================

  TEST(Rodin_Cast, TypeSafety)
  {
    // These should compile correctly
    Cast<int, double> intToDouble;
    Cast<double, int> doubleToInt;
    Cast<BaseClass*, DerivedClass*> ptrCast;
    
    // Verify types are distinct
    static_assert(!std::is_same_v<decltype(intToDouble), decltype(doubleToInt)>);
    static_assert(!std::is_same_v<decltype(intToDouble), decltype(ptrCast)>);
  }

  //=== Documentation Example Test ===========================================

  TEST(Rodin_Cast, DocumentationExample)
  {
    // Test the example from the documentation comment
    Cast<int, double> converter;
    int sourceObject = 42;
    double result = converter(sourceObject);
    
    EXPECT_DOUBLE_EQ(result, 42.0);
  }
}