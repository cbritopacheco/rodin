#include <gtest/gtest.h>

#include <Rodin/Types.h>

using namespace Rodin;

namespace Rodin::Tests::Unit
{
  //=== Basic Type Tests =====================================================

  TEST(Rodin_Types, BasicTypes)
  {
    // Test that basic types are defined and have expected properties
    Integer i = 42;
    Boolean b = true;
    Float f = 3.14f;
    Double d = 2.71828;
    Index idx = 100;
    Real r = 1.414;
    
    EXPECT_EQ(i, 42);
    EXPECT_TRUE(b);
    EXPECT_FLOAT_EQ(f, 3.14f);
    EXPECT_DOUBLE_EQ(d, 2.71828);
    EXPECT_EQ(idx, 100);
    EXPECT_DOUBLE_EQ(r, 1.414);
  }

  TEST(Rodin_Types, TypeSizes)
  {
    // Verify expected sizes for the types
    EXPECT_EQ(sizeof(Integer), sizeof(int));
    EXPECT_EQ(sizeof(Boolean), sizeof(bool));
    EXPECT_EQ(sizeof(Float), sizeof(float));
    EXPECT_EQ(sizeof(Double), sizeof(double));
    EXPECT_EQ(sizeof(Index), sizeof(std::size_t));
    EXPECT_EQ(sizeof(Real), sizeof(Double));
  }

  TEST(Rodin_Types, ComplexType)
  {
    Complex c1(1.0, 2.0);
    Complex c2 = Complex(3.0, 4.0);
    
    EXPECT_DOUBLE_EQ(c1.real(), 1.0);
    EXPECT_DOUBLE_EQ(c1.imag(), 2.0);
    EXPECT_DOUBLE_EQ(c2.real(), 3.0);
    EXPECT_DOUBLE_EQ(c2.imag(), 4.0);
    
    // Test complex arithmetic
    Complex sum = c1 + c2;
    EXPECT_DOUBLE_EQ(sum.real(), 4.0);
    EXPECT_DOUBLE_EQ(sum.imag(), 6.0);
  }

  //=== Container Type Tests =================================================

  TEST(Rodin_Types, ListType)
  {
    List<Integer> intList;
    intList.push_back(1);
    intList.push_back(2);
    intList.push_back(3);
    
    EXPECT_EQ(intList.size(), 3);
    EXPECT_EQ(intList.front(), 1);
    EXPECT_EQ(intList.back(), 3);
  }

  TEST(Rodin_Types, DequeType)
  {
    Deque<Double> doubleDeque;
    doubleDeque.push_back(1.1);
    doubleDeque.push_front(2.2);
    doubleDeque.push_back(3.3);
    
    EXPECT_EQ(doubleDeque.size(), 3);
    EXPECT_DOUBLE_EQ(doubleDeque[0], 2.2);
    EXPECT_DOUBLE_EQ(doubleDeque[1], 1.1);
    EXPECT_DOUBLE_EQ(doubleDeque[2], 3.3);
  }

  TEST(Rodin_Types, StackType)
  {
    Stack<Integer> intStack;
    intStack.push(10);
    intStack.push(20);
    intStack.push(30);
    
    EXPECT_EQ(intStack.size(), 3);
    EXPECT_EQ(intStack.top(), 30);
    
    intStack.pop();
    EXPECT_EQ(intStack.top(), 20);
    EXPECT_EQ(intStack.size(), 2);
  }

  TEST(Rodin_Types, FlatSetType)
  {
    FlatSet<Integer> intSet;
    intSet.insert(5);
    intSet.insert(1);
    intSet.insert(3);
    intSet.insert(1);  // Duplicate
    
    EXPECT_EQ(intSet.size(), 3);  // Should only have 3 unique elements
    
    // FlatSet maintains sorted order
    auto it = intSet.begin();
    EXPECT_EQ(*it++, 1);
    EXPECT_EQ(*it++, 3);
    EXPECT_EQ(*it++, 5);
  }

  TEST(Rodin_Types, UnorderedSetType)
  {
    UnorderedSet<Integer> intSet;
    intSet.insert(10);
    intSet.insert(20);
    intSet.insert(30);
    intSet.insert(20);  // Duplicate
    
    EXPECT_EQ(intSet.size(), 3);
    EXPECT_TRUE(intSet.find(10) != intSet.end());
    EXPECT_TRUE(intSet.find(20) != intSet.end());
    EXPECT_TRUE(intSet.find(30) != intSet.end());
    EXPECT_TRUE(intSet.find(40) == intSet.end());
  }

  TEST(Rodin_Types, MapType)
  {
    Map<Integer, Double> intDoubleMap;
    intDoubleMap[1] = 1.1;
    intDoubleMap[2] = 2.2;
    intDoubleMap[3] = 3.3;
    
    EXPECT_EQ(intDoubleMap.size(), 3);
    EXPECT_DOUBLE_EQ(intDoubleMap[1], 1.1);
    EXPECT_DOUBLE_EQ(intDoubleMap[2], 2.2);
    EXPECT_DOUBLE_EQ(intDoubleMap[3], 3.3);
  }

  TEST(Rodin_Types, UnorderedMapType)
  {
    UnorderedMap<Integer, Double> intDoubleMap;
    intDoubleMap[10] = 10.1;
    intDoubleMap[20] = 20.2;
    intDoubleMap[30] = 30.3;
    
    EXPECT_EQ(intDoubleMap.size(), 3);
    EXPECT_DOUBLE_EQ(intDoubleMap[10], 10.1);
    EXPECT_DOUBLE_EQ(intDoubleMap[20], 20.2);
    EXPECT_DOUBLE_EQ(intDoubleMap[30], 30.3);
  }

  TEST(Rodin_Types, FlatMapType)
  {
    FlatMap<Integer, Double> intDoubleMap;
    intDoubleMap[5] = 5.5;
    intDoubleMap[1] = 1.1;
    intDoubleMap[3] = 3.3;
    
    EXPECT_EQ(intDoubleMap.size(), 3);
    EXPECT_DOUBLE_EQ(intDoubleMap[1], 1.1);
    EXPECT_DOUBLE_EQ(intDoubleMap[3], 3.3);
    EXPECT_DOUBLE_EQ(intDoubleMap[5], 5.5);
    
    // FlatMap maintains sorted order
    auto it = intDoubleMap.begin();
    EXPECT_EQ(it++->first, 1);
    EXPECT_EQ(it++->first, 3);
    EXPECT_EQ(it++->first, 5);
  }

  //=== Index-specific Types =================================================

  TEST(Rodin_Types, IndexSetType)
  {
    IndexSet indices;
    indices.insert(100);
    indices.insert(50);
    indices.insert(200);
    indices.insert(50);  // Duplicate
    
    EXPECT_EQ(indices.size(), 3);
    
    // Should be sorted
    auto it = indices.begin();
    EXPECT_EQ(*it++, 50);
    EXPECT_EQ(*it++, 100);
    EXPECT_EQ(*it++, 200);
  }

  TEST(Rodin_Types, IndexMapType)
  {
    IndexMap<Double> indexDoubleMap;
    indexDoubleMap[100] = 10.0;
    indexDoubleMap[50] = 5.0;
    indexDoubleMap[200] = 20.0;
    
    EXPECT_EQ(indexDoubleMap.size(), 3);
    EXPECT_DOUBLE_EQ(indexDoubleMap[50], 5.0);
    EXPECT_DOUBLE_EQ(indexDoubleMap[100], 10.0);
    EXPECT_DOUBLE_EQ(indexDoubleMap[200], 20.0);
  }

  //=== Utility Types ========================================================

  TEST(Rodin_Types, BitSetTypes)
  {
    BitSet<8> bits8;
    bits8.set(0);
    bits8.set(3);
    bits8.set(7);
    
    EXPECT_TRUE(bits8[0]);
    EXPECT_FALSE(bits8[1]);
    EXPECT_FALSE(bits8[2]);
    EXPECT_TRUE(bits8[3]);
    EXPECT_FALSE(bits8[4]);
    EXPECT_FALSE(bits8[5]);
    EXPECT_FALSE(bits8[6]);
    EXPECT_TRUE(bits8[7]);
    
    BitSet2 bits2;
    bits2.set(0);
    EXPECT_TRUE(bits2[0]);
    EXPECT_FALSE(bits2[1]);
  }

  TEST(Rodin_Types, OptionalType)
  {
    Optional<Integer> optInt;
    EXPECT_FALSE(optInt.has_value());
    
    optInt = 42;
    EXPECT_TRUE(optInt.has_value());
    EXPECT_EQ(optInt.value(), 42);
    EXPECT_EQ(*optInt, 42);
    
    Optional<Double> optDouble = 3.14;
    EXPECT_TRUE(optDouble.has_value());
    EXPECT_DOUBLE_EQ(*optDouble, 3.14);
  }

  //=== User-defined Literal Tests ===========================================

#if __cpp_size_t_suffix < 202011L
  TEST(Rodin_Types, UserDefinedLiteral_UZ)
  {
    // Test the _UZ suffix for size_t values
    auto value1 = 42_UZ;
    auto value2 = 0_UZ;
    auto value3 = 1000000_UZ;
    
    EXPECT_EQ(value1, static_cast<std::size_t>(42));
    EXPECT_EQ(value2, static_cast<std::size_t>(0));
    EXPECT_EQ(value3, static_cast<std::size_t>(1000000));
    
    // Verify the type
    static_assert(std::is_same_v<decltype(value1), std::size_t>);
    static_assert(std::is_same_v<decltype(value2), std::size_t>);
    static_assert(std::is_same_v<decltype(value3), std::size_t>);
  }
#endif

  //=== Template Instantiation Tests ========================================

  TEST(Rodin_Types, TemplateInstantiation)
  {
    // Test that template types can be instantiated with different types
    List<Boolean> boolList;
    boolList.push_back(true);
    boolList.push_back(false);
    EXPECT_EQ(boolList.size(), 2);
    
    FlatSet<Real> realSet;
    realSet.insert(1.1);
    realSet.insert(2.2);
    EXPECT_EQ(realSet.size(), 2);
    
    Map<Index, Complex> complexMap;
    complexMap[0] = Complex(1.0, 2.0);
    complexMap[1] = Complex(3.0, 4.0);
    EXPECT_EQ(complexMap.size(), 2);
  }
}