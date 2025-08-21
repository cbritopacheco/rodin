/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/FormLanguage/List.h"
#include "Rodin/FormLanguage/Base.h"
#include "Rodin/Types.h"

using namespace Rodin;
using namespace Rodin::FormLanguage;

// Simple test class that inherits from Base
class TestBaseObject : public Base
{
  public:
    TestBaseObject(int value) : m_value(value) {}
    TestBaseObject* copy() const noexcept override 
    { 
      return new TestBaseObject(m_value); 
    }
    int getValue() const { return m_value; }
    void setValue(int value) { m_value = value; }
  private:
    int m_value;
};

class FormLanguageListTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// Test basic List construction and functionality
TEST_F(FormLanguageListTest, BasicConstruction)
{
  // Test empty list construction
  List<TestBaseObject> list;
  EXPECT_EQ(list.size(), 0);
  EXPECT_TRUE(list.empty());
}

// Test adding elements to list
TEST_F(FormLanguageListTest, AddElements)
{
  List<TestBaseObject> list;
  
  // Add some test objects
  list.emplace_back(std::make_unique<TestBaseObject>(1));
  list.emplace_back(std::make_unique<TestBaseObject>(2));
  list.emplace_back(std::make_unique<TestBaseObject>(3));
  
  EXPECT_EQ(list.size(), 3);
  EXPECT_FALSE(list.empty());
}

// Test element access
TEST_F(FormLanguageListTest, ElementAccess)
{
  List<TestBaseObject> list;
  
  // Add test objects
  list.emplace_back(std::make_unique<TestBaseObject>(10));
  list.emplace_back(std::make_unique<TestBaseObject>(20));
  list.emplace_back(std::make_unique<TestBaseObject>(30));
  
  // Test operator[] access
  EXPECT_EQ(list[0].getValue(), 10);
  EXPECT_EQ(list[1].getValue(), 20);
  EXPECT_EQ(list[2].getValue(), 30);
  
  // Test at() access
  EXPECT_EQ(list.at(0).getValue(), 10);
  EXPECT_EQ(list.at(1).getValue(), 20);
  EXPECT_EQ(list.at(2).getValue(), 30);
  
  // Test front() and back()
  EXPECT_EQ(list.front().getValue(), 10);
  EXPECT_EQ(list.back().getValue(), 30);
}

// Test const element access
TEST_F(FormLanguageListTest, ConstElementAccess)
{
  List<TestBaseObject> list;
  list.emplace_back(std::make_unique<TestBaseObject>(100));
  list.emplace_back(std::make_unique<TestBaseObject>(200));
  
  const List<TestBaseObject>& constList = list;
  
  // Test const operator[] access
  EXPECT_EQ(constList[0].getValue(), 100);
  EXPECT_EQ(constList[1].getValue(), 200);
  
  // Test const at() access
  EXPECT_EQ(constList.at(0).getValue(), 100);
  EXPECT_EQ(constList.at(1).getValue(), 200);
  
  // Test const front() and back()
  EXPECT_EQ(constList.front().getValue(), 100);
  EXPECT_EQ(constList.back().getValue(), 200);
}

// Test iterators
TEST_F(FormLanguageListTest, Iterators)
{
  List<TestBaseObject> list;
  list.emplace_back(std::make_unique<TestBaseObject>(1));
  list.emplace_back(std::make_unique<TestBaseObject>(2));
  list.emplace_back(std::make_unique<TestBaseObject>(3));
  
  // Test begin() and end()
  auto it = list.begin();
  EXPECT_EQ((*it).getValue(), 1);
  ++it;
  EXPECT_EQ((*it).getValue(), 2);
  ++it;
  EXPECT_EQ((*it).getValue(), 3);
  ++it;
  EXPECT_EQ(it, list.end());
  
  // Test range-based for loop
  std::vector<int> values;
  for (auto& obj : list)
  {
    values.push_back(obj.getValue());
  }
  EXPECT_EQ(values, std::vector<int>({1, 2, 3}));
}

// Test const iterators
TEST_F(FormLanguageListTest, ConstIterators)
{
  List<TestBaseObject> list;
  list.emplace_back(std::make_unique<TestBaseObject>(10));
  list.emplace_back(std::make_unique<TestBaseObject>(20));
  
  const List<TestBaseObject>& constList = list;
  
  // Test cbegin() and cend()
  auto it = constList.begin();
  EXPECT_EQ((*it).getValue(), 10);
  ++it;
  EXPECT_EQ((*it).getValue(), 20);
  ++it;
  EXPECT_EQ(it, constList.end());
  
  // Test const range-based for loop
  std::vector<int> values;
  for (const auto& obj : constList)
  {
    values.push_back(obj.getValue());
  }
  EXPECT_EQ(values, std::vector<int>({10, 20}));
}

// Test iterator operations
TEST_F(FormLanguageListTest, IteratorOperations)
{
  List<TestBaseObject> list;
  list.emplace_back(std::make_unique<TestBaseObject>(5));
  list.emplace_back(std::make_unique<TestBaseObject>(15));
  list.emplace_back(std::make_unique<TestBaseObject>(25));
  
  // Test pre-increment
  auto it = list.begin();
  EXPECT_EQ((*it).getValue(), 5);
  ++it;
  EXPECT_EQ((*it).getValue(), 15);
  
  // Test post-increment
  auto it2 = list.begin();
  auto oldIt = it2++;
  EXPECT_EQ((*oldIt).getValue(), 5);
  EXPECT_EQ((*it2).getValue(), 15);
  
  // Test equality/inequality
  auto it3 = list.begin();
  auto it4 = list.begin();
  EXPECT_TRUE(it3 == it4);
  EXPECT_FALSE(it3 != it4);
  
  ++it4;
  EXPECT_FALSE(it3 == it4);
  EXPECT_TRUE(it3 != it4);
}

// Test list modification
TEST_F(FormLanguageListTest, ListModification)
{
  List<TestBaseObject> list;
  list.emplace_back(std::make_unique<TestBaseObject>(1));
  list.emplace_back(std::make_unique<TestBaseObject>(2));
  
  // Modify elements through iterators
  for (auto& obj : list)
  {
    obj.setValue(obj.getValue() * 10);
  }
  
  EXPECT_EQ(list[0].getValue(), 10);
  EXPECT_EQ(list[1].getValue(), 20);
  
  // Modify elements through direct access
  list[0].setValue(100);
  EXPECT_EQ(list[0].getValue(), 100);
}

// Test list copying (since it inherits from Base)
TEST_F(FormLanguageListTest, ListCopying)
{
  List<TestBaseObject> list;
  list.emplace_back(std::make_unique<TestBaseObject>(42));
  list.emplace_back(std::make_unique<TestBaseObject>(84));
  
  // Test copy through Base interface
  std::unique_ptr<Base> copyBase(list.copy());
  auto* copyList = dynamic_cast<List<TestBaseObject>*>(copyBase.get());
  
  ASSERT_NE(copyList, nullptr);
  EXPECT_EQ(copyList->size(), 2);
  EXPECT_EQ((*copyList)[0].getValue(), 42);
  EXPECT_EQ((*copyList)[1].getValue(), 84);
  
  // Verify it's a deep copy
  list[0].setValue(999);
  EXPECT_EQ((*copyList)[0].getValue(), 42);  // Should remain unchanged
}

// Test empty list operations
TEST_F(FormLanguageListTest, EmptyListOperations)
{
  List<TestBaseObject> emptyList;
  
  // Test size and empty
  EXPECT_EQ(emptyList.size(), 0);
  EXPECT_TRUE(emptyList.empty());
  
  // Test iterators on empty list
  EXPECT_EQ(emptyList.begin(), emptyList.end());
  
  // Test range-based for loop on empty list (should not execute)
  int count = 0;
  for (auto& obj : emptyList)
  {
    count++;
    // This should never execute
    (void)obj;  // Suppress unused variable warning
  }
  EXPECT_EQ(count, 0);
}

// Test unique_ptr management
TEST_F(FormLanguageListTest, UniquePointerManagement)
{
  List<TestBaseObject> list;
  
  // Test that unique_ptr ownership is transferred correctly
  auto ptr1 = std::make_unique<TestBaseObject>(123);
  auto* rawPtr1 = ptr1.get();
  
  list.emplace_back(std::move(ptr1));
  EXPECT_EQ(ptr1.get(), nullptr);  // Ownership transferred
  EXPECT_EQ(list[0].getValue(), 123);
  EXPECT_EQ(&list[0], rawPtr1);  // Same object
  
  // Test that the list owns the objects
  EXPECT_EQ(list.size(), 1);
}

// Test large list performance (basic)
TEST_F(FormLanguageListTest, LargeListBasic)
{
  List<TestBaseObject> list;
  
  // Add many elements
  const int numElements = 1000;
  for (int i = 0; i < numElements; ++i)
  {
    list.emplace_back(std::make_unique<TestBaseObject>(i));
  }
  
  EXPECT_EQ(list.size(), numElements);
  
  // Verify all elements
  for (int i = 0; i < numElements; ++i)
  {
    EXPECT_EQ(list[i].getValue(), i);
  }
  
  // Test iteration performance
  int sum = 0;
  for (const auto& obj : list)
  {
    sum += obj.getValue();
  }
  
  // Expected sum = 0 + 1 + 2 + ... + (n-1) = n*(n-1)/2
  int expectedSum = numElements * (numElements - 1) / 2;
  EXPECT_EQ(sum, expectedSum);
}