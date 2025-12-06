/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/FormLanguage/List.h"
#include "Rodin/FormLanguage/Base.h"

// Simple test element class that can be stored in FormLanguage::List
class TestElement : public Rodin::FormLanguage::Base
{
  public:
    TestElement() : m_value(0) {}

    explicit TestElement(int value) : m_value(value) {}

    TestElement(const TestElement& other)
      : Base(other), m_value(other.m_value)
    {}

    TestElement(TestElement&& other)
      : Base(std::move(other)), m_value(other.m_value)
    {}

    TestElement& operator=(const TestElement& other)
    {
      if (this != &other)
      {
        m_value = other.m_value;
      }
      return *this;
    }

    TestElement& operator=(TestElement&& other)
    {
      m_value = other.m_value;
      return *this;
    }

    int getValue() const { return m_value; }
    void setValue(int value) { m_value = value; }

    virtual TestElement* copy() const noexcept override
    {
      return new TestElement(*this);
    }

  private:
    int m_value;
};

namespace Rodin::Tests::Unit
{
  //=== Construction Tests ===================================================

  TEST(FormLanguage_List, DefaultConstruction)
  {
    Rodin::FormLanguage::List<TestElement> list;
    EXPECT_TRUE(list.empty());
    EXPECT_EQ(list.size(), 0);
  }

  TEST(FormLanguage_List, CopyConstruction)
  {
    Rodin::FormLanguage::List<TestElement> list1;
    TestElement elem1(42);
    TestElement elem2(99);

    list1.add(elem1);
    list1.add(elem2);

    Rodin::FormLanguage::List<TestElement> list2(list1);

    EXPECT_EQ(list2.size(), 2);
    EXPECT_EQ(list2.at(0).getValue(), 42);
    EXPECT_EQ(list2.at(1).getValue(), 99);

    // Verify deep copy - modifying list2 should not affect list1
    list2.at(0).setValue(100);
    EXPECT_EQ(list1.at(0).getValue(), 42);
  }

  TEST(FormLanguage_List, MoveConstruction)
  {
    Rodin::FormLanguage::List<TestElement> list1;
    TestElement elem(42);
    list1.add(elem);

    Rodin::FormLanguage::List<TestElement> list2(std::move(list1));

    EXPECT_EQ(list2.size(), 1);
    EXPECT_EQ(list2.at(0).getValue(), 42);
  }

  //=== Assignment Tests =====================================================

  TEST(FormLanguage_List, CopyAssignment)
  {
    Rodin::FormLanguage::List<TestElement> list1;
    TestElement elem1(42);
    TestElement elem2(99);

    list1.add(elem1);
    list1.add(elem2);

    Rodin::FormLanguage::List<TestElement> list2;
    list2 = list1;

    EXPECT_EQ(list2.size(), 2);
    EXPECT_EQ(list2.at(0).getValue(), 42);
    EXPECT_EQ(list2.at(1).getValue(), 99);
  }

  TEST(FormLanguage_List, MoveAssignment)
  {
    Rodin::FormLanguage::List<TestElement> list1;
    TestElement elem(42);
    list1.add(elem);

    Rodin::FormLanguage::List<TestElement> list2;
    list2 = std::move(list1);

    EXPECT_EQ(list2.size(), 1);
    EXPECT_EQ(list2.at(0).getValue(), 42);
  }


  //=== Add and Access Tests =================================================

  TEST(FormLanguage_List, AddSingleElement)
  {
    Rodin::FormLanguage::List<TestElement> list;
    TestElement elem(42);

    list.add(elem);

    EXPECT_FALSE(list.empty());
    EXPECT_EQ(list.size(), 1);
    EXPECT_EQ(list.at(0).getValue(), 42);
  }

  TEST(FormLanguage_List, AddMultipleElements)
  {
    Rodin::FormLanguage::List<TestElement> list;
    TestElement elem1(10);
    TestElement elem2(20);
    TestElement elem3(30);

    list.add(elem1).add(elem2).add(elem3);

    EXPECT_EQ(list.size(), 3);
    EXPECT_EQ(list.at(0).getValue(), 10);
    EXPECT_EQ(list.at(1).getValue(), 20);
    EXPECT_EQ(list.at(2).getValue(), 30);
  }

  TEST(FormLanguage_List, AddAnotherList)
  {
    Rodin::FormLanguage::List<TestElement> list1;
    TestElement elem1(10);
    TestElement elem2(20);
    list1.add(elem1).add(elem2);

    Rodin::FormLanguage::List<TestElement> list2;
    TestElement elem3(30);
    list2.add(elem3);

    list2.add(list1);

    EXPECT_EQ(list2.size(), 3);
    EXPECT_EQ(list2.at(0).getValue(), 30);
    EXPECT_EQ(list2.at(1).getValue(), 10);
    EXPECT_EQ(list2.at(2).getValue(), 20);
  }

  TEST(FormLanguage_List, AtAccessConst)
  {
    Rodin::FormLanguage::List<TestElement> list;
    TestElement elem(42);
    list.add(elem);

    const Rodin::FormLanguage::List<TestElement>& constList = list;
    EXPECT_EQ(constList.at(0).getValue(), 42);
  }

  //=== Clear Test ===========================================================

  TEST(FormLanguage_List, Clear)
  {
    Rodin::FormLanguage::List<TestElement> list;
    TestElement elem1(10);
    TestElement elem2(20);
    TestElement elem3(30);

    list.add(elem1).add(elem2).add(elem3);
    EXPECT_EQ(list.size(), 3);

    list.clear();

    EXPECT_TRUE(list.empty());
    EXPECT_EQ(list.size(), 0);
  }

  //=== Iterator Tests =======================================================

  TEST(FormLanguage_List, IteratorBeginEnd)
  {
    Rodin::FormLanguage::List<TestElement> list;
    TestElement elem1(10);
    TestElement elem2(20);
    TestElement elem3(30);

    list.add(elem1).add(elem2).add(elem3);

    auto it = list.begin();
    EXPECT_EQ((*it).getValue(), 10);

    ++it;
    EXPECT_EQ((*it).getValue(), 20);

    ++it;
    EXPECT_EQ((*it).getValue(), 30);

    ++it;
    EXPECT_EQ(it, list.end());
  }

  TEST(FormLanguage_List, IteratorPostIncrement)
  {
    Rodin::FormLanguage::List<TestElement> list;
    TestElement elem1(10);
    TestElement elem2(20);

    list.add(elem1).add(elem2);

    auto it = list.begin();
    auto it2 = it++;

    EXPECT_EQ((*it2).getValue(), 10);
    EXPECT_EQ((*it).getValue(), 20);
  }

  TEST(FormLanguage_List, IteratorRangeBasedFor)
  {
    Rodin::FormLanguage::List<TestElement> list;
    TestElement elem1(10);
    TestElement elem2(20);
    TestElement elem3(30);

    list.add(elem1).add(elem2).add(elem3);

    int sum = 0;
    for (auto& elem : list)
    {
      sum += elem.getValue();
    }

    EXPECT_EQ(sum, 60);
  }

  TEST(FormLanguage_List, ConstIteratorBeginEnd)
  {
    Rodin::FormLanguage::List<TestElement> list;
    TestElement elem1(10);
    TestElement elem2(20);

    list.add(elem1).add(elem2);

    const Rodin::FormLanguage::List<TestElement>& constList = list;

    auto it = constList.begin();
    EXPECT_EQ((*it).getValue(), 10);

    ++it;
    EXPECT_EQ((*it).getValue(), 20);

    ++it;
    EXPECT_EQ(it, constList.end());
  }

  TEST(FormLanguage_List, ConstIteratorCBeginCEnd)
  {
    Rodin::FormLanguage::List<TestElement> list;
    TestElement elem1(10);
    TestElement elem2(20);

    list.add(elem1).add(elem2);

    auto it = list.cbegin();
    EXPECT_EQ((*it).getValue(), 10);

    ++it;
    EXPECT_EQ((*it).getValue(), 20);

    ++it;
    EXPECT_EQ(it, list.cend());
  }

  TEST(FormLanguage_List, ConstIteratorRangeBasedFor)
  {
    Rodin::FormLanguage::List<TestElement> list;
    TestElement elem1(10);
    TestElement elem2(20);
    TestElement elem3(30);

    list.add(elem1).add(elem2).add(elem3);

    const Rodin::FormLanguage::List<TestElement>& constList = list;

    int sum = 0;
    for (const auto& elem : constList)
    {
      sum += elem.getValue();
    }

    EXPECT_EQ(sum, 60);
  }

  //=== Copy Method Test =====================================================

  TEST(FormLanguage_List, CopyMethod)
  {
    Rodin::FormLanguage::List<TestElement> list;
    TestElement elem1(42);
    TestElement elem2(99);

    list.add(elem1).add(elem2);

    std::unique_ptr<Rodin::FormLanguage::List<TestElement>> copiedList(list.copy());

    EXPECT_EQ(copiedList->size(), 2);
    EXPECT_EQ(copiedList->at(0).getValue(), 42);
    EXPECT_EQ(copiedList->at(1).getValue(), 99);

    // Verify deep copy
    copiedList->at(0).setValue(100);
    EXPECT_EQ(list.at(0).getValue(), 42);
  }

  //=== Empty List Edge Cases ================================================

  TEST(FormLanguage_List, EmptyListIterators)
  {
    Rodin::FormLanguage::List<TestElement> list;

    EXPECT_EQ(list.begin(), list.end());
    EXPECT_EQ(list.cbegin(), list.cend());
  }

  TEST(FormLanguage_List, EmptyListCopy)
  {
    Rodin::FormLanguage::List<TestElement> list1;
    Rodin::FormLanguage::List<TestElement> list2(list1);

    EXPECT_TRUE(list2.empty());
    EXPECT_EQ(list2.size(), 0);
  }

  TEST(FormLanguage_List, ClearEmptyList)
  {
    Rodin::FormLanguage::List<TestElement> list;
    list.clear();

    EXPECT_TRUE(list.empty());
    EXPECT_EQ(list.size(), 0);
  }
}
