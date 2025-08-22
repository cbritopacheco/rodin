#include <gtest/gtest.h>

#include <Rodin/Copyable.h>

using namespace Rodin;

namespace Rodin::Tests::Unit
{
  // Test implementation of Copyable for testing purposes
  class TestCopyable : public Copyable
  {
    public:
      TestCopyable(int value) : m_value(value) {}

      TestCopyable* copy() const noexcept override
      {
        return new TestCopyable(m_value);
      }

      int getValue() const { return m_value; }

    private:
      int m_value;
  };

  // Test implementation that throws during copy (should not happen with noexcept)
  class TestCopyableWithValue : public Copyable
  {
    public:
      TestCopyableWithValue(const std::string& value) : m_value(value) {}

      TestCopyableWithValue* copy() const noexcept override
      {
        return new TestCopyableWithValue(m_value);
      }

      const std::string& getValue() const { return m_value; }

    private:
      std::string m_value;
  };

  //=== Basic Functionality Tests ============================================

  TEST(Rodin_Copyable, BasicCopy)
  {
    TestCopyable original(42);
    std::unique_ptr<Copyable> copied(original.copy());

    ASSERT_NE(copied.get(), nullptr);

    // Downcast to access test-specific methods
    TestCopyable* typedCopy = dynamic_cast<TestCopyable*>(copied.get());
    ASSERT_NE(typedCopy, nullptr);
    EXPECT_EQ(typedCopy->getValue(), 42);
  }

  TEST(Rodin_Copyable, CopyIndependence)
  {
    TestCopyable original(100);
    std::unique_ptr<TestCopyable> copied(original.copy());

    ASSERT_NE(copied.get(), nullptr);
    EXPECT_EQ(copied->getValue(), 100);

    // Verify they are different objects
    EXPECT_NE(&original, copied.get());
  }

  TEST(Rodin_Copyable, PolymorphicCopy)
  {
    std::unique_ptr<Copyable> original = std::make_unique<TestCopyable>(77);
    std::unique_ptr<Copyable> copied(original->copy());

    ASSERT_NE(copied.get(), nullptr);

    // Verify polymorphic behavior
    TestCopyable* originalTyped = dynamic_cast<TestCopyable*>(original.get());
    TestCopyable* copiedTyped = dynamic_cast<TestCopyable*>(copied.get());

    ASSERT_NE(originalTyped, nullptr);
    ASSERT_NE(copiedTyped, nullptr);
    EXPECT_EQ(originalTyped->getValue(), copiedTyped->getValue());
  }

  TEST(Rodin_Copyable, StringValueCopy)
  {
    TestCopyableWithValue original("test_string");
    std::unique_ptr<TestCopyableWithValue> copied(original.copy());

    ASSERT_NE(copied.get(), nullptr);
    EXPECT_EQ(copied->getValue(), "test_string");
    EXPECT_NE(&original, copied.get());
  }

  TEST(Rodin_Copyable, EmptyStringCopy)
  {
    TestCopyableWithValue original("");
    std::unique_ptr<TestCopyableWithValue> copied(original.copy());

    ASSERT_NE(copied.get(), nullptr);
    EXPECT_EQ(copied->getValue(), "");
  }

  TEST(Rodin_Copyable, NoExceptSpecification)
  {
    // Verify that copy() is declared as noexcept
    TestCopyable original(1);

    // This should compile without issues since copy() is noexcept
    EXPECT_NO_THROW({
      std::unique_ptr<TestCopyable> copied(original.copy());
    });
  }
}
