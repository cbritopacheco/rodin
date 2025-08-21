#include <gtest/gtest.h>

#include <Rodin/Moveable.h>

using namespace Rodin;

namespace Rodin::Tests::Unit
{
  // Test implementation of Moveable for testing purposes
  class TestMoveable : public Moveable
  {
    public:
      TestMoveable(int value) : m_value(value), m_moved(false) {}

      Moveable* move() noexcept override
      {
        m_moved = true;
        return this;
      }

      int getValue() const { return m_value; }
      bool wasMoved() const { return m_moved; }

    private:
      int m_value;
      bool m_moved;
  };

  // Test implementation with resource management
  class TestMoveableWithResource : public Moveable
  {
    public:
      TestMoveableWithResource(std::unique_ptr<int> resource) 
        : m_resource(std::move(resource)), m_moved(false) {}

      Moveable* move() noexcept override
      {
        m_moved = true;
        return this;
      }

      int* getResource() const { return m_resource.get(); }
      bool wasMoved() const { return m_moved; }
      bool hasResource() const { return m_resource != nullptr; }

    private:
      std::unique_ptr<int> m_resource;
      bool m_moved;
  };

  //=== Basic Functionality Tests ============================================

  TEST(Rodin_Moveable, BasicMove)
  {
    TestMoveable original(42);
    EXPECT_FALSE(original.wasMoved());

    Moveable* moved = original.move();

    EXPECT_EQ(moved, &original);  // Should return pointer to self
    EXPECT_TRUE(original.wasMoved());
    EXPECT_EQ(original.getValue(), 42);  // Value should remain
  }

  TEST(Rodin_Moveable, PolymorphicMove)
  {
    std::unique_ptr<Moveable> moveable = std::make_unique<TestMoveable>(99);
    TestMoveable* typedOriginal = dynamic_cast<TestMoveable*>(moveable.get());

    ASSERT_NE(typedOriginal, nullptr);
    EXPECT_FALSE(typedOriginal->wasMoved());

    Moveable* moved = moveable->move();

    EXPECT_EQ(moved, moveable.get());
    EXPECT_TRUE(typedOriginal->wasMoved());
    EXPECT_EQ(typedOriginal->getValue(), 99);
  }

  TEST(Rodin_Moveable, MoveWithResource)
  {
    auto resource = std::make_unique<int>(123);
    int* resourcePtr = resource.get();

    TestMoveableWithResource original(std::move(resource));
    EXPECT_TRUE(original.hasResource());
    EXPECT_EQ(*original.getResource(), 123);
    EXPECT_EQ(original.getResource(), resourcePtr);
    EXPECT_FALSE(original.wasMoved());

    Moveable* moved = original.move();

    EXPECT_EQ(moved, &original);
    EXPECT_TRUE(original.wasMoved());
    EXPECT_TRUE(original.hasResource());  // Resource should still be there
    EXPECT_EQ(*original.getResource(), 123);
  }

  TEST(Rodin_Moveable, MultipleMoves)
  {
    TestMoveable original(77);

    // First move
    Moveable* moved1 = original.move();
    EXPECT_EQ(moved1, &original);
    EXPECT_TRUE(original.wasMoved());

    // Second move should still work
    Moveable* moved2 = original.move();
    EXPECT_EQ(moved2, &original);
    EXPECT_TRUE(original.wasMoved());
  }

  TEST(Rodin_Moveable, NoExceptSpecification)
  {
    // Verify that move() is declared as noexcept
    TestMoveable original(1);

    // This should compile without issues since move() is noexcept
    EXPECT_NO_THROW({
      Moveable* moved = original.move();
      (void)moved;  // Suppress unused variable warning
    });
  }

  TEST(Rodin_Moveable, OwnershipTransfer)
  {
    // Test that the caller assumes ownership responsibility
    TestMoveable original(555);
    Moveable* moved = original.move();

    // The returned pointer should be the same object
    EXPECT_EQ(moved, &original);

    // The object should be marked as moved
    EXPECT_TRUE(original.wasMoved());

    // The caller now has responsibility for the object
    // (This is conceptual - in this test we don't actually transfer ownership
    // since we're using stack allocation, but the interface contract is tested)
  }
}
