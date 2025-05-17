#include <gtest/gtest.h>

#include <thread>
#include <chrono>
#include <vector>

#include <Rodin/Threads/Mutable.h>

using Rodin::Threads::Mutable;

namespace Rodin::Tests::Unit
{
  // -----------------------------------------------------------------------------
  // Construction Tests
  // -----------------------------------------------------------------------------
  TEST(MutableTest, DefaultConstructor)
  {
    // For types like std::string, default construction yields an empty string.
    Mutable<std::string> m;
    EXPECT_EQ(m.read(), "");
  }

  TEST(MutableTest, ConstructFromConstResource)
  {
    Mutable<std::string> m("Hello");
    EXPECT_EQ(m.read(), "Hello");
  }

  TEST(MutableTest, ConstructFromRvalueResource)
  {
    Mutable<std::string> m(std::string("World"));
    EXPECT_EQ(m.read(), "World");
  }

  TEST(MutableTest, CopyConstructor)
  {
    Mutable<int> m1(42);
    Mutable<int> m2(m1);
    EXPECT_EQ(m2.read(), 42);
  }

  TEST(MutableTest, MoveConstructor)
  {
    Mutable<std::string> m1("Mutable");
    Mutable<std::string> m2(std::move(m1));
    EXPECT_EQ(m2.read(), "Mutable");
  }

  // -----------------------------------------------------------------------------
  // Assignment Tests
  // -----------------------------------------------------------------------------
  TEST(MutableTest, AssignmentFromResourceCopy)
  {
    Mutable<std::string> m;
    m = std::string("Test");
    EXPECT_EQ(m.read(), "Test");
  }

  TEST(MutableTest, AssignmentFromResourceMove)
  {
    Mutable<std::string> m;
    m = std::string("MoveTest");
    EXPECT_EQ(m.read(), "MoveTest");
  }

  TEST(MutableTest, CopyAssignmentOperator)
  {
    Mutable<int> m1(55);
    Mutable<int> m2;
    m2 = m1;
    EXPECT_EQ(m2.read(), 55);
  }

  TEST(MutableTest, MoveAssignmentOperator)
  {
    Mutable<std::string> m1("Data");
    Mutable<std::string> m2;
    m2 = std::move(m1);
    EXPECT_EQ(m2.read(), "Data");
  }

  // -----------------------------------------------------------------------------
  // Direct Accessor Tests
  // -----------------------------------------------------------------------------
  TEST(MutableTest, DirectAccessors)
  {
    Mutable<int> m(10);
    // Use non-const write() to modify the resource.
    m.write() = 20;
    EXPECT_EQ(m.read(), 20);
  }

  // -----------------------------------------------------------------------------
  // Lambda Accessor Tests
  // -----------------------------------------------------------------------------
  TEST(MutableTest, ReadWithLambda)
  {
    Mutable<int> m(100);
    int captured = 0;
    m.read([&captured](const int& value) {
      captured = value;
    });
    EXPECT_EQ(captured, 100);
  }

  TEST(MutableTest, WriteWithLambda)
  {
    Mutable<int> m(50);
    m.write([](int& value) {
      value += 25;
    });
    EXPECT_EQ(m.read(), 75);
  }

  // -----------------------------------------------------------------------------
  // Thread Safety Tests (only compiled if RODIN_THREAD_SAFE is defined)
  // -----------------------------------------------------------------------------
#ifdef RODIN_THREAD_SAFE
  TEST(MutableTest, ThreadSafeConcurrentWrite)
  {
    // This test verifies that concurrent modifications are serialized.
    Mutable<int> m(0);
    const int numThreads = 10;
    const int incrementsPerThread = 1000;
    std::vector<std::thread> threads;

    for (int i = 0; i < numThreads; i++)
    {
      threads.emplace_back([&m]() {
        for (int j = 0; j < incrementsPerThread; j++)
        {
          m.write([](int& value) { value++; });
        }
      });
    }

    for (auto& t : threads)
    {
      t.join();
    }

    EXPECT_EQ(m.read(), numThreads * incrementsPerThread);
  }

#endif // RODIN_THREAD_SAFE
}
