#include <gtest/gtest.h>

#include <thread>
#include <chrono>
#include <vector>

#include <Rodin/Threads/Shared.h>

using Rodin::Threads::Shared;

namespace Rodin::Tests::Unit
{
  // -----------------------------------------------------------------------------
  // Construction Tests
  // -----------------------------------------------------------------------------
  TEST(SharedTest, DefaultConstructor)
  {
    // Using std::string since default construction of int is indeterminate.
    Shared<std::string> s;
    // A default-constructed std::string is empty.
    EXPECT_EQ(s.read(), "");
  }

  TEST(SharedTest, ConstructFromConstResource)
  {
    Shared<std::string> s("Hello");
    EXPECT_EQ(s.read(), "Hello");
  }

  TEST(SharedTest, ConstructFromRvalueResource)
  {
    Shared<std::string> s(std::string("World"));
    EXPECT_EQ(s.read(), "World");
  }

  TEST(SharedTest, CopyConstruction)
  {
    Shared<int> s1(42);
    Shared<int> s2(s1);
    EXPECT_EQ(s2.read(), 42);
  }

  TEST(SharedTest, MoveConstruction)
  {
    Shared<int> s1(77);
    Shared<int> s2(std::move(s1));
    // s2 should now hold the value.
    EXPECT_EQ(s2.read(), 77);
  }

  // -----------------------------------------------------------------------------
  // Assignment Tests
  // -----------------------------------------------------------------------------
  TEST(SharedTest, CopyAssignmentOperator)
  {
    Shared<std::string> s1("Alpha");
    Shared<std::string> s2("Beta");
    s2 = s1;
    EXPECT_EQ(s2.read(), "Alpha");
  }

  TEST(SharedTest, MoveAssignmentOperator)
  {
    Shared<std::string> s1("Gamma");
    Shared<std::string> s2("Delta");
    s2 = std::move(s1);
    EXPECT_EQ(s2.read(), "Gamma");
  }

  // -----------------------------------------------------------------------------
  // Direct Accessor Tests
  // -----------------------------------------------------------------------------
  TEST(SharedTest, DirectWriteAccessor)
  {
    Shared<int> s(10);
    // Modify the underlying resource directly.
    s.write() = 25;
    EXPECT_EQ(s.read(), 25);
  }

  TEST(SharedTest, ReadAccessorReturnsConstReference)
  {
    Shared<int> s(33);
    const int& value = s.read();
    EXPECT_EQ(value, 33);
  }

  // -----------------------------------------------------------------------------
  // Lambda Accessor Tests
  // -----------------------------------------------------------------------------
  TEST(SharedTest, ReadWithLambda)
  {
    Shared<int> s(100);
    int captured = 0;
    s.read([&captured](const int& value) {
      captured = value;
    });
    EXPECT_EQ(captured, 100);
  }

  TEST(SharedTest, WriteWithLambda)
  {
    Shared<int> s(50);
    s.write([](int& value) {
      value += 20;
    });
    EXPECT_EQ(s.read(), 70);
  }

  // -----------------------------------------------------------------------------
  // Thread Safety Tests (only compiled when RODIN_THREAD_SAFE is defined)
  // -----------------------------------------------------------------------------
#ifdef RODIN_THREAD_SAFE
  TEST(SharedTest, ConcurrentRead)
  {
    // Create a Shared<int> initialized to 42.
    Shared<int> s(42);
    const int numThreads = 10;
    std::vector<std::thread> readers;
    std::atomic<int> sum{0};

    // Launch multiple threads that read the value concurrently.
    for (int i = 0; i < numThreads; i++)
    {
      readers.emplace_back([&s, &sum]() {
        s.read([&sum](const int& value) {
          sum += value;
        });
      });
    }
    for (auto& t : readers)
    {
      t.join();
    }
    EXPECT_EQ(sum.load(), 42 * numThreads);
  }

  TEST(SharedTest, ThreadSafeConcurrentWrite)
  {
    // Test that exclusive write access serializes modifications.
    Shared<int> s(0);
    const int numThreads = 10;
    const int incrementsPerThread = 100;
    std::vector<std::thread> writers;

    for (int i = 0; i < numThreads; i++)
    {
      writers.emplace_back([&s]() {
        for (int j = 0; j < incrementsPerThread; j++)
        {
          s.write([](int& value) {
            value++;  // increment safely under exclusive lock
          });
        }
      });
    }
    for (auto& t : writers)
    {
      t.join();
    }
    EXPECT_EQ(s.read(), numThreads * incrementsPerThread);
  }
#endif
}

