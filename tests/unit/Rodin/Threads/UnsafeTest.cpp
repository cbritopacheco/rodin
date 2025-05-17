#include <gtest/gtest.h>

#include <string>
#include <thread>
#include <atomic>
#include <chrono>

#include <Rodin/Threads/Unsafe.h>

using Rodin::Threads::Unsafe;

namespace Rodin::Tests::Unit
{

#if defined(RODIN_THREAD_SAFE) && !defined(NDEBUG)
  // In thread-safe debug builds the exception is thrown.
  using UnsafeAccessExceptionInt =
    Rodin::Threads::UnsafeAccessException<Unsafe<int>, const char*>;
#endif

  // -----------------------------------------------------------------------------
  // Basic Functionality Tests
  // -----------------------------------------------------------------------------
  TEST(UnsafeTest, DefaultConstructorAndRead)
  {
    // Create a default-constructed Unsafe<int> and then assign a value.
    Unsafe<int> u;
    u = 10;
    EXPECT_EQ(u.read(), 10);
  }

  TEST(UnsafeTest, CopyConstructorAndCopyAssignment)
  {
    Unsafe<int> u1(100);
    Unsafe<int> u2(u1); // copy constructor
    EXPECT_EQ(u2.read(), 100);

    u2 = 200; // copy assignment operator
    EXPECT_EQ(u2.read(), 200);
  }

  TEST(UnsafeTest, MoveConstructorAndMoveAssignment)
  {
    // Test with a nontrivial resource (std::string).
    Unsafe<std::string> u1(std::string("Hello"));
    Unsafe<std::string> u2(std::move(u1)); // move constructor
    EXPECT_EQ(u2.read(), "Hello");

    u2 = std::string("World"); // move assignment operator
    EXPECT_EQ(u2.read(), "World");
  }

  TEST(UnsafeTest, WriteModification)
  {
    // Test modifying the encapsulated resource using the write() member.
    Unsafe<int> u(0);
    u.write([](int& value) {
      value = 12345;
    });
    EXPECT_EQ(u.read(), 12345);
  }

  TEST(UnsafeTest, AssignmentOperators)
  {
    // Test assignment from a resource and from another Unsafe.
    Unsafe<int> u;
    u = 55;
    EXPECT_EQ(u.read(), 55);
    u = 66;
    EXPECT_EQ(u.read(), 66);

    Unsafe<int> u2(77);
    u = u2;
    EXPECT_EQ(u.read(), 77);

    Unsafe<int> u3(88);
    u = std::move(u3);
    EXPECT_EQ(u.read(), 88);
  }

  // -----------------------------------------------------------------------------
  // Thread Safety / Unsafe Access Exception Tests
  // -----------------------------------------------------------------------------
#if defined(RODIN_THREAD_SAFE) && !defined(NDEBUG)
  // In thread-safe debug builds, reentrant or concurrent writes are not allowed.
  TEST(UnsafeTest, ReentrantWriteThrows)
  {
    Unsafe<int> u(0);
    // Attempt to call write() from within an active write() call.
    EXPECT_THROW({
      u.write([&u](int& value) {
        // This nested call should detect that a write is already in progress.
        u.write([](int& innerValue) { innerValue = 999; });
      });
    }, Alert::Exception);
  }

  TEST(UnsafeTest, ConcurrentWriteThrows)
  {
    // Use two threads to simulate concurrent write attempts.
    Unsafe<int> u(0);
    std::atomic<bool> exceptionCaught{false};

    std::thread t1([&u]() {
      // Long-running write operation.
      u.write([](int& value) {
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        value = 42;
      });
    });

    // Ensure t1 starts and holds the write lock.
    std::this_thread::sleep_for(std::chrono::milliseconds(50));

    std::thread t2([&u, &exceptionCaught]() {
      try {
        u.write([](int& value) { value = 100; });
      }
      catch(const Alert::Exception& /*e*/) {
        exceptionCaught = true;
      }
    });

    t1.join();
    t2.join();
    EXPECT_TRUE(exceptionCaught);
  }
#endif

#if !defined(RODIN_THREAD_SAFE) || defined(NDEBUG)
  // In non-thread-safe mode (or in release builds), no exception is thrown.
  TEST(UnsafeTest, ReentrantWriteNoThrowWhenNoThreadSafety)
  {
    Unsafe<int> u(0);
    EXPECT_NO_THROW({
      u.write([&u](int& value) {
        u.write([](int& innerValue) { innerValue = 777; });
      });
    });
    EXPECT_EQ(u.read(), 777);
  }
#endif

}
