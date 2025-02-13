#include <gtest/gtest.h>

#include <thread>
#include <chrono>
#include <vector>

#include <Rodin/Threads/ThreadPool.h>

using Rodin::Threads::ThreadPool;

namespace Rodin::Tests::Unit
{
  //
  // Test basic construction and thread count retrieval
  //
  TEST(ThreadPoolTest, ConstructionAndGetThreadCount)
  {
    const size_t threadCount = 4;
    ThreadPool pool(threadCount);
    EXPECT_EQ(pool.getThreadCount(), threadCount);
  }

  //
  // Test resetting the thread pool with a new thread count
  //
  TEST(ThreadPoolTest, ResetThreadPool)
  {
    ThreadPool pool(2);
    EXPECT_EQ(pool.getThreadCount(), 2);
    pool.reset(6);
    EXPECT_EQ(pool.getThreadCount(), 6);
  }

  //
  // Test submitting a simple task using submit()
  //
  TEST(ThreadPoolTest, SubmitTask)
  {
    ThreadPool pool(4);
    // Submit a lambda that returns 42.
    auto future = pool.submit([]() -> int { return 42; });
    EXPECT_EQ(future.get(), 42);
  }

  //
  // Test that waitForTasks() correctly waits for submitted tasks to complete.
  // We submit a task that sleeps briefly before returning a value.
  //
  TEST(ThreadPoolTest, WaitForTasks)
  {
    ThreadPool pool(4);
    auto future = pool.submit([]() -> int {
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
      return 123;
    });
    // Wait for all tasks (should wait for the submitted task).
    pool.waitForTasks();
    EXPECT_EQ(future.get(), 123);
  }

  //
  // If the global thread pool is enabled (i.e. RODIN_MULTITHREADED is defined),
  // test that getGlobalThreadPool() returns a valid pool and that tasks can be submitted.
  //
#ifdef RODIN_MULTITHREADED
  TEST(ThreadPoolTest, GlobalThreadPool)
  {
    auto& globalPool = Rodin::Threads::getGlobalThreadPool();
    // Expect that the global pool has at least one thread.
    EXPECT_GT(globalPool.getThreadCount(), 0u);

    auto future = globalPool.submit([]() -> int { return 999; });
    EXPECT_EQ(future.get(), 999);
  }
#endif
}
