#ifndef RODIN_THREADPOOL_H
#define RODIN_THREADPOOL_H

#include <thread>
#include <BS_thread_pool.hpp>

#define RODIN_THREADPOOL_GLOBALTHREADPOOL_CONCURRENCY std::thread::hardware_concurrency()

namespace Rodin::Threads
{
  class ThreadPool
  {
    public:
      ThreadPool(size_t numThreads)
        : m_pool(numThreads)
      {}

      template <class ... Args>
      void pushLoop(Args... args)
      {
        m_pool.push_loop(std::forward<Args>(args)...);
      }

      void waitForTasks()
      {
        m_pool.wait_for_tasks();
      }

    private:
      BS::thread_pool m_pool;
  };

#ifdef RODIN_MULTITHREADED
  static ThreadPool globalThreadPool(RODIN_THREADPOOL_GLOBALTHREADPOOL_CONCURRENCY);
#endif
}

#endif

