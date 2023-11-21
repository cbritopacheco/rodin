#ifndef RODIN_THREADPOOL_H
#define RODIN_THREADPOOL_H

#include "BS_thread_pool.hpp"

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
}

#endif

