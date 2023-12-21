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
      template <class T>
      using MultiFuture = BS::multi_future<T>;

      ThreadPool(size_t numThreads)
        : m_pool(numThreads)
      {}

      template <class ... Args>
      auto submit(Args... args)
      {
        return m_pool.submit(std::forward<Args>(args)...);
      }

      template <class ... Args>
      auto pushLoop(Args... args)
      {
        return m_pool.push_loop(std::forward<Args>(args)...);
      }

      void waitForTasks()
      {
        m_pool.wait_for_tasks();
      }

      inline
      auto getThreadCount() const
      {
        return m_pool.get_thread_count();
      }

    private:
      BS::thread_pool m_pool;
  };

#ifdef RODIN_MULTITHREADED
  namespace Internal
  {
    static ThreadPool globalThreadPool(RODIN_THREADPOOL_GLOBALTHREADPOOL_CONCURRENCY);
  }

  inline
  static ThreadPool& getGlobalThreadPool()
  {
    return Internal::globalThreadPool;
  }
#endif
}

#endif

