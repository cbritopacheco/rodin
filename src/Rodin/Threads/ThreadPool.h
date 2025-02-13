#ifndef RODIN_THREADPOOL_H
#define RODIN_THREADPOOL_H

#include <thread>
#include <BS_thread_pool.hpp>

/**
 * @def RODIN_THREADPOOL_GLOBALTHREADPOOL_CONCURRENCY
 * @brief Defines the default number of threads for the global thread pool.
 *
 * This macro uses std::thread::hardware_concurrency() to set the default
 * thread count for the global thread pool, which represents the number of
 * concurrent threads supported by the hardware.
 */
#define RODIN_THREADPOOL_GLOBALTHREADPOOL_CONCURRENCY std::thread::hardware_concurrency()

namespace Rodin::Threads
{
  /**
   * @brief A thread pool wrapper for managing concurrent tasks.
   *
   * The ThreadPool class wraps around BS::thread_pool and provides a simplified
   * interface for submitting tasks, pushing loops of tasks, and waiting for all
   * tasks to complete.
   */
  class ThreadPool
  {
    public:
      /**
       * @brief Alias for a multi_future type returned by task submissions.
       *
       * @tparam T The type of the value produced by the tasks.
       */
      template <class T>
      using MultiFuture = BS::multi_future<T>;

      /**
       * @brief Constructs a ThreadPool with a specified number of threads.
       *
       * @param threadCount The number of threads to create in the pool.
       */
      ThreadPool(size_t threadCount)
        : m_pool(threadCount)
      {}

      /**
       * @brief Resets the thread pool with a new thread count.
       *
       * This function reinitializes the thread pool to use the specified number
       * of threads.
       *
       * @param threadCount The new number of threads for the pool.
       * @return A reference to this ThreadPool instance.
       */
      ThreadPool& reset(size_t threadCount)
      {
        m_pool.reset(threadCount);
        return *this;
      }

      /**
       * @brief Submits a task to the thread pool.
       *
       * Forwards the provided arguments to BS::thread_pool::submit.
       *
       * @tparam Args The types of the arguments to be forwarded.
       * @param args The arguments for the task to be submitted.
       * @return A future or multi_future representing the submitted task.
       */
      template <class ... Args>
      auto submit(Args... args)
      {
        return m_pool.submit(std::forward<Args>(args)...);
      }

      /**
       * @brief Pushes a loop of tasks to the thread pool.
       *
       * Forwards the provided arguments to BS::thread_pool::push_loop.
       *
       * @tparam Args The types of the arguments to be forwarded.
       * @param args The arguments for the loop of tasks to be pushed.
       * @return A future or multi_future representing the pushed tasks.
       */
      template <class ... Args>
      void pushLoop(Args... args)
      {
        return m_pool.push_loop(std::forward<Args>(args)...);
      }

      /**
       * @brief Waits for all tasks in the thread pool to complete.
       */
      void waitForTasks()
      {
        m_pool.wait_for_tasks();
      }

      /**
       * @brief Retrieves the current number of threads in the pool.
       *
       * @return The number of threads currently used by the thread pool.
       */
      auto getThreadCount() const
      {
        return m_pool.get_thread_count();
      }

    private:
      /// The underlying BS::thread_pool instance.
      BS::thread_pool m_pool;
  };

#ifdef RODIN_MULTITHREADED
  namespace Internal
  {
    /**
     * @brief Global thread pool instance.
     *
     * This static instance of ThreadPool is initialized with the number of threads
     * equal to the hardware concurrency reported by std::thread.
     */
    static ThreadPool globalThreadPool(RODIN_THREADPOOL_GLOBALTHREADPOOL_CONCURRENCY);
  }

  /**
   * @brief Retrieves the global thread pool.
   *
   * This function returns a reference to a global ThreadPool instance which can be used
   * throughout the application for managing concurrent tasks.
   *
   * @return A reference to the global ThreadPool.
   */
  inline
  static ThreadPool& getGlobalThreadPool()
  {
    return Internal::globalThreadPool;
  }
#endif
}

#endif

