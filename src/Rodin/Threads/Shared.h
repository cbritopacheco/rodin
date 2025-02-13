#ifndef RODIN_SHARED_H
#define RODIN_SHARED_H

#include <shared_mutex>

#include "Rodin/Configure.h"

namespace Rodin::Threads
{
  /**
   * @brief A thread-safe wrapper for shared resources.
   *
   * The Shared class encapsulates a resource and provides thread-safe read and
   * write access. When RODIN_THREAD_SAFE is defined, access to the resource is
   * guarded by a std::shared_mutex. In non-thread-safe builds, the class
   * provides direct access to the underlying resource.
   *
   * @tparam Resource The type of the resource to be encapsulated.
   */
  template <class Resource>
  class Shared
  {
    public:
      /**
       * @brief Default constructor.
       *
       * Constructs an empty Shared object.
       */
      constexpr
      Shared()
      {}

      /**
       * @brief Constructs a Shared object from a const reference.
       *
       * @param resource The resource to encapsulate.
       */
      constexpr
      Shared(const Resource& resource)
        : m_resource(resource)
      {}

      /**
       * @brief Constructs a Shared object by moving a resource.
       *
       * @param resource The resource to encapsulate.
       */
      constexpr
      Shared(Resource&& resource)
        : m_resource(std::move(resource))
      {}

      /**
       * @brief Move constructor.
       *
       * Transfers the resource from another Shared object.
       *
       * @param other The Shared object to move from.
       */
      constexpr
      Shared(Shared&& other)
        : m_resource(std::move(other.m_resource))
      {}

      /**
       * @brief Copy constructor.
       *
       * Creates a Shared object by copying the resource from another.
       *
       * @param other The Shared object to copy from.
       */
      constexpr
      Shared(const Shared& other)
        : m_resource(other.m_resource)
      {}

      /**
       * @brief Copy assignment operator from a Resource.
       *
       * This operation is deleted to prevent direct assignment from a resource.
       */
      constexpr
      Shared& operator=(const Resource&) = delete;

      /**
       * @brief Move assignment operator from a Resource.
       *
       * This operation is deleted to prevent direct move-assignment from a resource.
       */
      constexpr
      Shared& operator=(Resource&&) = delete;

      /**
       * @brief Copy assignment operator.
       *
       * Assigns the resource from another Shared object.
       * In thread-safe builds, this operation is protected by locking.
       *
       * @param other The Shared object to copy from.
       * @return A reference to this Shared object.
       */
      constexpr
      Shared& operator=(const Shared& other)
      {
        if (this != &other)
        {
#ifdef RODIN_THREAD_SAFE
          lock();
          m_resource = other.m_resource;
          unlock();
#else
          m_resource = other.m_resource;
#endif
        }
        return *this;
      }

      /**
       * @brief Move assignment operator.
       *
       * Transfers the resource from another Shared object.
       * In thread-safe builds, this operation is protected by locking.
       *
       * @param other The Shared object to move from.
       * @return A reference to this Shared object.
       */
      constexpr
      Shared& operator=(Shared&& other)
      {
#ifdef RODIN_THREAD_SAFE
        lock();
        other.lock();
        m_resource = std::move(other.m_resource);
        other.unlock();
        unlock();
#else
        m_resource = std::move(other.m_resource);
#endif
        return *this;
      }

      /**
       * @brief Provides non-const (write) access to the encapsulated resource.
       *
       * @return A reference to the encapsulated resource.
       */
      constexpr
      Resource& write()
      {
        return m_resource;
      }

      /**
       * @brief Provides read-only access to the encapsulated resource.
       *
       * @return A const reference to the encapsulated resource.
       */
      constexpr
      const Resource& read() const
      {
        return m_resource;
      }

      /**
       * @brief Executes a callable with read-only access to the resource.
       *
       * The provided callable is invoked with a const reference to the
       * resource. In thread-safe builds, the shared lock is acquired during
       * the execution of the callable.
       *
       * @tparam F The type of the callable.
       * @param f The callable to execute.
       * @return A const reference to this Shared object.
       */
      template <class F>
      constexpr
      const Shared& read(F&& f) const
      {
        static_assert(std::is_invocable_v<F, const Resource&>);
#ifdef RODIN_THREAD_SAFE
        m_lock.lock_shared();
        f(static_cast<const Resource&>(m_resource));
        m_lock.unlock_shared();
#else
        f(static_cast<const Resource&>(m_resource));
#endif
        return *this;
      }

      /**
       * @brief Executes a callable with write access to the resource.
       *
       * The provided callable is invoked with a non-const reference to the
       * resource. In thread-safe builds, the exclusive lock is acquired during
       * the execution of the callable.
       *
       * @tparam F The type of the callable.
       * @param f The callable to execute.
       * @return A reference to this Shared object.
       */
      template <class F>
      constexpr
      Shared& write(F&& f)
      {
        static_assert(std::is_invocable_v<F, Resource&>);
#ifdef RODIN_THREAD_SAFE
        lock();
        f(m_resource);
        unlock();
#else
        f(m_resource);
#endif
        return *this;
      }

#ifdef RODIN_THREAD_SAFE
    protected:
      /**
       * @brief Acquires an exclusive lock on the resource.
       *
       * This function is used internally to guard write operations.
       */
      constexpr
      void lock()
      {
        m_lock.lock();
      }

      /**
       * @brief Releases the exclusive lock on the resource.
       *
       * This function is used internally after write operations are completed.
       */
      constexpr
      void unlock()
      {
        m_lock.unlock();
      }
#endif

    private:
#ifdef RODIN_THREAD_SAFE
      /// Shared mutex used for thread-safe access.
      mutable std::shared_mutex m_lock;
#endif
      /// The encapsulated resource.
      Resource m_resource;
  };
}

#endif


