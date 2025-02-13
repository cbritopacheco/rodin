#ifndef RODIN_MUTABLE_H
#define RODIN_MUTABLE_H

#include "Rodin/Configure.h"

#include "Mutex.h"

namespace Rodin::Threads
{
  /**
   * @brief A thread-safe mutable wrapper for a resource.
   *
   * The Mutable class encapsulates a resource and provides both read and write access
   * with optional thread-safety. When RODIN_THREAD_SAFE is defined, all modifying and
   * accessing operations are protected by a lock of type Lock (defaulting to Mutex).
   *
   * @tparam Resource The type of the resource to be encapsulated.
   * @tparam Lock The type of the lock to use for thread-safety (default is Mutex).
   */
  template <class Resource, class Lock = Mutex>
  class Mutable
  {
    public:
      /**
       * @brief Default constructor.
       *
       * Constructs an empty Mutable object.
       */
      constexpr
      Mutable()
      {}

      /**
       * @brief Constructs a Mutable object from a const reference.
       *
       * @param resource The resource to be encapsulated.
       */
      constexpr
      Mutable(const Resource& resource)
        : m_resource(resource)
      {}

      /**
       * @brief Constructs a Mutable object by moving a resource.
       *
       * @param resource The resource to be encapsulated.
       */
      constexpr
      Mutable(Resource&& resource)
        : m_resource(std::move(resource))
      {}

      /**
       * @brief Move constructor.
       *
       * Transfers the encapsulated resource from another Mutable object.
       *
       * @param other The Mutable object to move from.
       */
      constexpr
      Mutable(Mutable&& other)
        : m_resource(std::move(other.m_resource))
      {}

      /**
       * @brief Copy constructor.
       *
       * Copies the encapsulated resource from another Mutable object.
       *
       * @param other The Mutable object to copy from.
       */
      constexpr
      Mutable(const Mutable& other)
        : m_resource(other.m_resource)
      {}

      /**
       * @brief Assignment operator from a Resource.
       *
       * Assigns a new value to the encapsulated resource.
       * In thread-safe builds, the assignment is performed under lock protection.
       *
       * @param resource The new value to assign to the resource.
       * @return A reference to this Mutable object.
       */
      constexpr
      Mutable& operator=(const Resource& resource)
      {
#ifdef RODIN_THREAD_SAFE
        lock();
        m_resource = resource;
        unlock();
#else
        m_resource = resource;
#endif
        return *this;
      }

      /**
       * @brief Move assignment operator from a Resource.
       *
       * Moves a new value into the encapsulated resource.
       * In thread-safe builds, the assignment is performed under lock protection.
       *
       * @param resource The resource to move into this object.
       * @return A reference to this Mutable object.
       */
      constexpr
      Mutable& operator=(Resource&& resource)
      {
#ifdef RODIN_THREAD_SAFE
        lock();
        m_resource = std::move(resource);
        unlock();
#else
        m_resource = std::move(resource);
#endif
        return *this;
      }

      /**
       * @brief Copy assignment operator.
       *
       * Assigns the resource from another Mutable object.
       * In thread-safe builds, the assignment is performed under lock protection.
       *
       * @param other The Mutable object to copy from.
       * @return A reference to this Mutable object.
       */
      constexpr
      Mutable& operator=(const Mutable& other)
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
       * Transfers the resource from another Mutable object into this one.
       * In thread-safe builds, the assignment is performed under lock protection.
       *
       * @param other The Mutable object to move from.
       * @return A reference to this Mutable object.
       */
      constexpr
      Mutable& operator=(Mutable&& other)
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
       * @return A reference to the resource.
       */
      constexpr
      Resource& write()
      {
        return m_resource;
      }

      /**
       * @brief Provides read-only access to the encapsulated resource.
       *
       * @return A const reference to the resource.
       */
      constexpr
      const Resource& read() const
      {
        return m_resource;
      }

      /**
       * @brief Executes a callable with read-only access to the resource.
       *
       * The callable is invoked with a const reference to the resource.
       * Note: The parameter name for the callable is intentionally omitted.
       *
       * @tparam F The type of the callable.
       * @param f The callable to execute.
       * @return A reference to this Mutable object.
       */
      template <class F>
      constexpr
      Mutable& read(F&& f) const
      {
        static_assert(std::is_invocable_v<F, const Resource&>);
        f(m_resource);
        return const_cast<Mutable&>(*this);
      }

      /**
       * @brief Executes a callable with write access to the resource.
       *
       * The callable is invoked with a non-const reference to the resource.
       * In thread-safe builds, the operation is performed under lock protection.
       *
       * @tparam F The type of the callable.
       * @param f The callable to execute.
       * @return A reference to this Mutable object.
       */
      template <class F>
      constexpr
      Mutable& write(F&& f)
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
       * @brief Acquires the lock for exclusive access.
       *
       * This function is used internally to protect write operations.
       */
      constexpr
      void lock()
      {
        m_lock.lock();
      }

      /**
       * @brief Releases the exclusive lock.
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
      /// The lock used to ensure thread-safety.
      Lock m_lock;
#endif
      /// The encapsulated resource.
      Resource m_resource;
  };
}

#endif

