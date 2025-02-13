#ifndef RODIN_UNSAFE_H
#define RODIN_UNSAFE_H

#include "Rodin/Configure.h"
#include "Rodin/Alert/MemberFunctionException.h"

namespace Rodin::Threads
{
  /**
   * @brief Exception thrown when unsafe concurrent access is detected.
   *
   * This exception is thrown by the Unsafe class when multiple threads attempt
   * to write to the encapsulated resource concurrently (in debug mode with
   * thread safety enabled).
   *
   * @tparam T The type of the object (resource) that triggered the exception.
   * @tparam FuncName The type used for representing the name of the member function.
   */
  template <class T, class FuncName>
  class UnsafeAccessException : public Alert::MemberFunctionException<T, FuncName>
  {
    public:
      /// Alias for the parent exception class.
      using Parent = Alert::MemberFunctionException<T, FuncName>;

      /**
       * @brief Constructs an UnsafeAccessException.
       *
       * Initializes the exception with the given object and function name.
       *
       * @param cls The object (resource) being accessed unsafely.
       * @param funcName The name of the member function where the unsafe
       * access occurred.
       */
      UnsafeAccessException(const T& cls, const FuncName& funcName)
        : Parent(cls, funcName)
      {
        *this << Alert::Text("Rodin::Threads::UnsafeAccessException").setBold().setUnderline() << ": ";
      }

      virtual void raise() const override
      {
        Parent::raise();
      }
  };

  /**
   * @brief A wrapper that manages unsafe access to a resource.
   *
   * The Unsafe class encapsulates a resource and provides controlled access to
   * it. In builds where RODIN_THREAD_SAFE is defined and in debug mode, write
   * operations are guarded to prevent concurrent modifications which could
   * lead to data races. In other configurations, it forwards operations
   * directly to the underlying resource.
   *
   * @tparam Resource The type of the resource to be encapsulated.
   */
  template <class Resource>
  class Unsafe
  {
    public:
#if defined(RODIN_THREAD_SAFE) && !defined(NDEBUG)
      /**
       * @brief Default constructor.
       *
       * Initializes the encapsulated resource using its default constructor
       * and marks it as not being written.
       */
      constexpr
      Unsafe()
        : m_writing(false)
      {}

      /**
       * @brief Constructs an Unsafe wrapper by copying a resource.
       *
       * @param resource The resource to be encapsulated.
       */
      constexpr
      Unsafe(const Resource& resource)
        : m_resource(resource),
          m_writing(false)
      {}

      /**
       * @brief Constructs an Unsafe wrapper by moving a resource.
       *
       * @param resource The resource to be encapsulated.
       */
      constexpr
      Unsafe(Resource&& resource)
        : m_resource(std::move(resource)),
          m_writing(false)
      {}

      /**
       * @brief Move constructor.
       *
       * Transfers ownership of the resource from another Unsafe object.
       *
       * @param other The Unsafe object to move from.
       */
      constexpr
      Unsafe(Unsafe&& other)
        : m_resource(std::move(other.m_resource)),
          m_writing(false)
      {}

      /**
       * @brief Copy constructor.
       *
       * Copies the encapsulated resource from another Unsafe object.
       *
       * @param other The Unsafe object to copy from.
       */
      constexpr
      Unsafe(const Unsafe& other)
        : m_resource(other.m_resource),
          m_writing(false)
      {}
#else
      /**
       * @brief Default constructor.
       *
       * Initializes the encapsulated resource using its default constructor.
       */
      constexpr
      Unsafe() = default;

      /**
       * @brief Constructs an Unsafe wrapper from a const reference to a
       * resource.
       *
       * @param resource The resource to be encapsulated.
       */
      constexpr
      Unsafe(const Resource& resource)
        : m_resource(resource)
      {}

      /**
       * @brief Constructs an Unsafe wrapper by moving a resource.
       *
       * @param resource The resource to be encapsulated.
       */
      constexpr
      Unsafe(Resource&& resource)
        : m_resource(std::move(resource))
      {}

      /**
       * @brief Move constructor.
       *
       * Transfers ownership of the resource from another Unsafe object.
       *
       * @param other The Unsafe object to move from.
       */
      constexpr
      Unsafe(Unsafe&& other)
        : m_resource(std::move(other.m_resource))
      {}

      /**
       * @brief Copy constructor.
       *
       * Copies the encapsulated resource from another Unsafe object.
       *
       * @param other The Unsafe object to copy from.
       */
      constexpr
      Unsafe(const Unsafe& other)
        : m_resource(other.m_resource)
      {}
#endif

      /**
       * @brief Assigns a new value to the encapsulated resource.
       *
       * This operator assigns a new value to the resource. In thread-safe
       * debug builds, the operation is protected by a lock to ensure that no
       * concurrent write occurs.
       *
       * @param resource The new value for the resource.
       * @return A reference to this Unsafe object.
       */
      constexpr
      Unsafe& operator=(const Resource& resource)
      {
#if defined(RODIN_THREAD_SAFE) && !defined(NDEBUG)
        lock();
        m_resource = resource;
        unlock();
#else
        m_resource = resource;
#endif
        return *this;
      }

      /**
       * @brief Assigns a new value to the encapsulated resource using move semantics.
       *
       * This operator moves a new value into the resource. In thread-safe
       * debug builds, the operation is protected by a lock to ensure exclusive
       * access.
       *
       * @param resource The new value for the resource.
       * @return A reference to this Unsafe object.
       */
      constexpr
      Unsafe& operator=(Resource&& resource)
      {
#if defined(RODIN_THREAD_SAFE) && !defined(NDEBUG)
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
       * Assigns the value of another Unsafe object to this one.
       *
       * @param other The Unsafe object to copy from.
       * @return A reference to this Unsafe object.
       */
      constexpr
      Unsafe& operator=(const Unsafe& other)
      {
        if (this != &other)
        {
#if defined(RODIN_THREAD_SAFE) && !defined(NDEBUG)
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
       * Transfers the resource from another Unsafe object into this one.
       *
       * @param other The Unsafe object to move from.
       * @return A reference to this Unsafe object.
       */
      constexpr
      Unsafe& operator=(Unsafe&& other)
      {
#if defined(RODIN_THREAD_SAFE) && !defined(NDEBUG)
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
       * @brief Modifies the encapsulated resource.
       *
       * Accepts a callable that takes a non-const reference to the resource
       * and applies modifications. In thread-safe debug builds, the operation
       * is protected by a lock to ensure exclusive access.
       *
       * @tparam F A callable type that can be invoked with a reference to Resource.
       * @param f The callable that modifies the resource.
       * @return A reference to this Unsafe object.
       */
      template <class F>
      constexpr
      Unsafe& write(F&& f)
      {
        static_assert(std::is_invocable_v<F, Resource&>);
#if defined(RODIN_THREAD_SAFE) && !defined(NDEBUG)
        lock();
        f(m_resource);
        unlock();
#else
        f(m_resource);
#endif
        return *this;
      }

#if defined(RODIN_THREAD_SAFE) && !defined(NDEBUG)
    protected:
      /**
       * @brief Locks the resource for writing.
       *
       * Checks if the resource is already being written to. If so, throws an
       * UnsafeAccessException.
       */
      constexpr
      void lock()
      {
        if (m_writing)
        {
          UnsafeAccessException(*this, __func__)
            << "Multiple threads writing to the resource at the same time."
            << Alert::Raise;
        }
        m_writing = true;
      }

      /**
       * @brief Unlocks the resource after writing.
       */
      constexpr
      void unlock()
      {
        m_writing = false;
      }
#endif

    private:
      /// The encapsulated resource.
      Resource m_resource;
#if defined(RODIN_THREAD_SAFE) && !defined(NDEBUG)
      /// Atomic flag indicating if the resource is currently being written to.
      std::atomic_bool m_writing;
#endif
  };
}

#endif

