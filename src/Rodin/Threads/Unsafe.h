#ifndef RODIN_UNSAFE_H
#define RODIN_UNSAFE_H

#include "Rodin/Configure.h"
#include "Rodin/Alert/MemberFunctionException.h"

namespace Rodin::Threads
{
  template <class T, class FuncName>
  class UnsafeAccessException : public Alert::MemberFunctionException<T, FuncName>
  {
    public:
      using Parent = Alert::MemberFunctionException<T, FuncName>;

      UnsafeAccessException(const T& cls, const FuncName& funcName)
        : Parent(cls, funcName)
      {
        *this << Alert::Text(
            "Rodin::Threads::UnsafeAccessException").setBold().setUnderline() << ": ";
      };
  };

  template <class Resource>
  class Unsafe
  {
    public:
#if defined(RODIN_THREAD_SAFE) && !defined(NDEBUG)
      constexpr
      Unsafe()
        : m_writing(false)
      {}

      constexpr
      Unsafe(const Resource& resource)
        : m_resource(resource),
          m_writing(false)
      {}

      constexpr
      Unsafe(Resource&& resource)
        : m_resource(std::move(resource)),
          m_writing(false)
      {}

      constexpr
      Unsafe(Unsafe&& other)
        : m_resource(std::move(other.m_resource)),
          m_writing(false)
      {}

      constexpr
      Unsafe(const Unsafe& other)
        : m_resource(other.m_resource),
          m_writing(false)
      {}
#else
      constexpr
      Unsafe() = default;

      constexpr
      Unsafe(const Resource& resource)
        : m_resource(resource)
      {}

      constexpr
      Unsafe(Resource&& resource)
        : m_resource(std::move(resource))
      {}

      constexpr
      Unsafe(Unsafe&& other)
        : m_resource(std::move(other.m_resource))
      {}

      constexpr
      Unsafe(const Unsafe& other)
        : m_resource(other.m_resource)
      {}
#endif

      inline
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

      inline
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

      inline
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

      inline
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

      inline
      constexpr
      const Resource& read() const
      {
        return m_resource;
      }

      template <class F>
      inline
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
      inline
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

      inline
      constexpr
      void unlock()
      {
        m_writing = false;
      }
#endif

    private:
      Resource m_resource;
#if defined(RODIN_THREAD_SAFE) && !defined(NDEBUG)
      std::atomic_bool m_writing;
#endif
  };
}

#endif

