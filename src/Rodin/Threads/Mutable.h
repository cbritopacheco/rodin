#ifndef RODIN_MUTABLE_H
#define RODIN_MUTABLE_H

#include "Rodin/Configure.h"

#include "Mutex.h"

namespace Rodin::Threads
{
  template <class Resource, class Lock = Mutex>
  class Mutable
  {
    public:
      constexpr
      Mutable()
      {}

      constexpr
      Mutable(const Resource& resource)
        : m_resource(resource)
      {}

      constexpr
      Mutable(Resource&& resource)
        : m_resource(std::move(resource))
      {}

      constexpr
      Mutable(Mutable&& other)
        : m_resource(std::move(other.m_resource))
      {}

      constexpr
      Mutable(const Mutable& other)
        : m_resource(other.m_resource)
      {}

      constexpr
      Mutable& operator=(const Resource& resource)
      {
#ifdef RODIN_THREAD_SAFE
        lock();
        m_resource = resource;
        unlock();
#endif
        return *this;
      }

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

      inline
      constexpr
      const Resource& read() const
      {
        return m_resource;
      }

      template <class F>
      inline
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
      inline
      constexpr
      void lock()
      {
        m_lock.lock();
      }

      inline
      constexpr
      void unlock()
      {
        m_lock.unlock();
      }
#endif

    private:
#ifdef RODIN_THREAD_SAFE
      Lock m_lock;
#endif
      Resource m_resource;
  };
}

#endif

