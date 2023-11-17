#ifndef RODIN_SHARED_H
#define RODIN_SHARED_H

#include <shared_mutex>

#include "Rodin/Configure.h"

namespace Rodin::Threads
{
  template <class Resource>
  class Shared
  {
    public:
      constexpr
      Shared()
      {}

      constexpr
      Shared(const Resource& resource)
        : m_resource(resource)
      {}

      constexpr
      Shared(Resource&& resource)
        : m_resource(std::move(resource))
      {}

      constexpr
      Shared(Shared&& other)
        : m_resource(std::move(other.m_resource))
      {}

      constexpr
      Shared(const Shared& other)
        : m_resource(other.m_resource)
      {}

      constexpr
      Shared& operator=(const Resource& resource)
      {
#ifdef RODIN_THREAD_SAFE
        lock();
        m_resource = resource;
        unlock();
#endif
        return *this;
      }

      constexpr
      Shared& operator=(Resource&& resource)
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

      inline
      constexpr
      const Resource& read() const
      {
        return m_resource;
      }

      template <class F>
      inline
      constexpr
      const Shared& read(F&& f) const
      {
        static_assert(std::is_invocable_v<F, const Resource&>);
#ifdef RODIN_THREAD_SAFE
        m_lock.lock_shared();
        f(m_resource);
        m_lock.unlock_shared();
#else
        f(m_resource);
#endif
        return *this;
      }

      template <class F>
      inline
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
      mutable std::shared_mutex m_lock;
#endif
      Resource m_resource;
  };
}

#endif


