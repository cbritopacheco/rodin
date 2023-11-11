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

    protected:
      inline
      constexpr
      Mutable& lock()
      {
        m_lock.lock();
        return *this;
      }

      inline
      constexpr
      Mutable& unlock()
      {
        m_lock.unlock();
        return *this;
      }

    private:
      Lock m_lock;
      Resource m_resource;
  };
}

#endif

