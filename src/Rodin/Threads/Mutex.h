#ifndef RODIN_MUTEX_H
#define RODIN_MUTEX_H

#include <mutex>

#include "Rodin/Configure.h"

namespace Rodin::Threads
{
  using Mutex = std::mutex;

  template <class Mutex>
  class LockGuard : public std::lock_guard<Mutex>
  {
    public:
      using Parent = std::lock_guard<Mutex>;
      using Parent::Parent;
  };

  template <class Mutex>
  LockGuard(Mutex&) -> LockGuard<Mutex>;

  template <class Lock>
  class GuardCriticalSection
  {
    public:
      template <class F>
      GuardCriticalSection(Lock& lock, F&& f)
      {
        lock.lock();
        f();
        lock.unlock();
      }
  };

  template <class Resource, class Lock = Mutex>
  class Mutable
  {
    public:
      template <class ... Args>
      Mutable(Args&&... args)
        : m_resource(std::forward<Args>(args)...)
      {}

      inline
      constexpr
      const Resource& read() const
      {
        return m_resource.get();
      }

      template <class F>
      inline
      constexpr
      Mutable& write(F&& f)
      {
        lock();
        f(m_resource);
        unlock();
        return *this;
      }

      inline
      constexpr
      Mutable& lock() const
      {
        m_lock.lock();
        return *this;
      }

      inline
      constexpr
      Mutable& unlock() const
      {
        m_lock.unlock();
        return *this;
      }

    private:
      mutable Lock m_lock;
      mutable Resource m_resource;
  };
}

#endif
