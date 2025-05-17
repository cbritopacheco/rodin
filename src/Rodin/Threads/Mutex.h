#ifndef RODIN_MUTEX_H
#define RODIN_MUTEX_H

#include <mutex>

#include "Rodin/Configure.h"

namespace Rodin::Threads
{
  /**
   * @brief Alias for std::mutex.
   *
   * This alias defines Mutex as an alternative name for std::mutex, which is
   * used throughout Rodin for mutual exclusion in multithreaded contexts.
   */
  using Mutex = std::mutex;
}

#endif
