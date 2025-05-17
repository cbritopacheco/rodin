#include "ThreadPool.h"

namespace Rodin::Threads
{
  namespace Private
  {
    size_t getRodinNumThreads()
    {
      const char* env = std::getenv("RODIN_NUM_THREADS");
      if (env && *env)
        return std::stoi(env);
      else
        return 1;
    }
  }
}
