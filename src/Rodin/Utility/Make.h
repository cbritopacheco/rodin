#ifndef RODIN_UTILITY_MAKE_H
#define RODIN_UTILITY_MAKE_H

#include <functional>

namespace Rodin::Utility
{
  template <class T>
  struct Make
  {
    template <class ... Params>
    T operator()(Params&&... params)
    {
      return T(std::forward<Params>(params)...);
    }
  };
}

#endif
