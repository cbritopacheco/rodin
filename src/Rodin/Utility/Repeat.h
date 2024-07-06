#ifndef RODIN_UTILITY_REPEAT_H
#define RODIN_UTILITY_REPEAT_H

#include <utility>

#include "Rodin/Tuple.h"

namespace Rodin::Utility
{
  template <size_t N, class T>
  struct Repeat;

  template <class T>
  struct Repeat<0, T>
  {
    using Type = Tuple<>;
  };

  template <class T>
  struct Repeat<1, T>
  {
    using Type = Tuple<T>;
  };

  template <size_t N, class T>
  struct Repeat
  {
    using Type =
      decltype(
        std::declval<Tuple<T>>().concatenate(std::declval<typename Repeat<N - 1, T>::Type>()));
  };
}

#endif
