#ifndef RODIN_UTILITY_UNWRAPREFERENCE_H
#define RODIN_UTILITY_UNWRAPREFERENCE_H

#include <functional>

namespace Rodin::Utility
{
  template<class T>
  struct UnwrapReference { using Type = T; };

  template<class U>
  struct UnwrapReference<std::reference_wrapper<U>> { using Type = U&; };

  template<class T>
  struct UnwrapRefDecay : UnwrapReference<std::decay_t<T>> {};
}

#endif
