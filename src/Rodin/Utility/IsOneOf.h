#ifndef RODIN_UTILITY_ISONEOF_H
#define RODIN_UTILITY_ISONEOF_H

#include <type_traits>

namespace Rodin::Utility
{
  template <typename ...>
  struct IsOneOf
  {
    static constexpr bool Value = false;
  };

  template <typename F, typename S, typename ... T>
  struct IsOneOf<F, S, T...> {
    static constexpr bool Value =
      std::is_same<F, S>::Value || IsOneOf<F, T...>::Value;
  };
}

#endif
