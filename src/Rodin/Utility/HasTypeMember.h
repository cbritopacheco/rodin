#ifndef RODIN_UTILITY_HASTYPEMEMBER_H
#define RODIN_UTILITY_HASTYPEMEMBER_H

#include <utility>
#include <type_traits>

namespace Rodin::Utility
{
  template <typename T, typename = void>
  struct HasTypeMember
  {
    static constexpr const bool Value = false;
  };

  template <class T>
  struct HasTypeMember<T, std::void_t<typename T::Type>>
  {
    static constexpr const bool Value = true;
  };
}


#endif
