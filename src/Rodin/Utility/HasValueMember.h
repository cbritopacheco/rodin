#ifndef RODIN_UTILITY_HASVALUEMEMBER_H
#define RODIN_UTILITY_HASVALUEMEMBER_H

#include <type_traits>

namespace Rodin::Utility
{
  template <typename T, typename = void>
  struct HasValueMember
  {
    static constexpr const bool Value = false;
  };

  template <class T>
  struct HasValueMember<T, std::void_t<typename T::Value>>
  {
    static constexpr const bool Value = true;
  };
}

#endif
