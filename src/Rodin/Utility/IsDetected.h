#ifndef RODIN_UTILITY_ISDETECTED_H
#define RODIN_UTILITY_ISDETECTED_H

#include <type_traits>

namespace Rodin::Utility
{
  namespace Internal
  {
    template <class Default, class AlwaysVoid, template <class...> class Op, class... Args>
    struct Detector
    {
      using ValueType = std::false_type;
      using Type = Default;
    };

    template <class Default, template<class...> class Op, class... Args>
    struct Detector<Default, std::void_t<Op<Args...>>, Op, Args...>
    {
      using ValueType = std::true_type;
      using Type = Op<Args...>;
    };
  } // namespace detail

  struct NoneSuch
  {};

  template <template <class...> class Op, class... Args>
  using IsDetected = typename Internal::Detector<NoneSuch, void, Op, Args...>::ValueType;
}

#endif
