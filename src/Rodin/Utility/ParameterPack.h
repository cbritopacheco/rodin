/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_PARAMETERPACK_H
#define RODIN_UTILITY_PARAMETERPACK_H

#include <type_traits>

namespace Rodin::Utility
{
  namespace Internal
  {
    template <std::size_t N, class ... Types>
    struct AtImpl;

    template <class First, class ... Rest>
    struct AtImpl<0, First, Rest...>
    {
      using Type = First;
    };

    template <std::size_t N, class First, class ... Rest>
    struct AtImpl<N, First, Rest...>
    {
      using Type = typename AtImpl<N - 1, Rest...>::Type;
    };

    template <template <class> class Predicate, class ...>
    struct AllImpl;

    template <template <class> class Predicate, class T1>
    struct AllImpl<Predicate, T1>
    {
      static constexpr bool Value = Predicate<T1>::Value;
    };

    template <template <class> class Predicate, class T1, class T2, class ... Ts>
    struct AllImpl<Predicate, T1, T2, Ts...>
    {
      static constexpr bool Value = Predicate<T1>::Value && AllImpl<Predicate, T2, Ts...>::Value;
    };
  }

  template <class ... Params>
  class ParameterPack
  {
    public:
      template <std::size_t Index>
      using At = typename Internal::AtImpl<Index, Params...>::Type;

      template <template <class> class Predicate>
      using All = typename Internal::AllImpl<Predicate, Params...>;

      static constexpr size_t Size = sizeof...(Params);
  };
}

#endif
