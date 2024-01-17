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
  }

  template <class ... Params>
  class ParameterPack
  {
    public:
      template <std::size_t Index>
      using At = typename Internal::AtImpl<Index, Params...>::Type;
  };
}

#endif
