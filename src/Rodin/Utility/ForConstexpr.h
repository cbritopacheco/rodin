#ifndef RODIN_UTILITY_FORCONSTEXPR_H
#define RODIN_UTILITY_FORCONSTEXPR_H

#include <functional>
#include <type_traits>


namespace Rodin::Utility::Internal
{
}

namespace Rodin::Utility
{
  template <class F, class... Args>
  constexpr void For(F&& f, Args&&... args)
  {
    (f(std::forward<Args>(args)), ...);
  }

  namespace Internal
  {
    template <size_t N>
    struct Index
    {
      static constexpr const size_t value = N;
      constexpr operator size_t() const { return N; }
    };

    template <size_t ... Is, class F>
    constexpr void ForIndexImpl(F&& f, std::index_sequence<Is...>)
    {
      (std::forward<F>(f)(Index<Is>{}), ...);
    }
  }

  template <size_t N, class F>
  constexpr void ForIndex(F&& f)
  {
    Internal::ForIndexImpl(std::forward<F>(f), std::make_index_sequence<N>{});
  }
}

#endif
