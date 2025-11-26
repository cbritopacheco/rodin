/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_FORCONSTEXPR_H
#define RODIN_UTILITY_FORCONSTEXPR_H

/**
 * @file
 * @brief Defines compile-time loop utilities for template metaprogramming.
 *
 * This file provides utilities for performing compile-time iterations
 * over parameter packs and index sequences, enabling powerful template
 * metaprogramming patterns.
 */

#include <utility>

namespace Rodin::Utility
{
  /**
   * @brief Performs the application of a function over all arguments.
   * @ingroup UtilityModule
   * @tparam F The callable type.
   * @tparam Args The argument types.
   * @param f The function to apply.
   * @param args The arguments to apply the function to.
   *
   * This function uses fold expressions to apply the callable @p f
   * to each argument in @p args in sequence. Useful for compile-time
   * iteration over heterogeneous argument packs.
   *
   * Example usage:
   * @code{.cpp}
   * auto printer = [](auto x) { std::cout << x << " "; };
   * Utility::For(printer, 1, 2.5, "hello");
   * // Output: 1 2.5 hello
   * @endcode
   */
  template <class F, class... Args>
  constexpr void For(F&& f, Args&&... args)
  {
    (f(std::forward<Args>(args)), ...);
  }

  namespace Internal
  {
    /**
     * @brief Compile-time index wrapper.
     * @ingroup UtilityModule
     * @tparam N The index value.
     *
     * Wraps an index value as a type with a static constexpr member,
     * allowing compile-time access to loop indices within constexpr contexts.
     */
    template <size_t N>
    struct Index
    {
      static constexpr const size_t value = N;  ///< The index value
      constexpr operator size_t() const { return N; }  ///< Implicit conversion to size_t
    };

    template <size_t ... Is, class F>
    constexpr void ForIndexImpl(F&& f, std::index_sequence<Is...>)
    {
      (std::forward<F>(f)(Index<Is>{}), ...);
    }
  }

  /**
   * @brief Executes a function N times with compile-time index access.
   * @ingroup UtilityModule
   * @tparam N Number of times to execute the function.
   * @tparam F Callable type.
   * @param f Function to execute; receives an object with a `v.value` member for each iteration.
   *
   * This utility enables compile-time loops where the loop index is
   * available as a compile-time constant.
   *
   * Example usage:
   * @code{.cpp}
   * Utility::ForIndex<5>(
   *   [&](auto i){ std::cout << i.value << std::endl; });
   * // Output: 0 1 2 3 4
   * @endcode
   */
  template <size_t N, class F>
  constexpr void ForIndex(F&& f)
  {
    Internal::ForIndexImpl(std::forward<F>(f), std::make_index_sequence<N>{});
  }
}

#endif
