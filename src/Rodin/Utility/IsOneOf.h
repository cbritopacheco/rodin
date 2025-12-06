/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_ISONEOF_H
#define RODIN_UTILITY_ISONEOF_H

/**
 * @file
 * @brief Defines the IsOneOf type trait for checking type membership in a list.
 */

#include <type_traits>

namespace Rodin::Utility
{
  /**
   * @brief Type trait to check if a type is one of a specified list of types.
   * @ingroup UtilityModule
   *
   * The primary template handles the base case of an empty type list,
   * for which the result is always false.
   *
   * Example usage:
   * @code{.cpp}
   * static_assert(IsOneOf<int, float, int, double>::Value == true);
   * static_assert(IsOneOf<char, float, int, double>::Value == false);
   * @endcode
   */
  template <typename ...>
  struct IsOneOf
  {
    static constexpr bool Value = false;  ///< False for empty type list
  };

  /**
   * @brief Specialization for non-empty type lists.
   * @ingroup UtilityModule
   * @tparam F The type to search for.
   * @tparam S The first type in the remaining list.
   * @tparam T The rest of the types in the list.
   *
   * Recursively checks if F matches S, or if F is one of the remaining types T...
   */
  template <typename F, typename S, typename ... T>
  struct IsOneOf<F, S, T...>
  {
    /// True if F matches S or is one of T...
    static constexpr bool Value =
      std::is_same<F, S>::value || IsOneOf<F, T...>::Value;
  };
}

#endif
