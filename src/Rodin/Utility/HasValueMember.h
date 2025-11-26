/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_HASVALUEMEMBER_H
#define RODIN_UTILITY_HASVALUEMEMBER_H

/**
 * @file
 * @brief Defines the HasValueMember type trait for detecting nested Value members.
 */

#include <type_traits>

namespace Rodin::Utility
{
  /**
   * @brief Type trait to detect if a type has a nested Value member type.
   * @ingroup UtilityModule
   * @tparam T The type to check.
   *
   * This trait uses SFINAE to detect whether a type T has a nested
   * member type named `Value`. The primary template handles the case
   * where no such member exists.
   *
   * Example usage:
   * @code{.cpp}
   * struct WithValue { using Value = int; };
   * struct WithoutValue { };
   *
   * static_assert(HasValueMember<WithValue>::Value == true);
   * static_assert(HasValueMember<WithoutValue>::Value == false);
   * @endcode
   */
  template <typename T, typename = void>
  struct HasValueMember
  {
    static constexpr const bool Value = false;  ///< False when T has no Value member
  };

  /**
   * @brief Specialization for types that have a nested Value member type.
   * @ingroup UtilityModule
   * @tparam T The type being checked.
   */
  template <class T>
  struct HasValueMember<T, std::void_t<typename T::Value>>
  {
    static constexpr const bool Value = true;  ///< True when T has a Value member type
  };
}

#endif
