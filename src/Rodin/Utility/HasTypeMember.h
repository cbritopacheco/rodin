/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_HASTYPEMEMBER_H
#define RODIN_UTILITY_HASTYPEMEMBER_H

/**
 * @file
 * @brief Defines the HasTypeMember type trait for detecting nested Type members.
 */

#include <utility>
#include <type_traits>

namespace Rodin::Utility
{
  /**
   * @brief Type trait to detect if a type has a nested Type member.
   * @ingroup UtilityModule
   * @tparam T The type to check.
   *
   * This trait uses SFINAE to detect whether a type T has a nested
   * member type named `Type`. The primary template handles the case
   * where no such member exists.
   *
   * Example usage:
   * @code{.cpp}
   * struct WithType { using Type = int; };
   * struct WithoutType { };
   *
   * static_assert(HasTypeMember<WithType>::Value == true);
   * static_assert(HasTypeMember<WithoutType>::Value == false);
   * @endcode
   */
  template <typename T, typename = void>
  struct HasTypeMember
  {
    static constexpr const bool Value = false;  ///< False when T has no Type member
  };

  /**
   * @brief Specialization for types that have a nested Type member.
   * @ingroup UtilityModule
   * @tparam T The type being checked.
   */
  template <class T>
  struct HasTypeMember<T, std::void_t<typename T::Type>>
  {
    static constexpr const bool Value = true;  ///< True when T has a Type member
  };
}

#endif
