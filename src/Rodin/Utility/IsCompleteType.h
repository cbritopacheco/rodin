/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_ISCOMPLETETYPE_H
#define RODIN_UTILITY_ISCOMPLETETYPE_H

/**
 * @file
 * @brief Defines the IsCompleteType type trait for detecting complete types.
 */

#include <cstdlib>
#include <type_traits>

namespace Rodin::Utility
{
  namespace Internal
  {
    /**
     * @brief Helper function for IsCompleteType detection (complete type case).
     * @ingroup UtilityModule
     *
     * This overload is selected when sizeof(T) is valid, indicating T is complete.
     */
    template <class T, std::size_t = sizeof(T)>
    std::true_type IsCompleteTypeImpl(T *);

    /**
     * @brief Helper function for IsCompleteType detection (incomplete type case).
     * @ingroup UtilityModule
     *
     * This fallback overload is selected when sizeof(T) is not valid.
     */
    std::false_type IsCompleteTypeImpl(...);
  }

  /**
   * @brief Type trait to detect if a type is complete (fully defined).
   * @ingroup UtilityModule
   * @tparam T The type to check.
   *
   * A type is complete if its size is known (i.e., it has been fully defined).
   * Incomplete types include forward declarations and array types of unknown bound.
   *
   * This trait is useful for SFINAE-based conditional compilation when
   * different behavior is needed for complete vs. incomplete types.
   *
   * Example usage:
   * @code{.cpp}
   * struct Complete { int x; };
   * struct Incomplete;
   *
   * static_assert(IsCompleteType<Complete>::value == true);
   * static_assert(IsCompleteType<Incomplete>::value == false);
   * @endcode
   */
  template <class T>
  using IsCompleteType = decltype(Internal::IsCompleteTypeImpl(std::declval<T*>()));

}

#endif
