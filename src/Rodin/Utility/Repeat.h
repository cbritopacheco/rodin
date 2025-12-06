/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_REPEAT_H
#define RODIN_UTILITY_REPEAT_H

/**
 * @file
 * @brief Defines the Repeat metafunction for creating tuples with repeated types.
 */

#include <utility>

#include "Rodin/Tuple.h"

namespace Rodin::Utility
{
  /**
   * @brief Metafunction to create a Tuple with N repetitions of type T.
   * @ingroup UtilityModule
   * @tparam N The number of repetitions.
   * @tparam T The type to repeat.
   *
   * Repeat creates a Tuple containing N copies of type T. This is useful
   * for generating homogeneous tuple types at compile time.
   *
   * Example usage:
   * @code{.cpp}
   * using ThreeInts = Repeat<3, int>::Type;
   * // ThreeInts is Tuple<int, int, int>
   * @endcode
   */
  template <size_t N, class T>
  struct Repeat;

  /**
   * @brief Specialization for N = 0.
   * @ingroup UtilityModule
   * @tparam T The type (not used in result).
   *
   * Creates an empty Tuple.
   */
  template <class T>
  struct Repeat<0, T>
  {
    using Type = Tuple<>;  ///< Empty tuple for zero repetitions
  };

  /**
   * @brief Specialization for N = 1.
   * @ingroup UtilityModule
   * @tparam T The type to include.
   *
   * Creates a single-element Tuple containing T.
   */
  template <class T>
  struct Repeat<1, T>
  {
    using Type = Tuple<T>;  ///< Single-element tuple
  };

  /**
   * @brief General case for N > 1.
   * @ingroup UtilityModule
   * @tparam N The number of repetitions (N > 1).
   * @tparam T The type to repeat.
   *
   * Recursively concatenates Tuple<T> with Repeat<N-1, T>::Type.
   */
  template <size_t N, class T>
  struct Repeat
  {
    using Type =
      decltype(
        std::declval<Tuple<T>>().concatenate(std::declval<typename Repeat<N - 1, T>::Type>()));
  };
}

#endif
