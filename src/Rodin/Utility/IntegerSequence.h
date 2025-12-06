/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_INTEGERSEQUENCE_H
#define RODIN_UTILITY_INTEGERSEQUENCE_H

/**
 * @file
 * @brief Defines the IntegerSequence class for compile-time integer sequences.
 */

namespace Rodin::Utility
{
  /**
   * @brief Represents a compile-time sequence of integers.
   * @ingroup UtilityModule
   * @tparam T The integer type (e.g., int, size_t).
   * @tparam Is The sequence of integer values.
   *
   * IntegerSequence is similar to std::integer_sequence and provides
   * a way to represent and manipulate sequences of integers at compile time.
   * This is useful for template metaprogramming operations like
   * tuple indexing and parameter pack expansion.
   *
   * Example usage:
   * @code{.cpp}
   * using Seq = IntegerSequence<int, 0, 1, 2, 3, 4>;
   * @endcode
   */
  template <typename T, T... Is>
  class IntegerSequence
  {
  };
}

#endif
