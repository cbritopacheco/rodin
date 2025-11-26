/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_FALSE_H
#define RODIN_UTILITY_FALSE_H

/**
 * @file
 * @brief Defines the False compile-time constant for static_assert in template contexts.
 */

namespace Rodin::Utility
{
  /**
   * @brief Dependent false value for use in static_assert within template contexts.
   * @ingroup UtilityModule
   * @tparam ... Variadic template parameters (unused, but create type dependency).
   *
   * This utility provides a way to create a `static_assert(false, ...)` that only
   * triggers when a template is actually instantiated, rather than when it is
   * parsed. This is essential for implementing proper error messages in template
   * specializations that should never be instantiated.
   *
   * Example usage:
   * @code{.cpp}
   * template <typename T>
   * struct MyStruct
   * {
   *     static_assert(Utility::False<T>, "This template should not be instantiated");
   * };
   *
   * template <>
   * struct MyStruct<int>
   * {
   *     // Valid specialization
   * };
   * @endcode
   */
  template <class...>
  inline constexpr bool False = false;
}

#endif