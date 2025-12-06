/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_DEPENDENTVALUE_H
#define RODIN_UTILITY_DEPENDENTVALUE_H

/**
 * @file
 * @brief Defines the DependentValue utility for type-dependent compile-time constants.
 */

namespace Rodin::Utility
{
  /**
   * @brief Creates a type-dependent compile-time constant value.
   * @ingroup UtilityModule
   * @tparam T The type of the value.
   * @tparam Value_ The compile-time constant value.
   * @tparam Args Dependent type parameters (unused, but create type dependency).
   *
   * DependentValue provides a way to create compile-time constants that
   * depend on template parameters, preventing premature evaluation during
   * template parsing. This is useful in conditional compilation within
   * template contexts.
   *
   * Example usage:
   * @code{.cpp}
   * template <typename T>
   * void foo()
   * {
   *     if constexpr (DependentValue<bool, false, T>::Value)
   *     {
   *         // This branch won't be evaluated unless explicitly instantiated
   *     }
   * }
   * @endcode
   */
  template <class T, T Value_, class ... Args>
  struct DependentValue
  {
    static constexpr bool Value = Value_;  ///< The dependent value
  };
}

#endif
