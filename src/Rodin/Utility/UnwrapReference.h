/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_UNWRAPREFERENCE_H
#define RODIN_UTILITY_UNWRAPREFERENCE_H

/**
 * @file
 * @brief Defines utilities for unwrapping reference_wrapper types.
 */

#include <functional>

namespace Rodin::Utility
{
  /**
   * @brief Metafunction to unwrap std::reference_wrapper types.
   * @ingroup UtilityModule
   * @tparam T The type to potentially unwrap.
   *
   * For non-reference_wrapper types, the result is T itself.
   * For std::reference_wrapper<U>, the result is U&.
   *
   * Example usage:
   * @code{.cpp}
   * using Unwrapped1 = UnwrapReference<int>::Type;                         // int
   * using Unwrapped2 = UnwrapReference<std::reference_wrapper<int>>::Type; // int&
   * @endcode
   */
  template<class T>
  struct UnwrapReference
  {
    using Type = T;  ///< Identity for non-reference_wrapper types
  };

  /**
   * @brief Specialization for std::reference_wrapper.
   * @ingroup UtilityModule
   * @tparam U The type wrapped by reference_wrapper.
   */
  template <class U>
  struct UnwrapReference<std::reference_wrapper<U>>
  {
    using Type = U&;  ///< Reference to the wrapped type
  };

  /**
   * @brief Combines decay and reference unwrapping.
   * @ingroup UtilityModule
   * @tparam T The type to decay and unwrap.
   *
   * First applies std::decay_t to T, then unwraps any reference_wrapper.
   * This is useful for normalizing types in template contexts.
   */
  template <class T>
  struct UnwrapRefDecay : UnwrapReference<std::decay_t<T>> {};
}

#endif
