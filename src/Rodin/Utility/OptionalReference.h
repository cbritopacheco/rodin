/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_OPTIONALREFERENCE_H
#define RODIN_UTILITY_OPTIONALREFERENCE_H

/**
 * @file
 * @brief Defines the OptionalReference wrapper for optional reference semantics.
 */

#include <optional>
#include <functional>

namespace Rodin::Utility
{
  /**
   * @brief A wrapper providing optional reference semantics.
   * @ingroup UtilityModule
   * @tparam T The referenced type.
   *
   * OptionalReference provides a way to have an optional reference,
   * combining the functionality of std::optional with reference semantics.
   * This is useful when a function may or may not return a reference to an object.
   *
   * The class inherits from Optional<std::reference_wrapper<T>> and provides
   * convenient pointer-like and reference access operators.
   *
   * Example usage:
   * @code{.cpp}
   * int value = 42;
   * OptionalReference<int> optRef;  // Empty optional reference
   * optRef = value;  // Now refers to value
   * 
   * if (optRef)
   * {
   *     *optRef = 100;  // Modifies value through the reference
   * }
   * @endcode
   */
  template <typename T>
  class OptionalReference : public Optional<std::reference_wrapper<T>>
  {
    public:
      using Parent = Optional<std::reference_wrapper<T>>;  ///< Parent type alias
      using Parent::Parent;  ///< Inherit constructors

      /**
       * @brief Pointer access operator.
       * @return Pointer to the referenced object.
       */
      T* operator->()
      {
        return &(this->Parent::operator*().get());
      }

      /**
       * @brief Dereference operator.
       * @return Reference to the referenced object.
       */
      T& operator*()
      {
        return this->Parent::operator*().get();
      }
  };
}

#endif
