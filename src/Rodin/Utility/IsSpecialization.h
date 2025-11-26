/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_ISSPECIALIZATION_H
#define RODIN_UTILITY_ISSPECIALIZATION_H

/**
 * @file
 * @brief Defines the IsSpecialization type trait for detecting template specializations.
 */

#include <type_traits>

namespace Rodin::Utility
{
  /**
   * @brief Type trait to check if a type is a specialization of a given template.
   * @ingroup UtilityModule
   * @tparam Test The type to test.
   * @tparam Ref The reference template to check against.
   *
   * The primary template handles the case where Test is not a specialization
   * of the Ref template.
   *
   * Example usage:
   * @code{.cpp}
   * static_assert(IsSpecialization<std::vector<int>, std::vector>::Value == true);
   * static_assert(IsSpecialization<int, std::vector>::Value == false);
   * @endcode
   */
  template <class Test, template<class...> class Ref>
  struct IsSpecialization
  {
    static constexpr const bool Value = false;  ///< False when Test is not a specialization of Ref
  };

  /**
   * @brief Specialization for types that are specializations of the reference template.
   * @ingroup UtilityModule
   * @tparam Ref The reference template.
   * @tparam Args The template arguments used to specialize Ref.
   */
  template<template<class...> class Ref, class... Args>
  struct IsSpecialization<Ref<Args...>, Ref>
  {
    static constexpr const bool Value = true;  ///< True when Test is Ref<Args...>
  };
}

#endif
