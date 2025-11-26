/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_WRAP_H
#define RODIN_UTILITY_WRAP_H

/**
 * @file
 * @brief Defines the Wrap metafunction for wrapping tuple element types.
 */

#include <utility>

#include "Rodin/Tuple.h"

namespace Rodin::Utility
{
  /**
   * @brief Metafunction to wrap all types in a Tuple with an external template.
   * @ingroup UtilityModule
   * @tparam ... Primary template (undefined). Specializations implement the logic.
   *
   * Wrap takes a Tuple type and an external template, and produces a new
   * Tuple where each element type has been wrapped with the external template.
   */
  template <class ...>
  class Wrap;

  /**
   * @brief Specialization for Tuple types.
   * @ingroup UtilityModule
   * @tparam Ts The types in the tuple.
   *
   * Example usage:
   * @code
   * using Original = Tuple<int, double, char>;
   * using Wrapped = Wrap<Original>::Type<std::optional>;
   * // Wrapped is Tuple<std::optional<int>, std::optional<double>, std::optional<char>>
   * @endcode
   */
  template <class ... Ts>
  class Wrap<Tuple<Ts...>>
  {
    public:
      /**
       * @brief The resulting tuple with wrapped types.
       * @tparam External The template to wrap each type with.
       */
      template <template <class> class External>
      using Type = Tuple<External<Ts>...>;
  };

}

#endif