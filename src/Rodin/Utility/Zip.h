/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_ZIP_H
#define RODIN_UTILITY_ZIP_H

/**
 * @file
 * @brief Defines the Zip metafunction for combining two tuples element-wise.
 */

#include "Rodin/Tuple.h"

namespace Rodin::Utility
{
  /**
   * @brief Metafunction to combine two Tuples element-wise using a binary template.
   * @ingroup UtilityModule
   *
   * Zip takes two Tuple types of equal size and a binary template Pair,
   * and produces a new Tuple where each element is `Pair<Ti, Gi>` for corresponding
   * elements Ti from the first tuple and Gi from the second.
   */
  template <class ...>
  struct Zip;

  /**
   * @brief Specialization for two Tuples.
   * @ingroup UtilityModule
   * @tparam Ts Types in the first tuple.
   * @tparam Gs Types in the second tuple.
   *
   * Both tuples must have the same number of elements (enforced by static_assert).
   *
   * Example usage:
   * @code{.cpp}
   * using A = Tuple<int, double, char>;
   * using B = Tuple<bool, float, short>;
   * using Zipped = Zip<A, B>::Type<std::pair>;
   * // Zipped is Tuple<std::pair<int,bool>, std::pair<double,float>, std::pair<char,short>>
   * @endcode
   */
  template <class ... Ts, class ... Gs>
  struct Zip<Tuple<Ts...>, Tuple<Gs...>>
  {
    static_assert(sizeof...(Ts) == sizeof...(Gs),
      "Zip requires tuples of equal size");

    /**
     * @brief The resulting tuple of paired types.
     * @tparam Pair Binary template to combine corresponding elements.
     */
    template <template <class, class> class Pair>
    using Type = Tuple<Pair<Ts, Gs>...>;
  };
}

#endif
