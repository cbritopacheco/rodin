/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_PRODUCT_H
#define RODIN_UTILITY_PRODUCT_H

/**
 * @file
 * @brief Defines the Product metafunction for computing Cartesian products of tuple types.
 */

#include "Rodin/Tuple.h"

namespace Rodin::Utility
{
  /**
   * @brief Metafunction to compute the Cartesian product of two Tuple type-lists.
   * @ingroup UtilityModule
   *
   * Given two `Tuple<...>` types, `A = Tuple<A1, A2, ...>` and `B = Tuple<B1, B2, ...>`,
   * `Product<A, B>::Type<Pair>` will be `Tuple< Pair<Ai,Bj>... >` for all i,j.
   *
   * @tparam ... Primary template (undefined). Specializations below implement the logic.
   */
  template <class ...>
  struct Product;

  /**
   * @brief Specialization for Tuple with at least two elements (H1, H2, Hs…).
   * @ingroup UtilityModule
   *
   * For each G in the second tuple, creates Pair<H1,G> and Pair<H2,G>,
   * then concatenates that with the product of the remainder Hs… with Gs….
   *
   * @tparam H1 First head type of the first Tuple.
   * @tparam H2 Second head type of the first Tuple.
   * @tparam Hs Remaining types of the first Tuple.
   * @tparam Gs All types of the second Tuple.
   *
   * @note This relies on `Tuple<...>::concatenate(Tuple<...>)` to stitch the pieces together.
   */
  template <class H1, class H2, class ... Hs, class ... Gs>
  struct Product<Tuple<H1, H2, Hs...>, Tuple<Gs...>>
  {
    /**
     * @brief Resulting Cartesian product as a Tuple of Pair<..., ...>.
     * @tparam Pair Binary template taking one type from each input Tuple.
     *
     * Example:
     * @code{.cpp}
     * using A = Tuple<int, double, char>;
     * using B = Tuple<bool, long>;
     * using P = Product<A,B>::Type<std::pair>;
     * // yields Tuple<
     * //   std::pair<int,bool>, std::pair<int,long>,
     * //   std::pair<double,bool>, std::pair<double,long>,
     * //   std::pair<char,bool>, std::pair<char,long>
     * // >
     * @endcode
     */
    template <template <class, class> class Pair>
    using Type =
      decltype(
        std::declval<Tuple<Pair<H1, Gs>...>>().concatenate(
          std::declval<Tuple<Pair<H2, Gs>...>>()).concatenate(
            std::declval<typename Product<Tuple<Hs...>, Tuple<Gs...>>::template Type<Pair>>()));
  };

  /**
   * @brief Specialization for Tuple with exactly one element.
   * @ingroup UtilityModule
   *
   * Generates a single Tuple< Pair<H, G1>, Pair<H, G2>, ...>.
   *
   * @tparam H The lone type in the first Tuple.
   * @tparam Gs All types of the second Tuple.
   */
  template <class H, class ... Gs>
  struct Product<Tuple<H>, Tuple<Gs...>>
  {
    template <template <class, class> class Pair>
    using Type = Tuple<Pair<H, Gs>...>;
  };

  /**
   * @brief Specialization for empty first Tuple.
   * @ingroup UtilityModule
   * @tparam Gs All types of the second Tuple.
   */
  template <class ... Gs>
  struct Product<Tuple<>, Tuple<Gs...>>
  {
    template <template <class, class> class Pair>
    using Type = Tuple<>;
  };

  /**
   * @brief Specialization for empty second Tuple.
   * @ingroup UtilityModule
   * @tparam Gs All types of the first Tuple.
   */
  template <class ... Gs>
  struct Product<Tuple<Gs...>, Tuple<>>
  {
    template <template <class, class> class Pair>
    using Type = Tuple<>;
  };
}

#endif
