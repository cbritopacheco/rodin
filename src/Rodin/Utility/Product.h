/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_PRODUCT_H
#define RODIN_UTILITY_PRODUCT_H

#include "Rodin/Tuple.h"

namespace Rodin::Utility
{
  template <class ...>
  struct Product;

  template <class H1, class H2, class ... Hs, class ... Gs>
  struct Product<Tuple<H1, H2, Hs...>, Tuple<Gs...>>
  {
    template <template <class, class> class Pair>
    using Type =
      decltype(
        std::declval<Tuple<Pair<H1, Gs>...>>().concatenate(
          std::declval<Tuple<Pair<H2, Gs>...>>()).concatenate(
            std::declval<typename Product<Tuple<Hs...>, Tuple<Gs...>>::template Type<Pair>>()));
  };

  template <class H, class ... Gs>
  struct Product<Tuple<H>, Tuple<Gs...>>
  {
    template <template <class, class> class Pair>
    using Type = Tuple<Pair<H, Gs>...>;
  };

  template <class ... Gs>
  struct Product<Tuple<>, Tuple<Gs...>>
  {
    template <template <class, class> class Pair>
    using Type = Tuple<>;
  };

  template <class ... Gs>
  struct Product<Tuple<Gs...>, Tuple<>>
  {
    template <template <class, class> class Pair>
    using Type = Tuple<>;
  };
}

#endif


