/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_ZIP_H
#define RODIN_UTILITY_ZIP_H

#include "Rodin/Tuple.h"

namespace Rodin::Utility
{
  template <class ...>
  struct Zip;

  template <class ... Ts, class ... Gs>
  struct Zip<Tuple<Ts...>, Tuple<Gs...>>
  {
    static_assert(sizeof...(Ts) == sizeof...(Gs));

    template <template <class, class> class Pair>
    using Type = Tuple<Pair<Ts, Gs>...>;
  };
}

#endif



