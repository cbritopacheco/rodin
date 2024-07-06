/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_EXTRACT_H
#define RODIN_UTILITY_EXTRACT_H

#include "Rodin/Tuple.h"

namespace Rodin::Utility
{
  template <class ...>
  struct Extract;

  template <class T, class ... Ts>
  struct Extract<Tuple<T, Ts...>>
  {
    template <template <class> class Extractor>
    using Type = Tuple<typename Extractor<T>::Type, typename Extractor<Ts>::Type...>;
  };
}

#endif

