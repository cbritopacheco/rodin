/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_CORE_COMMON_HPP
#define RODIN_CORE_COMMON_HPP

#include <type_traits>
#include <cmath>

#include "Common.h"

namespace Rodin::Core
{
  template <class T>
  inline
  constexpr
  T abs(const T& x)
  {
    return x < T(0) ? -x : x;
  }

  template <class T>
  inline
  constexpr
  T sqrt(const T& x)
  {
    return std::sqrt(x);
  }
}

#endif
