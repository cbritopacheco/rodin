/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_CORE_COMMON_H
#define RODIN_CORE_COMMON_H

#include <cmath>
#include <type_traits>

namespace Rodin
{
  template <class T>
  constexpr
  std::enable_if_t<!std::is_void_v<decltype(std::abs(std::declval<T>()))>, T>
  abs(const T& x)
  {
    return std::abs(x);
  }

  template <class T>
  constexpr
  std::enable_if_t<!std::is_void_v<decltype(std::sqrt(std::declval<T>()))>, T>
  sqrt(const T& x)
  {
    return std::sqrt(x);
  }

  template <class T>
  constexpr
  std::enable_if_t<!std::is_void_v<decltype(std::isnan(std::declval<T>()))>, bool>
  isNaN(const T& x)
  {
    return std::isnan(x);
  }

  template <class T>
  constexpr
  std::enable_if_t<!std::is_void_v<decltype(std::isinf(std::declval<T>()))>, bool>
  isInf(const T& x)
  {
    return std::isinf(x);
  }
}

#endif
