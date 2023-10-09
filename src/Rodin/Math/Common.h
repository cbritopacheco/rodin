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

#include <boost/multi_array.hpp>

#include "Rodin/Types.h"

namespace Rodin::Math
{
  /**
  * @brief Computes the absolute value of a value of type T.
  * @param[in] x Value
  * @tparam T Type of value
  * @returns Absolute of value
  */
  template <class T>
  constexpr
  inline
  std::enable_if_t<!std::is_void_v<decltype(std::abs(std::declval<T>()))>, T>
  abs(const T& x)
  {
    return std::abs(x);
  }

  template <class Base, class Exponent>
  constexpr
  inline
  auto pow2(const Base& base)
  {
    return base * base;
  }

  template <class Base, class Exponent>
  constexpr
  inline
  auto pow(const Base& base, const Exponent& exponent)
  {
    return std::pow(base, exponent);
  }

  /**
  * @brief Computes the square root of a value of type T.
  * @param[in] x Value
  * @tparam T Type of value
  * @returns Square root of value
  */
  template <class T>
  constexpr
  inline
  std::enable_if_t<!std::is_void_v<decltype(std::sqrt(std::declval<T>()))>, T>
  sqrt(const T& x)
  {
    return std::sqrt(x);
  }

  /**
  * @brief Determines if the floating point number is not-a-number (NaN).
  * @param[in] x Value
  * @tparam T Type of value
  * @returns True if value is NaN, false otherwise.
  */
  template <class T>
  constexpr
  inline
  std::enable_if_t<!std::is_void_v<decltype(std::isnan(std::declval<T>()))>, bool>
  isNaN(const T& x)
  {
    return std::isnan(x);
  }

  /**
  * @brief Determines if the floating point number is positive or negative infinity.
  * @param[in] x Value
  * @tparam T Type of value
  * @returns True if value is Inf (or -Inf), false otherwise.
  */
  template <class T>
  constexpr
  inline
  std::enable_if_t<!std::is_void_v<decltype(std::isinf(std::declval<T>()))>, bool>
  isInf(const T& x)
  {
    return std::isinf(x);
  }

  template <class T>
  constexpr
  inline
  std::enable_if_t<!std::is_void_v<decltype(std::cos(std::declval<T>()))>, bool>
  cos(const T& x)
  {
    return std::cos(x);
  }

  template <class T>
  constexpr
  inline
  std::enable_if_t<!std::is_void_v<decltype(std::tan(std::declval<T>()))>, bool>
  tan(const T& x)
  {
   return std::tan(x);
  }

  template <typename T>
  constexpr
  inline
  T sgn(const T& x)
  {
    return (T(0) < x) - (x < T(0));
  }

  template <class T>
  constexpr
  inline
  T binom(const T& n, const T& k)
  {
    assert(T(0) <= n);
    assert(T(0) <= k);
    assert(k <= n);
    T res(1);
    for (T i = 0; i < std::min(k, n - k); ++i)
    {
      res *= (n - i);
      res /= (i + T(1));
    }
    return res;
  }

  template <class T>
  constexpr
  inline
  T factorial(const T& n)
  {
    assert(T(0) <= n);
    T res(1);
    for (T i = T(2); i <= n; ++i)
      res *= i;
    return res;
  }

  template <class T>
  constexpr
  inline
  T permutation(const T& n, const T& k)
  {
    assert(T(0) <= n);
    assert(T(0) <= k);
    assert(k <= n);
    T res(1);
    for (T i = 0; i < k; i++)
        res *= (n - i);
    return res;
  }
}

#endif
