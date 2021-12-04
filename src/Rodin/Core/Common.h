/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_CORE_COMMON_H
#define RODIN_CORE_COMMON_H

#include <map>
#include <vector>
#include <cmath>
#include <limits>
#include <optional>
#include <type_traits>

#include <Eigen/Core>

namespace Rodin::Core
{
  /*
   * Removed the definition because we don't need it right now and I
   * don't want Magnum as part of the Rodin::Core module.
   */
  template <class T>
  class Rad;
  // using Rad = Magnum::Math::Rad<T>;

  template <class T>
  inline
  constexpr
  T abs(const T& x);

  template <class T>
  inline
  constexpr
  T sqrt(const T& x);
}

#include "Common.hpp"

#endif
