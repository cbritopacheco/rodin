/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_EUCLIDEAN_H
#define RODIN_GEOMETRY_EUCLIDEAN_H

#include <utility>
#include <type_traits>

#include "Rodin/Math/Common.h"

#include "Concepts.h"

namespace Rodin::Geometry::Euclidean
{
  template <class Derived, class T>
  class Base
  {
    public:
      virtual ~Base() = default;
  };
}

#endif
