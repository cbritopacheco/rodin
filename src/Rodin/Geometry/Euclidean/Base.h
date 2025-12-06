/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_EUCLIDEAN_H
#define RODIN_GEOMETRY_EUCLIDEAN_H

/**
 * @file
 * @brief Base class for Euclidean geometric objects using CRTP pattern.
 */

#include <utility>
#include <type_traits>

#include "Rodin/Math/Common.h"

#include "Concepts.h"

namespace Rodin::Geometry::Euclidean
{
  /**
   * @brief Base class for Euclidean geometric objects using CRTP.
   *
   * This class serves as the base for all Euclidean geometry classes,
   * providing a common interface through the Curiously Recurring Template
   * Pattern (CRTP).
   *
   * @tparam Derived The derived class type
   * @tparam T Scalar type (e.g., float, double)
   *
   * @note This is an abstract base class with no functionality, serving
   * only to establish the type hierarchy for Euclidean geometric objects.
   *
   * @see Point2D, Circle, Line2D, LineSegment2D, Rectangle
   */
  template <class Derived, class T>
  class Base
  {
    public:
      /**
       * @brief Virtual destructor for polymorphic behavior.
       */
      virtual ~Base() = default;
  };
}

#endif
