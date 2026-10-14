/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ForwardDecls.h
 * @brief Euclidean geometry declarations and definitions.
 */
#ifndef RODIN_GEOMETRY_EUCLIDEAN_FORWARDDECLS_H
#define RODIN_GEOMETRY_EUCLIDEAN_FORWARDDECLS_H

#include <cstddef>
#include "Rodin/Types.h"

namespace Rodin::Geometry::Euclidean
{
  template <class Derived, class T>
  class Base;
  template <size_t Dimension, class T = Real>
  class Point;
  template <size_t Dimension, class T = Real>
  class LineSegment;
  template <size_t Dimension, class T = Real>
  class Line;
  template <size_t Dimension, class T = Real>
  class Sphere;
  template <size_t Dimension, class T = Real>
  class Hyperplane;
  template <size_t Dimension, class T = Real>
  class Box;
  template <class T>
  using Point2D = Point<2, T>;
  template <class T>
  using LineSegment2D = LineSegment<2, T>;
  template <class T>
  using Line2D = Line<2, T>;
  template <class T>
  using Circle = Sphere<2, T>;
  template <class T>
  using Rectangle = Box<2, T>;
}

#endif
