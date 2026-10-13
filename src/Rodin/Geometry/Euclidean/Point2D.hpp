/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Point2D.hpp
 * @brief Point primitive for 2D Euclidean geometry computations.
 */
#ifndef RODIN_GEOMETRY_EUCLIDEAN_POINT2D_HPP
#define RODIN_GEOMETRY_EUCLIDEAN_POINT2D_HPP

#include "Line2D.h"
#include "Circle.h"
#include "LineSegment2D.h"

#include "Point2D.h"

namespace Rodin::Geometry::Euclidean
{
  template <class T>
  constexpr Optional<LineSegment2D<T>> Point2D<T>::connect(const Point2D& other) const
  {
    if (*this == other)
      return {};
    else
      return LineSegment2D<T>(*this, other);
  }

  template <class T>
  constexpr
  Optional<LineSegment2D<T>> Point2D<T>::connect(const Line2D<T>& line) const
  {
    assert(line.a() != T{0} || line.b() != T{0});

    if (line(*this) > T{0} || line(*this) < T{0}) // Point is not touching the line
    {
      return LineSegment2D<T>(*this,
        Point2D<T>(
          {line.b() * (line.b() * this->x() - line.a() * this->y()) + line.a() * line.c(),
            line.a() * (-line.b() * this->x() + line.a() * this->y()) +
              line.b() * line.c()}) /
          (line.a() * line.a() + line.b() * line.b()));
    }
    else
    {
      return {};
    }
  }

  template <class T>
  constexpr
  Optional<LineSegment2D<T>> Point2D<T>::connect(const Circle<T>& other) const
  {
    const Math::SpatialVector<T> radial = this->asVector() - other.center().asVector();
    const auto length = radial.norm();
    if (length == T{0} || other(*this) == T{0})
      return {};
    return LineSegment2D<T>(
      *this, Point2D<T>(other.center().asVector() + other.radius() * radial / length));
  }
}

#endif
