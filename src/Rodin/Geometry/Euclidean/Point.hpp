/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Point.hpp
 * @brief Point connections in a compile-time Euclidean dimension.
 */
#ifndef RODIN_GEOMETRY_EUCLIDEAN_POINT_HPP
#define RODIN_GEOMETRY_EUCLIDEAN_POINT_HPP

#include "LineSegment.h"

#include "Point.h"

namespace Rodin::Geometry::Euclidean
{
  template <size_t Dimension, class T>
  constexpr Optional<LineSegment<Dimension, T>> Point<Dimension, T>::connect(
    const Point& other) const
  {
    if (*this == other)
      return {};
    else
      return LineSegment<Dimension, T>(*this, other);
  }

  template <size_t Dimension, class T>
  constexpr Optional<LineSegment<Dimension, T>> Point<Dimension, T>::connect(
    const Line<Dimension, T>& line) const
  {
    return connect(line.project(*this));
  }

  template <size_t Dimension, class T>
  constexpr Optional<LineSegment<Dimension, T>> Point<Dimension, T>::connect(
    const Sphere<Dimension, T>& other) const
  {
    const Vector radial = this->asVector() - other.center().asVector();
    const auto length = radial.norm();
    if (length == T{0} || other(*this) == T{0})
      return {};
    return LineSegment<Dimension, T>(*this,
      Point<Dimension, T>(other.center().asVector() + other.radius() * radial / length));
  }
}

#endif
