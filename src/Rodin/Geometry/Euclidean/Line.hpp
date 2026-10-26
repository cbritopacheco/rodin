/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Line.hpp
 * @brief Euclidean geometry declarations and definitions.
 */
#ifndef RODIN_GEOMETRY_EUCLIDEAN_LINE_HPP
#define RODIN_GEOMETRY_EUCLIDEAN_LINE_HPP

#include "Line.h"
#include "Sphere.h"

namespace Rodin::Geometry::Euclidean
{
  template <size_t Dimension, class T>
  typename Line<Dimension, T>::Intersections Line<Dimension, T>::intersect(
    const Sphere<Dimension, T>& sphere) const
  {
    if constexpr (Dimension == 1)
    {
      auto first = sphere.center(), second = sphere.center();
      first[0] += sphere.radius();
      second[0] -= sphere.radius();
      return std::make_pair(first, second);
    }
    else
    {
      const auto nearest = project(sphere.center());
      const T heightSquared =
        sphere.radius() * sphere.radius() - (nearest - sphere.center()).squaredNorm();
      if (heightSquared < T(0))
        return std::nullopt;
      if (heightSquared == T(0))
        return nearest;
      const Vector offset = std::sqrt(heightSquared) * direction();
      return std::make_pair(nearest + offset, nearest - offset);
    }
  }

  template <size_t Dimension, class T>
  Optional<LineSegment<Dimension, T>> Line<Dimension, T>::connect(
    const Sphere<Dimension, T>& sphere) const
  {
    const auto nearest = project(sphere.center());
    const auto radial = nearest - sphere.center();
    const T length = radial.norm();
    if (length <= sphere.radius())
      return {};
    return LineSegment<Dimension, T>(
      nearest, sphere.center() + (sphere.radius() / length) * radial);
  }
}

#endif
