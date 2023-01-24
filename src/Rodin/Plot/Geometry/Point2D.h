/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_CORE_GEOMETRY_POINT2D_H
#define RODIN_CORE_GEOMETRY_POINT2D_H

#include <optional>

#include <Magnum/Math/Vector2.h>

#include "Rodin/Core/Common.h"

#include "ForwardDecls.h"

namespace Rodin::Plot::Geometry
{
  template <class T>
  class Point2D : public Magnum::Math::Vector2<T>, public Base<T, Point2D<T>>
  {
   public:
    using Magnum::Math::Vector2<T>::Vector2;

    /**
     * Computes the line segment where the start point is the current point
     * and the end point is the other point.
     *
     * @param[in] other Point to connect.
     * @returns `std::nullopt` if the point is equal to the other point, else
     * the line segment connecting both.
     */
    inline
    constexpr
    std::optional<LineSegment2D<T>> connect(const Point2D<T>& other) const;

    /**
     * Computes the line segment where the start point is the current point
     * and the end point is the nearest point on the line.
     *
     * @param[in] other Line to connect.
     * @returns `std::nullopt` if the point touches the line, otherwise a
     * line segment.
     */
    inline
    constexpr
    std::optional<LineSegment2D<T>> connect(const Line2D<T>& other) const;

    /**
     * Computes the line segment where the start point is the current point and
     * the end point is the nearest point on the circle.
     *
     * @param[in] other Circle to connect.
     * @returns `std::nullopt` if the point touches the circle, otherwise a
     * line segment.
     */
    inline
    constexpr
    std::optional<LineSegment2D<T>> connect(const Circle<T>& other) const;
  };
}

#include "Point2D.hpp"

#endif
