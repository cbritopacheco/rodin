/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_EUCLIDEAN_POINT2D_H
#define RODIN_GEOMETRY_EUCLIDEAN_POINT2D_H

#include <optional>

#include <Eigen/Core>

#include "ForwardDecls.h"

namespace Rodin::Geometry::Euclidean
{
  template <class T>
  class Point2D : public Eigen::Vector2<T>, public Base<Point2D<T>, T>
  {
    public:
      using Eigen::Vector2<T>::Vector2;
      using Eigen::Vector2<T>::operator=;

      inline
      Eigen::Vector2<T>& asVector()
      {
        return static_cast<Eigen::Vector2<T>&>(*this);
      }

      inline
      const Eigen::Vector2<T>& asVector() const
      {
        return static_cast<const Eigen::Vector2<T>&>(*this);
      }

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
      std::optional<LineSegment2D<T>> connect(const Point2D& other) const;

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
