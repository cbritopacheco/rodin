/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_EUCLIDEAN_POINT2D_H
#define RODIN_GEOMETRY_EUCLIDEAN_POINT2D_H

/**
 * @file
 * @brief Point in 2D Euclidean space.
 */

#include <optional>

#include <Eigen/Core>

#include "ForwardDecls.h"

namespace Rodin::Geometry::Euclidean
{
  /**
   * @brief Point in 2D Euclidean space.
   *
   * Represents a point @f$ p = (x, y) @f$ in two-dimensional Euclidean space.
   * This class inherits from Eigen::Vector2 to provide vector arithmetic
   * operations.
   *
   * @tparam T Scalar type (e.g., float, double)
   *
   * # Operations
   *
   * The class supports:
   * - Vector arithmetic (addition, subtraction, scaling) via Eigen::Vector2
   * - Connection to other geometric objects (Line2D, Circle, other points)
   * - Distance computations
   *
   * @note This class uses multiple inheritance from Eigen::Vector2 and Base,
   * providing both vector operations and CRTP-based geometry interface.
   *
   * @see Line2D, Circle, LineSegment2D
   */
  template <class T>
  class Point2D : public Eigen::Vector2<T>, public Base<Point2D<T>, T>
  {
    public:
      using Eigen::Vector2<T>::Vector2;
      using Eigen::Vector2<T>::operator=;

      /**
       * @brief Gets mutable reference to underlying Eigen vector.
       * @returns Reference to Eigen::Vector2
       */
      inline
      Eigen::Vector2<T>& asVector()
      {
        return static_cast<Eigen::Vector2<T>&>(*this);
      }

      /**
       * @brief Gets const reference to underlying Eigen vector.
       * @returns Const reference to Eigen::Vector2
       */
      inline
      const Eigen::Vector2<T>& asVector() const
      {
        return static_cast<const Eigen::Vector2<T>&>(*this);
      }

      /**
       * @brief Computes the line segment connecting this point to another.
       *
       * Creates a line segment where the start point is the current point
       * and the end point is the other point.
       *
       * @param[in] other Point to connect to
       * @returns Optional line segment connecting the points
       * @retval std::nullopt If the points are equal (no segment can be formed)
       * @retval LineSegment2D If the points are distinct
       *
       * @see LineSegment2D
       */
      inline
      constexpr
      Optional<LineSegment2D<T>> connect(const Point2D& other) const;

      /**
       * @brief Computes the line segment connecting this point to nearest point on a line.
       *
       * Creates a line segment where the start point is the current point
       * and the end point is the nearest point on the given line.
       *
       * @param[in] other Line to connect to
       * @returns Optional line segment
       * @retval std::nullopt If the point lies on the line
       * @retval LineSegment2D Perpendicular segment from point to line
       *
       * @see Line2D, LineSegment2D
       */
      inline
      constexpr
      Optional<LineSegment2D<T>> connect(const Line2D<T>& other) const;

      /**
       * @brief Computes the line segment connecting this point to nearest point on a circle.
       *
       * Creates a line segment where the start point is the current point and
       * the end point is the nearest point on the circle's circumference.
       *
       * @param[in] other Circle to connect to
       * @returns Optional line segment
       * @retval std::nullopt If the point lies on the circle
       * @retval LineSegment2D Segment from point to nearest point on circle
       *
       * @see Circle, LineSegment2D
       */
      inline
      constexpr
      Optional<LineSegment2D<T>> connect(const Circle<T>& other) const;
  };
}

#include "Point2D.hpp"

#endif
