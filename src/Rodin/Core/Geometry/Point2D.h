/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_CORE_GEOMETRY_POINT2D_H
#define RODIN_CORE_GEOMETRY_POINT2D_H

#include <optional>

#include "Rodin/Core/Common.h"
#include "Rodin/Core/ForwardDecls.h"
#include "Rodin/Core/Geometry/Geometry.h"

namespace Rodin::Core::Geometry
{
  template <class T>
  class Point2D : public Eigen::Vector2<T>, public Geometry<T, Point2D<T>>
  {
    public:
      typedef Eigen::Vector2<T> Base;

      constexpr
      Point2D(void);

      constexpr
      Point2D(T x, T y);

      constexpr
      Point2D(std::initializer_list<T> list);

      constexpr
      Point2D(const Eigen::Vector2<T>& other);

      template<typename OtherDerived>
      constexpr
      Point2D(const Eigen::MatrixBase<OtherDerived>& other);

      template<typename OtherDerived>
      constexpr
      Point2D& operator=(const Eigen::MatrixBase<OtherDerived>& other);

      constexpr
      Point2D& operator=(std::initializer_list<T> list);

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

#include "Point2D.ipp"

#endif
