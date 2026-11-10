/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_EUCLIDEAN_LINESEGMENT_H
#define RODIN_GEOMETRY_EUCLIDEAN_LINESEGMENT_H

/**
 * @file
 * @brief Line segment with a compile-time Euclidean dimension.
 */

#include <algorithm>
#include "Rodin/Math/Vector.h"

#include "ForwardDecls.h"
#include "Base.h"

namespace Rodin::Geometry::Euclidean
{
  /**
   * @brief Segment @f$ s(t)=p+t(q-p) @f$ between distinct endpoints.
   * The closed segment corresponds to @f$ 0\leq t\leq1 @f$; parameter
   * evaluation also accepts extrapolation. Storage has exactly Dimension components. Intersection returns a point only when unique, including
   * a shared endpoint of collinear segments; positive-length overlaps return
   * no unique point.
   * @tparam Dimension Ambient coordinate count.
   * @tparam T Real coordinate scalar type.
   */
  template <size_t Dimension, class T>
  class LineSegment : public Base<LineSegment<Dimension, T>, T>
  {
    public:
      /**
       * @brief Constructs a line segment from two points.
       * @param[in] start Start point @f$ p @f$
       * @param[in] end End point @f$ q @f$
       */
      constexpr LineSegment(
        const Point<Dimension, T>& start, const Point<Dimension, T>& end);

      /**
       * @brief Evaluates the parametric representation at parameter @f$ t @f$.
       *
       * Given @f$ t \in \mathbb{R} @f$, computes the parametric point
       * @f$ s(t) @f$ of the line segment parametrization.
       * @f[
       *   s(t) = p + t \left(q - p\right)
       * @f]
       * where @f$p@f$ and @f$q@f$ are the start and end of the line
       * segment, respectively.
       *
       * @param[in] t Parameter value
       * @returns @f$ s(t) @f$
       *
       * @note For @f$ t \in [0, 1] @f$, the result lies on the segment.
       * Values outside this range extrapolate beyond the segment endpoints.
       */
      inline constexpr Point<Dimension, T> operator()(const T& t) const;

      /**
       * @brief Computes the length of the segment.
       * @returns Euclidean distance @f$ \|q - p\| @f$
       */
      inline constexpr T length() const;

      /**
       * @brief Gets the unit direction vector.
       * @returns Unit vector from `start()` to `end()`
       *
       * Computes the following quantity.
       * @f[
       *   \frac{q - p}{\|q - p\|}
       * @f]
       */
      inline constexpr typename Point<Dimension, T>::Vector direction() const;

      /**
       * @brief Gets the start point.
       * @returns Start point @f$ p @f$
       */
      inline constexpr Point<Dimension, T> start() const;

      /**
       * @brief Gets the end point.
       * @returns End point @f$ q @f$
       */
      inline constexpr Point<Dimension, T> end() const;

      /**
       * @brief Sets the start point.
       * @param[in] start New start point
       * @returns Reference to this segment for method chaining
       */
      inline constexpr LineSegment& setStart(const Point<Dimension, T>& start);

      /**
       * @brief Sets the end point.
       * @param[in] end New end point
       * @returns Reference to this segment for method chaining
       */
      inline constexpr LineSegment& setEnd(const Point<Dimension, T>& end);

      /**
       * @brief Returns the reversed segment.
       * @returns New segment with swapped start and end points
       */
      inline constexpr LineSegment reverse() const;

      /**
       * @brief Computes the intersection point with another segment.
       *
       * Determines if two line segments intersect and returns the intersection
       * point if it exists.
       *
       * @param[in] other Other line segment
       * @returns Optional intersection point
       * @retval std::nullopt If the intersection has no unique point
       * @retval Point The unique intersection point
       */
      inline constexpr Optional<Point<Dimension, T>> intersect(
        const LineSegment<Dimension, T>& other) const;

      /// @brief Projects onto the closed segment. @param point Point. @returns Closest segment point.
      Point<Dimension, T> project(const Point<Dimension, T>& point) const
      {
        if constexpr (Dimension == 0)
          return point;
        else if constexpr (Dimension == 1)
        {
          return Point<Dimension, T>(std::clamp(
            point[0], std::min(m_start[0], m_end[0]), std::max(m_start[0], m_end[0])));
        }
        else
        {
          const auto displacement = m_end - m_start;
          const T t = displacement.dot(point - m_start) / displacement.squaredNorm();
          return operator()(std::clamp(t, T(0), T(1)));
        }
      }

      /// @brief Computes point distance. @param point Point. @returns Euclidean distance.
      T distance(const Point<Dimension, T>& point) const
      {
        return (point - project(point)).norm();
      }

    private:
      Point<Dimension, T> m_start, ///< Segment start point
        m_end; ///< Segment end point
  };
}

#include "LineSegment.hpp"

#endif
