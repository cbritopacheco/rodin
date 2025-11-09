/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_EUCLIDEAN_SEGMENT2D_H
#define RODIN_GEOMETRY_EUCLIDEAN_SEGMENT2D_H

/**
 * @file
 * @brief Line segment in 2D Euclidean space.
 */

#include <optional>

#include "ForwardDecls.h"
#include "Base.h"

namespace Rodin::Geometry::Euclidean
{
  /**
   * @brief Line segment in 2D Euclidean space.
   *
   * Represents a finite line segment connecting two points @f$ (x_0, y_0) @f$
   * and @f$ (x_1, y_1) @f$ in two-dimensional space.
   *
   * @tparam T Scalar type (e.g., float, double)
   *
   * # Parametric Form
   *
   * The segment can be parameterized as:
   * @f[
   *   (x(t), y(t)) = (x_0, y_0) + t \left((x_1, y_1) - (x_0, y_0)\right)
   * @f]
   * where @f$ t \in [0, 1] @f$ maps to points on the segment, though the
   * parameterization accepts any real @f$ t @f$.
   *
   * @see Point2D, Line2D
   */
  template <class T>
  class LineSegment2D : public Base<LineSegment2D<T>, T>
  {
   public:
    /**
     * @brief Constructs a line segment from two points.
     * @param[in] start Start point @f$ (x_0, y_0) @f$
     * @param[in] end End point @f$ (x_1, y_1) @f$
     */
    constexpr
    LineSegment2D(const Point2D<T>& start, const Point2D<T>& end);

    /**
     * @brief Evaluates the parametric representation at parameter @f$ t @f$.
     *
     * Given @f$ t \in \mathbb{R} @f$, computes the parametric point
     * @f$ (x(t), y(t)) @f$ of the line segment parametrization:
     * @f[
     *   (x(t), y(t)) = (x_0, y_0) + t \left((x_1, y_1) - (x_0, y_0)\right)
     * @f]
     * where @f$(x_0, y_0)@f$ and @f$(x_1, y_1)@f$ are the start and end of the line
     * segment, respectively.
     *
     * @param[in] t Parameter value
     * @returns @f$ (x(t), y(t)) @f$
     *
     * @note For @f$ t \in [0, 1] @f$, the result lies on the segment.
     * Values outside this range extrapolate beyond the segment endpoints.
     */
    inline
    constexpr
    Point2D<T> operator()(const T& t) const;

    /**
     * @brief Computes the length of the segment.
     * @returns Euclidean distance @f$ \|(x_1, y_1) - (x_0, y_0)\| @f$
     */
    inline
    constexpr
    T length() const;

    /**
     * @brief Gets the unit direction vector.
     * @returns Unit vector from `start()` to `end()`
     *
     * Computes:
     * @f[
     *   \frac{(x_1, y_1) - (x_0, y_0)}{\|(x_1, y_1) - (x_0, y_0)\|}
     * @f]
     */
    inline
    constexpr
    Eigen::Vector2<T> direction() const;

    /**
     * @brief Gets the start point.
     * @returns Start point @f$ (x_0, y_0) @f$
     */
    inline
    constexpr
    Point2D<T> start() const;

    /**
     * @brief Gets the end point.
     * @returns End point @f$ (x_1, y_1) @f$
     */
    inline
    constexpr
    Point2D<T> end() const;

    /**
     * @brief Sets the start point.
     * @param[in] start New start point
     * @returns Reference to this segment for method chaining
     */
    inline
    constexpr
    LineSegment2D& setStart(const Point2D<T>& start);

    /**
     * @brief Sets the end point.
     * @param[in] end New end point
     * @returns Reference to this segment for method chaining
     */
    inline
    constexpr
    LineSegment2D& setEnd(const Point2D<T>& end);

    /**
     * @brief Returns the reversed segment.
     * @returns New segment with swapped start and end points
     */
    inline
    constexpr
    LineSegment2D reverse() const;

    /**
     * @brief Computes the intersection point with another segment.
     *
     * Determines if two line segments intersect and returns the intersection
     * point if it exists.
     *
     * @param[in] other Other line segment
     * @returns Optional intersection point
     * @retval std::nullopt If segments do not intersect or are parallel
     * @retval Point2D The intersection point if it exists
     */
    inline
    constexpr
    Optional<Point2D<T>> intersect(const LineSegment2D<T>& other)
    const;

   private:
    Point2D<T>  m_start,  ///< Segment start point
            m_end;    ///< Segment end point
  };
}

#include "LineSegment2D.hpp"

#endif

