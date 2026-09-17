/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Forward declarations for Euclidean geometry classes.
 */

#ifndef RODIN_GEOMETRY_EUCLIDEAN_FORWARDDECLS_H
#define RODIN_GEOMETRY_EUCLIDEAN_FORWARDDECLS_H

namespace Rodin::Geometry::Euclidean
{
  /**
   * @brief Base class for Euclidean geometric objects using CRTP.
   * @tparam T Scalar type (e.g., float, double)
   * @tparam Derived Derived class type
   * @see <a href="class_rodin_1_1_geometry_1_1_euclidean_1_1_point2_d.html">Point2D</a>
   * @see <a href="class_rodin_1_1_geometry_1_1_euclidean_1_1_circle.html">Circle</a>
   * @see <a href="class_rodin_1_1_geometry_1_1_euclidean_1_1_line2_d.html">Line2D</a>
   * @see <a href="class_rodin_1_1_geometry_1_1_euclidean_1_1_line_segment2_d.html">LineSegment2D</a>
   * @see <a href="class_rodin_1_1_geometry_1_1_euclidean_1_1_rectangle.html">Rectangle</a>
   */
  template <class T, class Derived>
  class Base;

  /**
   * @brief Circle in 2D Euclidean space.
   * @tparam T Scalar type
   * @see <a href="class_rodin_1_1_geometry_1_1_euclidean_1_1_line2_d.html">Line2D</a>
   * @see <a href="class_rodin_1_1_geometry_1_1_euclidean_1_1_point2_d.html">Point2D</a>
   */
  template <class T>
  class Circle;

  /**
   * @brief Infinite line in 2D Euclidean space.
   * @tparam T Scalar type
   * @see <a href="class_rodin_1_1_geometry_1_1_euclidean_1_1_line_segment2_d.html">LineSegment2D</a>
   * @see <a href="class_rodin_1_1_geometry_1_1_euclidean_1_1_point2_d.html">Point2D</a>
   * @see <a href="class_rodin_1_1_geometry_1_1_euclidean_1_1_circle.html">Circle</a>
   */
  template <class T>
  class Line2D;

  /**
   * @brief Point in 2D Euclidean space.
   * @tparam T Scalar type
   * @see <a href="class_rodin_1_1_geometry_1_1_euclidean_1_1_line2_d.html">Line2D</a>
   * @see <a href="class_rodin_1_1_geometry_1_1_euclidean_1_1_circle.html">Circle</a>
   * @see <a href="class_rodin_1_1_geometry_1_1_euclidean_1_1_line_segment2_d.html">LineSegment2D</a>
   */
  template <class T>
  class Point2D;

  /**
   * @brief Line segment in 2D Euclidean space.
   * @tparam T Scalar type
   * @see <a href="class_rodin_1_1_geometry_1_1_euclidean_1_1_line2_d.html">Line2D</a>
   * @see <a href="class_rodin_1_1_geometry_1_1_euclidean_1_1_point2_d.html">Point2D</a>
   */
  template <class T>
  class LineSegment2D;

  /**
   * @brief Axis-aligned rectangle in 2D Euclidean space.
   * @tparam T Scalar type
   * @see <a href="class_rodin_1_1_geometry_1_1_euclidean_1_1_point2_d.html">Point2D</a>
   */
  template <class T>
  class Rectangle;
}

#endif
