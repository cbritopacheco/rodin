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
   * @see Point2D, Circle, Line2D, LineSegment2D, Rectangle
   */
  template <class T, class Derived>
  class Base;

  /**
   * @brief Circle in 2D Euclidean space.
   * @tparam T Scalar type
   * @see Line2D, Point2D
   */
  template <class T>
  class Circle;

  /**
   * @brief Infinite line in 2D Euclidean space.
   * @tparam T Scalar type
   * @see LineSegment2D, Point2D, Circle
   */
  template <class T>
  class Line2D;

  /**
   * @brief Point in 2D Euclidean space.
   * @tparam T Scalar type
   * @see Line2D, Circle, LineSegment2D
   */
  template <class T>
  class Point2D;

  /**
   * @brief Line segment in 2D Euclidean space.
   * @tparam T Scalar type
   * @see Line2D, Point2D
   */
  template <class T>
  class LineSegment2D;

  /**
   * @brief Axis-aligned rectangle in 2D Euclidean space.
   * @tparam T Scalar type
   * @see Point2D
   */
  template <class T>
  class Rectangle;
}

#endif

