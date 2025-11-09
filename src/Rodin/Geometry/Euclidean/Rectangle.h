/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_EUCLIDEAN_RECTANGLE_H
#define RODIN_GEOMETRY_EUCLIDEAN_RECTANGLE_H

/**
 * @file
 * @brief Axis-aligned rectangle in 2D Euclidean space.
 */

#include "ForwardDecls.h"

namespace Rodin::Geometry::Euclidean
{
  /**
   * @brief Axis-aligned rectangle in 2D Euclidean space.
   *
   * Represents an axis-aligned rectangle defined by its bottom-left and
   * top-right corners. The rectangle's sides are parallel to the coordinate
   * axes.
   *
   * @tparam T Scalar type (e.g., float, double)
   *
   * # Definition
   *
   * A rectangle is defined by two points:
   * - Bottom-left corner: @f$ (x_{\min}, y_{\min}) @f$
   * - Top-right corner: @f$ (x_{\max}, y_{\max}) @f$
   *
   * All points @f$ (x, y) @f$ satisfying:
   * @f[
   *   x_{\min} \leq x \leq x_{\max} \quad \text{and} \quad y_{\min} \leq y \leq y_{\max}
   * @f]
   * are contained in the rectangle.
   *
   * @see Point2D
   */
  template <class T>
  class Rectangle
  {
    public:
      /**
       * @brief Corner enumeration for rectangle corners.
       */
      enum Corner
      {
        BottomLeft,   ///< Bottom-left corner @f$ (x_{\min}, y_{\min}) @f$
        BottomRight,  ///< Bottom-right corner @f$ (x_{\max}, y_{\min}) @f$
        TopLeft,      ///< Top-left corner @f$ (x_{\min}, y_{\max}) @f$
        TopRight      ///< Top-right corner @f$ (x_{\max}, y_{\max}) @f$
      };

      /**
       * @brief Constructs a rectangle from two opposite corners.
       * @param[in] bottomLeft Bottom-left corner @f$ (x_{\min}, y_{\min}) @f$
       * @param[in] topRight Top-right corner @f$ (x_{\max}, y_{\max}) @f$
       */
      constexpr
      Rectangle(const Point2D<T>& bottomLeft, const Point2D<T>& topRight);

      /**
       * @brief Tests if a point is contained in the rectangle.
       * @param[in] p Point to test
       * @returns True if point is inside or on the boundary
       *
       * A point @f$ (x, y) @f$ is contained if:
       * @f[
       *   x_{\min} \leq x \leq x_{\max} \quad \text{and} \quad y_{\min} \leq y \leq y_{\max}
       * @f]
       */
      inline
      constexpr
      bool contains(const Point2D<T>& p) const;

      /**
       * @brief Tests if a point is contained in the rectangle.
       * @param[in] x X-coordinate
       * @param[in] y Y-coordinate
       * @returns True if point @f$ (x, y) @f$ is inside or on the boundary
       */
      inline
      constexpr
      bool contains(const T& x, const T& y) const
      {
        return contains(Point2D<T>{x, y});
      }

      /**
       * @brief Gets a specific corner of the rectangle.
       * @tparam c Corner identifier (BottomLeft, BottomRight, TopLeft, TopRight)
       * @returns Corner point
       */
      template <Corner c>
      inline
      constexpr
      Point2D<T> getCorner() const;

      /**
       * @brief Computes the height of the rectangle.
       * @returns Height @f$ y_{\max} - y_{\min} @f$
       */
      inline
      constexpr
      T height() const;

      /**
       * @brief Computes the width of the rectangle.
       * @returns Width @f$ x_{\max} - x_{\min} @f$
       */
      inline
      constexpr
      T width() const;

      /**
       * @brief Computes the area of the rectangle.
       * @returns Area @f$ (x_{\max} - x_{\min})(y_{\max} - y_{\min}) @f$
       */
      inline
      constexpr
      T area() const;

    private:
      Point2D<T>  m_bottomLeft,  ///< Bottom-left corner
                  m_topRight;    ///< Top-right corner
  };
}

#include "Rectangle.hpp"

#endif
