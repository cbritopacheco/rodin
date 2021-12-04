/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_CORE_GEOMETRY_LINE2D_H
#define RODIN_CORE_GEOMETRY_LINE2D_H

#include <variant>

#include "Rodin/Core/Common.h"
#include "Rodin/Core/ForwardDecls.h"
#include "Rodin/Core/Geometry/Geometry.h"

namespace Rodin::Core::Geometry
{
  /**
   * Represents a two-dimensional line in general form.
   * A line \f$ L \f$ is defined as the set of all points \f$ (x, y) \f$ which
   * satisfy the equation:
   * \f[
   * ax + by - c = 0
   * \f]
   */
  template <class T>
  class Line2D : public Geometry<T, Line2D<T>>
  {
    public:
      /**
       * Constructs a line from its coefficients.
       * @param[in] a The \f$ a \f$ coefficient in the equation
       * @param[in] b The \f$ b \f$ coefficient in the equation
       * @param[in] c The \f$ c \f$ coefficient in the equation
       * @warning The behaviour of the line is undefined if constructed with
       * `a == 0 && b == 0`.
       */
      constexpr
      Line2D(T a, T b, T c);

      constexpr
      Line2D(std::initializer_list<T> list);

      /**
       * Constructs a line from two different points.
       * @param[in] p First point
       * @param[in] q Second point
       * @warning The behaviour is undefined if `p == q`.
       */
      constexpr
      Line2D(const Point2D<T>& p, const Point2D<T>& q);

      constexpr
      Line2D(const LineSegment2D<T>& ls);

      /**
       * Evaluates a number \f$ x \f$ using the slope-intercept form of
       * the line:
       * \f[
       * f(x) = mx + k
       * \f]
       * where
       * \f$ m \f$ is the slope, and \f$ k \f$ is the y-intercept.
       * @param[in] x Number to evaluate
       * @returns Result of evaluation
       */
      inline
      constexpr
      std::optional<T> operator()(T x) const;

      /**
       * Evaluates a point \f$ (x, y) \f$ by the following function:
       * \f[
       * f(x, y) = ax + by - c
       * \f]
       *
       * @param[in] p Point to evaluate
       * @returns Result of evaluation
       *
       * @b Example
       *
       * This method can be used to test whether the specified point lies on
       * the upper or lower half spaces defined by the line. For example,
       * @code{.cpp}
       * if (line(p) > 0)
       * {
       *   // Point lies in the upper half space
       * }
       * else
       * {
       *   // Point lies in the lower half space
       * }
       * @endcode
       */
      inline
      constexpr
      T operator()(const Point2D<T>& p) const;

      /**
       * @returns The \f$ a \f$ coefficient of the line.
       */
      inline
      constexpr
      T a() const;

      /**
       * @returns The \f$ b \f$ coefficient of the line.
       */
      inline
      constexpr
      T b() const;

      /**
       * @returns The \f$ c \f$ coefficient of the line.
       */
      inline
      constexpr
      T c() const;

      /**
       * Computes the slope of the line.
       * @retval std::nullopt If the line is vertical.
       * @retval T The slope of the line.
       */
      inline
      constexpr
      std::optional<T> slope() const;

      inline
      constexpr
      std::optional<T> xIntercept() const;

      inline
      constexpr
      std::optional<T> yIntercept() const;

      inline
      constexpr
      bool isVertical() const;

      inline
      constexpr
      bool isHorizontal() const;

      /**
       * Computes the intersection point with another line.
       *
       * If the two lines are given by
       * \f{eqnarray*}{
       *  a_1 x + b_1 y - c_1 &= 0\\
       *  a_2 x + b_2 y - c_2 &= 0
       * \f}
       * then the intersection point (if any) is given by:
       * \f[
       * (x, y) = \left(
       *  \dfrac{c_1 b_2 - b_1 c_2}{a_1 b_2 - a_2 b_1},
       *  \dfrac{a_1 c_2 - a_2 c_1}{a_1 b_2 - a_2 b_1} \right)
       * \f]
       *
       * @param[in] other Line considered
       * @retval std::nullopt If the lines do not intersect
       * @retval Point2D If the lines intersect
       */
      inline
      constexpr
      std::optional<Point2D<T>> intersect(const Line2D<T>& other) const;

      /**
       * Returns the intersection point(s) between the line and the specified
       * circle.
       *
       * @param[in] other Circle considered
       * @retval std::nullopt If no intersection points are found
       * @retval Point2D If one intersection point is found
       * @retval std::pair<Point2D, Point2D> If two intersection points
       * are found
       */
      inline
      constexpr
      std::variant<
        std::nullopt_t,
        Point2D<T>,
        std::pair<Point2D<T>, Point2D<T>>>
      intersect(const Circle<T>& other) const;

      /**
       * Computes the line segment where the start and end points are the
       * points closest to each other in the line and circle, respectively.
       *
       * @param[in] other Line to connect.
       * @retval std::nullopt If the line touches the circle.
       * @retval LineSegment2D If there exists a unique line segment between
       * the two.
       * @note This is equivalent to the result of `point.connect(line).reverse();`
       *
       * @see Point2D<T>::connect()
       */
      inline
      constexpr
      std::optional<LineSegment2D<T>> connect(const Circle<T>& other) const;

      /**
       * Computes the distance between parallel lines.
       *
       * @param[in] other Parallel line
       * @returns The distance between lines.
       * @warning The lines are assumed parallel with respect to each other, no
       * check is made to ensure this is the case.
       */
      inline
      constexpr
      T distance(const Line2D<T>& other) const;

    private:
      T m_a,
        m_b,
        m_c;
  };
}

#include "Line2D.ipp"
#endif
