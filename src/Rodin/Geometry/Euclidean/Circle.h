/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_EUCLIDEAN_CIRCLE_H
#define RODIN_GEOMETRY_EUCLIDEAN_CIRCLE_H

/**
 * @file
 * @brief Circle in 2D Euclidean space.
 */

#include <variant>

#include "Rodin/Math/Rad.h"
#include "ForwardDecls.h"
#include "Base.h"

namespace Rodin::Geometry::Euclidean
{
  /**
   * @brief Circle in 2D Euclidean space.
   *
   * Represents a circle defined by the equation:
   * @f[
   *   (x - a)^2 + (y - b)^2 = r^2
   * @f]
   * where @f$ (a, b) @f$ is the center and @f$ r @f$ is the radius.
   *
   * @tparam T Scalar type (e.g., float, double)
   *
   * # Parametric Form
   *
   * The circle can also be parameterized using angle @f$ \theta \in [0, 2\pi] @f$:
   * @f{eqnarray*}{
   *   x(\theta) &= a + r \cos\theta \\
   *   y(\theta) &= b + r \sin\theta
   * @f}
   *
   * # Operations
   *
   * The class supports:
   * - Point evaluation using implicit and parametric forms
   * - Tangent line computation
   * - Intersection with other geometric objects
   * - Distance computations
   *
   * @see Point2D, Line2D, LineSegment2D
   */
  template <class T>
  class Circle : public Base<Circle<T>, T>
  {
   public:
    /**
     * @brief Constructs a circle with given center and radius.
     * @param[in] center Center point @f$ (a, b) @f$
     * @param[in] radius Radius @f$ r @f$
     */
    constexpr
    Circle(const Point2D<T>& center = Point2D<T>({0, 0}), const T& radius = 1);

    /**
     * @brief Evaluates the implicit circle function at a point.
     *
     * Computes:
     * @f[
     *   f(x, y) = (x - a)^2 + (y - b)^2 - r^2
     * @f]
     *
     * @param[in] p Point @f$ (x, y) @f$ to evaluate
     * @returns Function value
     * @retval 0 If point is on the circle
     * @retval >0 If point is outside the circle
     * @retval <0 If point is inside the circle
     */
    inline
    constexpr
    T operator()(const Point2D<T>& p) const;

    /**
     * @brief Evaluates the parametric form of the circle at an angle.
     *
     * Computes the point on the circle using parametric equations:
     * @f{eqnarray*}{
     *   x &= a + r \cos\theta \\
     *   y &= b + r \sin\theta
     * @f}
     *
     * @param[in] angle Angle @f$ \theta \in [0, 2\pi] @f$
     * @returns Point @f$ (x, y) @f$ on the circle
     */
    inline
    constexpr
    Point2D<T> operator()(const Math::Rad& angle) const;

    /**
     * @brief Gets the circle's radius.
     * @returns Radius @f$ r @f$
     */
    inline
    constexpr
    T radius() const;

    /**
     * @brief Gets the circle's center.
     * @returns Center point @f$ (a, b) @f$
     */
    inline
    constexpr
    Point2D<T> center() const;

    /**
     * @brief Sets the circle's center.
     * @param[in] center New center point
     * @returns Reference to this circle for method chaining
     */
    constexpr
    Circle& setCenter(const Point2D<T>& center);

    /**
     * @brief Sets the circle's radius.
     * @param[in] radius New radius
     * @returns Reference to this circle for method chaining
     */
    constexpr
    Circle& setRadius(T radius);

    /**
     * @brief Computes the tangent line at a given angle.
     *
     * Returns the line tangent to the circle at the point determined
     * by the parametric angle @f$ \theta @f$.
     *
     * @param[in] angle Angle @f$ \theta @f$ determining the tangent point
     * @returns Tangent line at @f$ (a + r\cos\theta, b + r\sin\theta) @f$
     */
    inline
    constexpr
    Line2D<T> tangent(const Math::Rad& angle);

    private:
      Point2D<T>    m_center;  ///< Circle center
      T             m_radius;  ///< Circle radius
  };
}

#include "Circle.hpp"
#endif
