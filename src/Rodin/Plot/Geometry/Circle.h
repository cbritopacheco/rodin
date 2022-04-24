/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_CORE_GEOMETRY_CIRCLE_H
#define RODIN_CORE_GEOMETRY_CIRCLE_H

#include <variant>

#include <Magnum/Math/Angle.h>

#include "Rodin/Core/Common.h"

#include "ForwardDecls.h"
#include "Base.h"

namespace Rodin::Plot::Geometry
{
  /**
   * \f[
   * (x - a)^2 + (y - b)^2 = r^2
   * \f]
   */
  template <class T>
  class Circle : public Base<T, Circle<T>>
  {
    public:
      constexpr
      Circle(const Point2D<T>& center = Point2D<T>({0, 0}), const T& radius = 1);

      /**
       * This function evaluates the point \f$(x, y) \f$
       * \f[
       * f(x, y) = (x - a)^2 + (y - b)^2 - r^2
       * \f]
       */
      inline
      constexpr
      T operator()(const Point2D<T>& p) const;

      /**
       * This function evaluates the angle \f$\theta \in [0, 2\pi]\f$ using the
       * parametric form of the circle:
       * \f{eqnarray*}{
       * x &= a + r \cos{\theta}\\
       * y &= b + r \sin{\theta}
       * \f}
       * @param[in] angle Angle in the range \f$[0, 2\pi]\f$, which the ray
       *                  from \f$(a, b)\f$ to \f$(x, y)\f$ makes with the
       *                  positive x-axis.
       * @returns The point \f$ (x, y) \f$ on the circle.
       */
      inline
      constexpr
      Point2D<T> operator()(const Magnum::Math::Rad<T>& angle) const;

      /**
       * @returns The radius of the circle.
       */
      inline
      constexpr
      T radius() const;

      /**
       * @returns The center of the circle.
       */
      inline
      constexpr
      Point2D<T> center() const;

      /**
       * Sets the center of the circle.
       * @param[in] center New center of the circle.
       * @returns Reference to self (for method chaining).
       */
      constexpr
      Circle& setCenter(const Point2D<T>& center);

      /**
       * Sets the radius of the circle.
       * @param[in] radius New radius of the circle.
       * @returns Reference to self (for method chaining).
       */
      constexpr
      Circle& setRadius(T radius);

      inline
      constexpr
      Line2D<T> tangent(const Magnum::Math::Rad<T>& angle);

    private:
      Point2D<T>    m_center;
      T             m_radius;
  };
}

#include "Circle.hpp"
#endif
