/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_EUCLIDEAN_SPHERE_H
#define RODIN_GEOMETRY_EUCLIDEAN_SPHERE_H

#include <cmath>
#include "Rodin/Math/Rad.h"
#include "Point.h"
#include "Hyperplane.h"

namespace Rodin::Geometry::Euclidean
{
  /**
   * @brief Sphere @f$ \|x-c\|^2=r^2 @f$ in a compile-time ambient dimension.
   * In one dimension the sphere comprises two endpoints; in two it is a circle.
   * @tparam Dimension Positive ambient dimension.
   * @tparam T Real coordinate scalar type.
   */
  template <size_t Dimension, class T>
  class Sphere : public Base<Sphere<Dimension, T>, T>
  {
    public:
      /// @brief Fixed-size coordinate vector in the ambient dimension.
      using Vector = typename Point<Dimension, T>::Vector;

      /// @brief Constructs a sphere. @param center Center. @param radius Positive radius.
      Sphere(const Point<Dimension, T>& center = Point<Dimension, T>(), T radius = T(1))
        : m_center(center),
          m_radius(radius)
      {
        if constexpr (Dimension == 0)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Sphere requires a positive ambient dimension" << Alert::Raise;
        }
        else
          setRadius(radius);
      }

      /// @brief Gets the ambient dimension. @returns Coordinate count.
      static constexpr size_t getDimension()
      {
        return Dimension;
      }

      /// @brief Evaluates the implicit sphere equation. @param point Point. @returns Squared distance minus squared radius.
      T operator()(const Point<Dimension, T>& point) const
      {
        return (point - m_center).squaredNorm() - m_radius * m_radius;
      }

      /// @brief Gets the radius. @returns Positive radius.
      T radius() const
      {
        return m_radius;
      }

      /// @brief Gets the center. @returns Center point.
      const Point<Dimension, T>& center() const
      {
        return m_center;
      }

      /// @brief Changes the radius. @param radius Positive radius. @returns This sphere.
      Sphere& setRadius(T radius)
      {
        if (!(radius > T(0)))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Sphere requires a positive radius" << Alert::Raise;
        }
        m_radius = radius;
        return *this;
      }

      /// @brief Changes the center. @param center New center. @returns This sphere.
      Sphere& setCenter(const Point<Dimension, T>& center)
      {
        m_center = center;
        return *this;
      }

      /// @brief Computes a radial surface point. @param direction Nonzero spatial direction. @returns Surface point.
      template <class Derived>
      Point<Dimension, T> pointAt(const Eigen::MatrixBase<Derived>& direction) const
      {
        if (direction.size() != Dimension || direction.cols() != 1)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Sphere requires a direction with " << Dimension << " coordinates"
            << Alert::Raise;
        }
        const Vector radial = direction.template cast<T>();
        if (!(radial.array() != T(0)).any())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Sphere requires a nonzero direction" << Alert::Raise;
        }
        if constexpr (Dimension == 1)
        {
          return Point<Dimension, T>(
            m_center[0] + (radial[0] > T(0) ? m_radius : -m_radius));
        }
        else
          return m_center + m_radius * (radial / radial.stableNorm());
      }

      /// @brief Computes distance to the surface. @param point Point. @returns Unsigned surface distance.
      T distance(const Point<Dimension, T>& point) const
      {
        return std::abs((point - m_center).norm() - m_radius);
      }

      /// @brief Computes the tangent hyperplane. @param direction Nonzero radial direction. @returns Tangent at the radial surface point.
      template <class Derived>
      Hyperplane<Dimension, T> tangent(const Eigen::MatrixBase<Derived>& direction) const
      {
        return Hyperplane<Dimension, T>(pointAt(direction), direction);
      }

      /// @brief Evaluates the planar angular parameterization. @param angle Angle. @returns Circle point.
      Point<Dimension, T> operator()(const Math::Rad& angle) const
        requires(Dimension == 2)
      {
        return m_center +
          Point<Dimension, T>(
            m_radius * std::cos(static_cast<T>(static_cast<Real>(angle))),
            m_radius * std::sin(static_cast<T>(static_cast<Real>(angle))));
      }

      /// @brief Computes a planar tangent line. @param angle Angle. @returns Tangent line.
      Line<Dimension, T> tangent(const Math::Rad& angle) const
        requires(Dimension == 2)
      {
        const auto point = operator()(angle);
        const auto radial = point - m_center;
        return Line<Dimension, T>(radial[0], radial[1], radial.dot(point.asVector()));
      }

    private:
      Point<Dimension, T> m_center;
      T m_radius;
  };
}

#include "Line.h"

#endif
