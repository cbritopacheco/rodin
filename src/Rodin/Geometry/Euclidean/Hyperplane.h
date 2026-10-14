/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_EUCLIDEAN_HYPERPLANE_H
#define RODIN_GEOMETRY_EUCLIDEAN_HYPERPLANE_H

#include "Point.h"

namespace Rodin::Geometry::Euclidean
{
  /**
   * @brief Affine hyperplane @f$ n\cdot(x-p)=0 @f$.
   * @tparam Dimension Positive ambient dimension.
   * @tparam T Real coordinate scalar type.
   */
  template <size_t Dimension, class T>
  class Hyperplane : public Base<Hyperplane<Dimension, T>, T>
  {
    public:
      using Vector = typename Point<Dimension, T>::Vector;
      /// @brief Constructs a hyperplane. @param origin Anchor point. @param normal Nonzero normal with Dimension components.
      template <class Derived>
      Hyperplane(
        const Point<Dimension, T>& origin, const Eigen::MatrixBase<Derived>& normal)
        : m_origin(origin),
          m_normal(Vector::Zero())
      {
        if (normal.size() != Dimension || normal.cols() != 1)
          Alert::MemberFunctionException(*this, __func__)
            << "Hyperplane requires a normal with " << Dimension << " coordinates"
            << Alert::Raise;
        m_normal = normal.template cast<T>();
        if (!(m_normal.array() != T(0)).any())
          Alert::MemberFunctionException(*this, __func__)
            << "Hyperplane requires a nonzero normal" << Alert::Raise;
      }
      /// @brief Gets the ambient dimension. @returns Coordinate count.
      static constexpr size_t getDimension()
      {
        return Dimension;
      }
      /// @brief Gets the anchor. @returns Point on the hyperplane.
      const Point<Dimension, T>& origin() const
      {
        return m_origin;
      }
      /// @brief Gets the normal. @returns Normal vector.
      const Vector& normal() const
      {
        return m_normal;
      }
      /// @brief Evaluates the implicit equation. @param point Point. @returns Signed equation value.
      T operator()(const Point<Dimension, T>& point) const
      {
        return m_normal.dot(point - m_origin);
      }
      /// @brief Projects orthogonally. @param point Point. @returns Closest hyperplane point.
      Point<Dimension, T> project(const Point<Dimension, T>& point) const
      {
        if constexpr (Dimension == 0)
          return point;
        else if constexpr (Dimension == 1)
          return m_origin;
        else
        {
          const Vector unit = m_normal / m_normal.stableNorm();
          return point - unit.dot(point - m_origin) * unit;
        }
      }
      /// @brief Computes distance. @param point Point. @returns Unsigned distance.
      T distance(const Point<Dimension, T>& point) const
      {
        return (point - project(point)).norm();
      }

    private:
      Point<Dimension, T> m_origin;
      Vector m_normal;
  };
}

#endif
