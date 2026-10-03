/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_TESTS_CONVERGENCE_HELMHOLTZ_H
#define RODIN_TESTS_CONVERGENCE_HELMHOLTZ_H

#include <cmath>
#include <cstdint>
#include "Convergence.h"

namespace Rodin::Tests::Convergence
{
  /** @brief Manufactured complex fields for
   * @f$-\Delta u-\tfrac14u=f@f$ in physical coordinates.
   * With @f$s=\sum_jx_j@f$, patches are
   * @f$u=1+2i@f$, @f$u=1+2i+(1+i/2)s@f$ and
   * @f$u=1+2i+(1+i/2)s^2@f$;
   * the smooth field is @f$u=e^{is}@f$.
   */
  class HelmholtzData
  {
    public:
      enum class Field
      {
        Constant,
        Affine,
        Quadratic,
        Smooth
      };

      HelmholtzData(size_t dimension, Field field)
        : m_dimension(dimension),
          m_field(field)
      {}

      auto getSolution() const
      {
        return Variational::ComplexFunction([data = *this](const Geometry::Point& p) {
          return data.getSolution(p.getPhysicalCoordinates());
        });
      }

      Complex getSolution(const Math::SpatialPoint& x) const
      {
        Real s = 0;
        for (size_t j = 0; j < m_dimension; ++j)
          s += x(j);
        if (m_field == Field::Smooth)
          return std::exp(Complex(0, s));
        if (m_field == Field::Constant)
          return Complex(1, 2);
        return Complex(1, 2) + Complex(1, 0.5) * (m_field == Field::Affine ? s : s * s);
      }

      auto getGradient() const
      {
        return [data = *this](const Geometry::Point& p) {
          return data.getGradient(p.getPhysicalCoordinates());
        };
      }

      Math::SpatialVector<Complex> getGradient(const Math::SpatialPoint& x) const
      {
        Real s = 0;
        for (size_t j = 0; j < m_dimension; ++j)
          s += x(j);
        Complex derivative(1, 0.5);
        if (m_field == Field::Constant)
          derivative = 0;
        else if (m_field == Field::Quadratic)
          derivative *= 2 * s;
        else if (m_field == Field::Smooth)
          derivative = Complex(0, 1) * getSolution(x);
        Math::SpatialVector<Complex> value(static_cast<std::uint8_t>(m_dimension));
        for (size_t j = 0; j < m_dimension; ++j)
          value(j) = derivative;
        return value;
      }

      auto getSource() const
      {
        const auto exact = getSolution();
        return Variational::ComplexFunction(
          [dim = m_dimension, field = m_field, exact](const Geometry::Point& p) {
            Complex laplacian(0);
            if (field == Field::Quadratic)
              laplacian = 2 * Real(dim) * Complex(1, 0.5);
            else if (field == Field::Smooth)
              laplacian = -Real(dim) * exact(p);
            return -laplacian - 0.25 * exact(p);
          });
      }

    private:
      size_t m_dimension;
      Field m_field;
  };
}

#endif
