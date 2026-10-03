/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_TESTS_CONVERGENCE_CONDUCTIVITY_H
#define RODIN_TESTS_CONVERGENCE_CONDUCTIVITY_H

#include <cstdint>
#include "Convergence.h"

namespace Rodin::Tests::Convergence
{
  /** Manufactured data for @f$-\nabla\cdot(\gamma\nabla u)@f$ with
   * @f$\gamma=1+\sum_i x_i@f$, or @f$\gamma=1@f$ for Poisson.
   * Constant, affine and quadratic physical patches supplement the smooth field.
   */
  class ConductivityData
  {
    public:
      enum class Field
      {
        Constant,
        Affine,
        Quadratic,
        Smooth
      };

      ConductivityData(size_t dimension, Field field)
        : m_dimension(dimension),
          m_field(field)
      {}

      auto getCoefficient(bool constant = false) const
      {
        return Variational::RealFunction(
          [dim = m_dimension, constant](const Geometry::Point& p) {
            Real value = 1;
            if (!constant)
              for (size_t d = 0; d < dim; ++d)
                value += p(d);
            return value;
          });
      }

      auto getSolution() const
      {
        return Variational::RealFunction(
          [dim = m_dimension, field = m_field](const Geometry::Point& p) {
            Real value = 1;
            if (field == Field::Constant)
              return value;
            if (field == Field::Smooth)
            {
              value = 1;
              for (size_t d = 0; d < dim; ++d)
                value *= std::sin(Math::Constants::pi() * p(d));
              return 1 + value;
            }
            for (size_t d = 0; d < dim; ++d)
              value += field == Field::Affine ? p(d) : p(d) * p(d);
            return value;
          });
      }

      auto getGradient() const
      {
        return Variational::VectorFunction(
          m_dimension, [dim = m_dimension, field = m_field](const Geometry::Point& p) {
            Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
            const Real pi = Math::Constants::pi();
            for (size_t d = 0; d < dim; ++d)
            {
              value(d) = field == Field::Constant ? 0
                : field == Field::Affine          ? 1
                                                  : 2 * p(d);
              if (field == Field::Smooth)
              {
                value(d) = pi * std::cos(pi * p(d));
                for (size_t j = 0; j < dim; ++j)
                  if (j != d)
                    value(d) *= std::sin(pi * p(j));
              }
            }
            return value;
          });
      }

      auto getSource(bool constantCoefficient = false) const
      {
        const auto gamma = getCoefficient(constantCoefficient);
        const auto gradient = getGradient();
        return Variational::RealFunction(
          [dim = m_dimension, field = m_field, gamma, gradient, constantCoefficient](
            const Geometry::Point& p) {
            Real laplacian = field == Field::Quadratic ? 2 * Real(dim) : 0;
            if (field == Field::Smooth)
            {
              const Real pi = Math::Constants::pi();
              laplacian = -Real(dim) * pi * pi;
              for (size_t d = 0; d < dim; ++d)
                laplacian *= std::sin(pi * p(d));
            }
            const auto derivative = gradient(p);
            Real source = -gamma(p) * laplacian;
            if (!constantCoefficient)
              for (size_t d = 0; d < dim; ++d)
                source -= derivative(d);
            return source;
          });
      }

    private:
      size_t m_dimension;
      Field m_field;
  };
}

#endif
