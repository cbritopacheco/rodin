/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_TESTS_CONVERGENCE_REACTIONDIFFUSION_H
#define RODIN_TESTS_CONVERGENCE_REACTIONDIFFUSION_H

#include <cassert>
#include <cmath>
#include <cstdint>

#include "Convergence.h"

namespace Rodin::Tests::Convergence
{
  /** @brief Data for two coupled diffusion fields with reaction matrix
   * @f$R=\begin{pmatrix}1&0.2\\0.2&1\end{pmatrix}@f$ and
   * diffusion coefficients @f$\kappa=(1,2)@f$.
   * Component indices are zero-based. Constant patches use @f$u_i=i+1@f$.
   * Polynomial patches use
   * @f$u_i=(i+1)(1+s)@f$ or @f$u_i=(i+1)(1+s^2)@f$;
   * smooth fields use @f$u_0=e^s, u_1=2e^{-s}@f$, with
   * @f$s=\sum_j x_j@f$. Sources follow from
   * @f$f_i=-\kappa_i\Delta u_i+u_i+0.2u_{1-i}@f$.
   */
  class ReactionDiffusionData
  {
    public:
      enum class Field
      {
        Constant,
        Affine,
        Quadratic,
        Smooth
      };

      ReactionDiffusionData(size_t dimension, Field field)
        : m_dimension(dimension),
          m_field(field)
      {}

      auto getSolution(size_t component) const
      {
        assert(component < 2);
        return Variational::RealFunction(
          [data = *this, component](const Geometry::Point& p) {
            return data.getSolution(p.getPhysicalCoordinates(), component);
          });
      }

      /** @brief Physical-coordinate evaluation for exact-domain comparisons. */
      Real getSolution(const Math::SpatialPoint& x, size_t component) const
      {
        assert(component < 2);
        Real s = 0;
        for (size_t j = 0; j < m_dimension; ++j)
          s += x(j);
        const Real c = Real(component + 1);
        if (m_field == Field::Constant)
          return c;
        if (m_field == Field::Smooth)
          return c * std::exp(component == 0 ? s : -s);
        return c * (1 + (m_field == Field::Affine ? s : s * s));
      }

      auto getGradient(size_t component) const
      {
        return Variational::VectorFunction(
          m_dimension, [data = *this, component](const Geometry::Point& p) {
            return data.getGradient(p.getPhysicalCoordinates(), component);
          });
      }

      Math::SpatialVector<Real> getGradient(
        const Math::SpatialPoint& x, size_t component) const
      {
        assert(component < 2);
        Real s = 0;
        for (size_t j = 0; j < m_dimension; ++j)
          s += x(j);
        Real derivative = Real(component + 1);
        if (m_field == Field::Constant)
          derivative = 0;
        else if (m_field == Field::Quadratic)
          derivative *= 2 * s;
        else if (m_field == Field::Smooth)
          derivative = (component == 0 ? 1 : -1) * getSolution(x, component);
        Math::SpatialVector<Real> value(static_cast<std::uint8_t>(m_dimension));
        for (size_t j = 0; j < m_dimension; ++j)
          value(j) = derivative;
        return value;
      }

      auto getSource(size_t component) const
      {
        const auto exact = getSolution(component);
        const auto other = getSolution(1 - component);
        return Variational::RealFunction([dim = m_dimension, field = m_field, component,
                                           exact, other](const Geometry::Point& p) {
          Real laplacian = 0;
          if (field == Field::Quadratic)
            laplacian = 2 * Real(dim) * Real(component + 1);
          else if (field == Field::Smooth)
            laplacian = Real(dim) * exact(p);
          return -Real(component + 1) * laplacian + exact(p) + 0.2 * other(p);
        });
      }

    private:
      size_t m_dimension;
      Field m_field;
  };
}

#endif
