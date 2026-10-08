/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_TESTS_CONVERGENCE_STOKES_H
#define RODIN_TESTS_CONVERGENCE_STOKES_H

#include "Convergence.h"

namespace Rodin::Tests::Convergence
{
  struct StokesErrors
  {
      ErrorNorms velocity;
      ErrorNorms pressure;
      Real divergence;
  };

  /** @brief Unit-box data for @f$-\Delta u+\nabla p=f@f$,
   * @f$\nabla\cdot u=0@f$. The default pressure has zero unit-box mean.
   * The velocity has only its first component nonzero, depending on the
   * second coordinate; polynomial patches and smooth rate data are
   * therefore divergence-free in both two and three dimensions.
   * For degree-refinement data, @f$u=\sin(\pi x_1)e_0@f$ and
   * @f$p=\cos(\pi x_0)@f$. Both are nonpolynomial, and the forcing is
   * @f$f=(\pi^2\sin(\pi x_1)-\pi\sin(\pi x_0))e_0@f$.
   * An optional nonzero shear axis replaces the default second coordinate.
   * It remains distinct from the velocity component, preserving zero divergence.
   * An optional constant pressure offset leaves the velocity, pressure gradient,
   * and body force unchanged. It changes the mean and the traction by
   * @f$-c n@f$, and is intended for natural-boundary verification rather than
   * the zero-mean fully prescribed-velocity workload.
   * Physical-coordinate overloads support exact-domain lift measurements.
   */
  class StokesData
  {
    public:
      enum class Field
      {
        Affine,
        Quadratic,
        Cubic,
        Quartic,
        Smooth
      };

      StokesData(
        size_t dimension, Field field, size_t shearAxis = 1, Real pressureOffset = 0)
        : m_dimension(dimension),
          m_field(field),
          m_shearAxis(shearAxis),
          m_pressureOffset(pressureOffset)
      {
        assert(dimension == 2 || dimension == 3);
        assert(shearAxis > 0 && shearAxis < dimension);
      }

      auto getVelocity() const
      {
        return Variational::VectorFunction(
          m_dimension, [data = *this](const Geometry::Point& x) {
            return data.getVelocity(x.getPhysicalCoordinates());
          });
      }

      Math::SpatialVector<Real> getVelocity(const Math::SpatialPoint& x) const
      {
        Math::SpatialVector<Real> u(static_cast<std::uint8_t>(m_dimension));
        u.setZero();
        const Real t = x(m_shearAxis);
        u(0) = m_field == Field::Smooth ? std::sin(Math::Constants::pi() * t)
          : m_field == Field::Quartic   ? t * t * t * t
          : m_field == Field::Affine    ? t
          : m_field == Field::Quadratic ? t * t
                                        : t * t * t;
        return u;
      }

      auto getPressure() const
      {
        return Variational::RealFunction([data = *this](const Geometry::Point& x) {
          return data.getPressure(x.getPhysicalCoordinates());
        });
      }

      Real getPressure(const Math::SpatialPoint& x) const
      {
        const Real pressure = m_field == Field::Smooth
          ? std::cos(Math::Constants::pi() * x(0))
          : m_field == Field::Quartic ? x(0) * x(0) * x(0) - Real(1) / 4
          : m_field == Field::Cubic   ? x(0) * x(0) - Real(1) / 3
                                      : x(0) - 0.5;
        return pressure + m_pressureOffset;
      }

      auto getForcing() const
      {
        return Variational::VectorFunction(m_dimension,
          [dim = m_dimension, field = m_field, axis = m_shearAxis](
            const Geometry::Point& x) {
            Math::SpatialVector<Real> f(static_cast<std::uint8_t>(dim));
            f.setZero();
            const Real pi = Math::Constants::pi();
            f(0) = field == Field::Smooth
              ? pi * pi * std::sin(pi * x(axis)) - pi * std::sin(pi * x(0))
              : field == Field::Quartic   ? -12 * x(axis) * x(axis) + 3 * x(0) * x(0)
              : field == Field::Affine    ? 1
              : field == Field::Quadratic ? -1
                                          : -6 * x(axis) + 2 * x(0);
            return f;
          });
      }

      auto getVelocityJacobian() const
      {
        return [data = *this](const Geometry::Point& x) {
          return data.getVelocityJacobian(x.getPhysicalCoordinates());
        };
      }

      Math::SpatialMatrix<Real> getVelocityJacobian(const Math::SpatialPoint& x) const
      {
        Math::SpatialMatrix<Real> j(
          static_cast<std::uint8_t>(m_dimension), static_cast<std::uint8_t>(m_dimension));
        j.setZero();
        const Real pi = Math::Constants::pi();
        const Real t = x(m_shearAxis);
        j(0, m_shearAxis) = m_field == Field::Smooth ? pi * std::cos(pi * t)
          : m_field == Field::Quartic                ? 4 * t * t * t
          : m_field == Field::Affine                 ? 1
          : m_field == Field::Quadratic              ? 2 * t
                                                     : 3 * t * t;
        return j;
      }

      /** @brief Physical viscous stress for manufactured traction data.
       * @f$\sigma(u,p)=\nu(\nabla u+\nabla u^T)-pI@f$.
       * This is pointwise data evaluation and has no collective semantics.
       */
      Math::SpatialMatrix<Real> getStress(
        const Math::SpatialPoint& x, Real viscosity = 1) const
      {
        const auto gradient = getVelocityJacobian(x);
        const Real pressure = getPressure(x);
        Math::SpatialMatrix<Real> stress(
          static_cast<std::uint8_t>(m_dimension), static_cast<std::uint8_t>(m_dimension));
        for (size_t i = 0; i < m_dimension; ++i)
        {
          for (size_t j = 0; j < m_dimension; ++j)
          {
            stress(i, j) = viscosity * (gradient(i, j) + gradient(j, i)) -
              (i == j ? pressure : Real(0));
          }
        }
        return stress;
      }

      auto getPressureGradient() const
      {
        return Variational::VectorFunction(
          m_dimension, [data = *this](const Geometry::Point& x) {
            return data.getPressureGradient(x.getPhysicalCoordinates());
          });
      }

      Math::SpatialVector<Real> getPressureGradient(const Math::SpatialPoint& x) const
      {
        Math::SpatialVector<Real> g(static_cast<std::uint8_t>(m_dimension));
        g.setZero();
        g(0) = m_field == Field::Smooth
          ? -Math::Constants::pi() * std::sin(Math::Constants::pi() * x(0))
          : m_field == Field::Quartic ? 3 * x(0) * x(0)
          : m_field == Field::Cubic   ? 2 * x(0)
                                      : 1;
        return g;
      }

    private:
      size_t m_dimension;
      Field m_field;
      size_t m_shearAxis;
      Real m_pressureOffset;
  };
}

#endif
