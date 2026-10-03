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
   * @f$\nabla\cdot u=0@f$, and @f$\int_\Omega p=0@f$.
   * The velocity has only its first component nonzero, depending on the
   * second coordinate; polynomial patches and smooth rate data are
   * therefore divergence-free in both two and three dimensions.
   * For degree-refinement data, @f$u=\sin(\pi x_1)e_0@f$ and
   * @f$p=\cos(\pi x_0)@f$. Both are nonpolynomial, and the forcing is
   * @f$f=(\pi^2\sin(\pi x_1)-\pi\sin(\pi x_0))e_0@f$.
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

      StokesData(size_t dimension, Field field)
        : m_dimension(dimension),
          m_field(field)
      {
        assert(dimension == 2 || dimension == 3);
      }

      auto getVelocity() const
      {
        return Variational::VectorFunction(
          m_dimension, [dim = m_dimension, field = m_field](const Geometry::Point& x) {
            Math::SpatialVector<Real> u(static_cast<std::uint8_t>(dim));
            u.setZero();
            u(0) = field == Field::Smooth ? std::sin(Math::Constants::pi() * x(1))
              : field == Field::Quartic   ? x(1) * x(1) * x(1) * x(1)
              : field == Field::Affine    ? x(1)
              : field == Field::Quadratic ? x(1) * x(1)
                                          : x(1) * x(1) * x(1);
            return u;
          });
      }

      auto getPressure() const
      {
        return Variational::RealFunction([field = m_field](const Geometry::Point& x) {
          return field == Field::Smooth ? std::cos(Math::Constants::pi() * x(0))
            : field == Field::Quartic   ? x(0) * x(0) * x(0) - Real(1) / 4
            : field == Field::Cubic     ? x(0) * x(0) - Real(1) / 3
                                        : x(0) - 0.5;
        });
      }

      auto getForcing() const
      {
        return Variational::VectorFunction(
          m_dimension, [dim = m_dimension, field = m_field](const Geometry::Point& x) {
            Math::SpatialVector<Real> f(static_cast<std::uint8_t>(dim));
            f.setZero();
            const Real pi = Math::Constants::pi();
            f(0) = field == Field::Smooth
              ? pi * pi * std::sin(pi * x(1)) - pi * std::sin(pi * x(0))
              : field == Field::Quartic   ? -12 * x(1) * x(1) + 3 * x(0) * x(0)
              : field == Field::Affine    ? 1
              : field == Field::Quadratic ? -1
                                          : -6 * x(1) + 2 * x(0);
            return f;
          });
      }

      auto getVelocityJacobian() const
      {
        return [dim = m_dimension, field = m_field](const Geometry::Point& x) {
          Math::SpatialMatrix<Real> j(
            static_cast<std::uint8_t>(dim), static_cast<std::uint8_t>(dim));
          j.setZero();
          const Real pi = Math::Constants::pi();
          j(0, 1) = field == Field::Smooth ? pi * std::cos(pi * x(1))
            : field == Field::Quartic      ? 4 * x(1) * x(1) * x(1)
            : field == Field::Affine       ? 1
            : field == Field::Quadratic    ? 2 * x(1)
                                           : 3 * x(1) * x(1);
          return j;
        };
      }

      auto getPressureGradient() const
      {
        return Variational::VectorFunction(
          m_dimension, [dim = m_dimension, field = m_field](const Geometry::Point& x) {
            Math::SpatialVector<Real> g(static_cast<std::uint8_t>(dim));
            g.setZero();
            g(0) = field == Field::Smooth
              ? -Math::Constants::pi() * std::sin(Math::Constants::pi() * x(0))
              : field == Field::Quartic ? 3 * x(0) * x(0)
              : field == Field::Cubic   ? 2 * x(0)
                                        : 1;
            return g;
          });
      }

    private:
      size_t m_dimension;
      Field m_field;
  };
}

#endif
