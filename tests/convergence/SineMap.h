/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TESTS_CONVERGENCE_SINE_MAP_H
#define RODIN_TESTS_CONVERGENCE_SINE_MAP_H
#include "Convergence.h"
namespace Rodin::Tests::Convergence
{
  /** @brief Analytic reference map, independent of geometry installation.
   * @f$\Phi(\xi)=\xi+a\sin(\pi\xi_0)e_{d-1}@f$ on the unit box.
   * The dimensionless amplitude must satisfy @f$|a|\pi<1@f$ in 1D.
   */
  class SineMap
  {
    public:
      explicit SineMap(Real amplitude = 0.1)
        : m_amplitude(amplitude)
      {}
      Math::SpatialPoint operator()(Math::SpatialPoint xi) const
      {
        xi(xi.size() - 1) += m_amplitude * std::sin(Math::Constants::pi() * xi(0));
        return xi;
      }
      Math::SpatialMatrix<Real> getJacobian(const Math::SpatialPoint& xi) const
      {
        Math::SpatialMatrix<Real> jacobian =
          Math::SpatialMatrix<Real>::Identity(xi.size(), xi.size());
        jacobian(xi.size() - 1, 0) +=
          m_amplitude * Math::Constants::pi() * std::cos(Math::Constants::pi() * xi(0));
        return jacobian;
      }

    private:
      Real m_amplitude;
  };
}
#endif
