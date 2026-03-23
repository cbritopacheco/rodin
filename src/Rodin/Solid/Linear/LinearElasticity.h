/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file LinearElasticity.h
 * @brief Linear elasticity constitutive model within the Solid module.
 *
 * Provides a self-contained linear elasticity class that wraps Hooke's law
 * for computing the Cauchy stress from infinitesimal strain.
 */
#ifndef RODIN_SOLID_LINEAR_LINEARELASTICITY_H
#define RODIN_SOLID_LINEAR_LINEARELASTICITY_H

#include "Rodin/Types.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SpatialMatrix.h"

#include "Hooke.h"

namespace Rodin::Solid
{
  /**
   * @brief Linear elasticity constitutive model.
   *
   * Wraps Hooke's law to provide stress computation from infinitesimal strain:
   * @f[
   *   \boldsymbol{\sigma} = \lambda \operatorname{tr}(\boldsymbol{\varepsilon}) \mathbf{I}
   *                       + 2\mu \boldsymbol{\varepsilon}
   * @f]
   *
   * ## Usage
   * @code
   * auto le = LinearElasticity::fromYoungPoisson(200.0, 0.3);
    * Math::SpatialMatrix<Real> sigma;
   * le.getStress(sigma, epsilon);
   * @endcode
   */
  class LinearElasticity
  {
    public:
      /**
       * @brief Constructs from Lamé parameters.
       * @param lambda First Lamé parameter
       * @param mu Shear modulus
       */
      LinearElasticity(Real lambda, Real mu)
        : m_hooke(lambda, mu)
      {}

      /**
       * @brief Constructs from an existing Hooke law.
       * @param hooke The Hooke law object
       */
      LinearElasticity(const Hooke& hooke)
        : m_hooke(hooke)
      {}

      LinearElasticity(const LinearElasticity&) = default;
      LinearElasticity(LinearElasticity&&) = default;

      /// @brief Gets the underlying Hooke law.
      const Hooke& getHooke() const { return m_hooke; }

      /// @brief Gets the first Lamé parameter.
      Real getLambda() const { return m_hooke.getLambda(); }

      /// @brief Gets the shear modulus.
      Real getMu() const { return m_hooke.getMu(); }

      /**
       * @brief Computes the Cauchy stress from infinitesimal strain.
       * @param[out] sigma Output stress tensor
       * @param[in] epsilon Infinitesimal strain tensor
       */
      void getStress(Math::SpatialMatrix<Real>& sigma, const Math::SpatialMatrix<Real>& epsilon) const
      {
        m_hooke.stress(sigma, epsilon);
      }

      /**
       * @brief Creates a LinearElasticity from Young's modulus and Poisson's ratio.
       * @param E Young's modulus
       * @param nu Poisson's ratio
       */
      static LinearElasticity fromYoungPoisson(Real E, Real nu)
      {
        return LinearElasticity(Hooke::fromYoungPoisson(E, nu));
      }

    private:
      Hooke m_hooke;
  };
}

#endif
