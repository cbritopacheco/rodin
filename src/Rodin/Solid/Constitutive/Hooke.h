/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Constitutive/Hooke.h
 * @brief Hooke's law for linear elasticity.
 *
 * Provides the isotropic elasticity tensor @f$ \mathbb{C} @f$ such that:
 * @f[
 *   \boldsymbol{\sigma} = \mathbb{C} : \boldsymbol{\varepsilon}
 *                       = \lambda (\operatorname{tr}\boldsymbol{\varepsilon})\mathbf{I}
 *                       + 2\mu \boldsymbol{\varepsilon}
 * @f]
 * where @f$ \boldsymbol{\varepsilon} = \tfrac{1}{2}(\nabla\mathbf{u}
 * + (\nabla\mathbf{u})^T) @f$ is the infinitesimal strain tensor.
 */
#ifndef RODIN_SOLID_CONSTITUTIVE_HOOKE_H
#define RODIN_SOLID_CONSTITUTIVE_HOOKE_H

#include "Rodin/Types.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SpatialMatrix.h"

namespace Rodin::Solid
{
  /**
   * @brief Isotropic Hooke's law for linear elasticity.
   *
   * Given the Lamé parameters @f$ \lambda @f$ and @f$ \mu @f$, computes:
   * @f[
   *   \boldsymbol{\sigma} = \lambda (\operatorname{tr}\boldsymbol{\varepsilon})\mathbf{I}
   *                       + 2\mu \boldsymbol{\varepsilon}
   * @f]
   *
   * ## Material parameter conversions
   * From Young's modulus @f$ E @f$ and Poisson's ratio @f$ \nu @f$:
   * @f[
   *   \lambda = \frac{E\nu}{(1+\nu)(1-2\nu)}, \quad
   *   \mu = \frac{E}{2(1+\nu)}
   * @f]
   */
  class Hooke
  {
    public:
      /**
       * @brief Constructs Hooke's law with the given elastic parameters.
       * @param lambda First Lamé parameter @f$ \lambda @f$
       * @param mu Second Lamé parameter (shear modulus) @f$ \mu @f$
       */
      Hooke(Real lameFirstParameter, Real shearModulus)
        : m_lambda(lameFirstParameter), m_mu(shearModulus)
      {}

      Hooke(const Hooke&) = default;
      Hooke(Hooke&&) = default;

      /// @brief Gets the first Lamé parameter.
      Real getLameFirstParameter() const { return m_lambda; }

      /// @brief Gets the shear modulus.
      Real getShearModulus() const { return m_mu; }

      /**
       * @brief Computes the Cauchy stress from the infinitesimal strain.
       * @param[out] sigma Output stress tensor
       * @param[in] epsilon Infinitesimal strain tensor @f$ \boldsymbol{\varepsilon} @f$
       */
      void getStress(Math::SpatialMatrix<Real>& sigma, const Math::SpatialMatrix<Real>& epsilon) const
      {
        const size_t d = epsilon.rows();
        Math::SpatialMatrix<Real> I;
        I.resize(static_cast<std::uint8_t>(d), static_cast<std::uint8_t>(d));
        I.setIdentity();
        sigma = m_lambda * epsilon.trace()
              * I
              + 2.0 * m_mu * epsilon;
      }

      /**
       * @brief Constructs Hooke's law from Young's modulus and Poisson's ratio.
       * @param E Young's modulus
       * @param nu Poisson's ratio
       * @returns Hooke object with computed Lamé parameters
       */
      static Hooke YoungPoisson(Real E, Real nu)
      {
        const Real lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        const Real mu = E / (2.0 * (1.0 + nu));
        return Hooke(lambda, mu);
      }

    private:
      Real m_lambda;
      Real m_mu;
  };
}

#endif
