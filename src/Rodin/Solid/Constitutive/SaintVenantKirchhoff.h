/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file SaintVenantKirchhoff.h
 * @brief Saint-Venant-Kirchhoff hyperelastic constitutive law.
 *
 * Implements the stored energy density:
 * @f[
 *   W(\mathbf{E}) = \frac{\lambda}{2}(\operatorname{tr}\mathbf{E})^2
 *                 + \mu \, \mathbf{E} : \mathbf{E}
 * @f]
 * where @f$ \mathbf{E} = \tfrac{1}{2}(\mathbf{F}^T\mathbf{F} - \mathbf{I}) @f$
 * is the Green-Lagrange strain tensor.
 *
 * The second Piola-Kirchhoff stress is:
 * @f[
 *   \mathbf{S} = \lambda (\operatorname{tr}\mathbf{E})\mathbf{I} + 2\mu\mathbf{E}
 * @f]
 *
 * The first Piola-Kirchhoff stress is:
 * @f[
 *   \mathbf{P} = \mathbf{F}\mathbf{S}
 * @f]
 */
#ifndef RODIN_SOLID_CONSTITUTIVE_SAINTVENANTKIRCHHOFF_H
#define RODIN_SOLID_CONSTITUTIVE_SAINTVENANTKIRCHHOFF_H

#include "Rodin/Types.h"
#include "Rodin/Math/Matrix.h"

#include "Rodin/Solid/Kinematics/KinematicState.h"

#include "HyperElasticLaw.h"

namespace Rodin::Solid
{
  /**
   * @brief Saint-Venant-Kirchhoff hyperelastic law.
   *
   * This is the simplest hyperelastic law, extending linear elasticity
   * to large deformations through the Green-Lagrange strain tensor.
   *
   * @see HyperElasticLaw
   */
  class SaintVenantKirchhoff final : public HyperElasticLaw<SaintVenantKirchhoff>
  {
    public:
      /// Precomputed cache for the Saint-Venant-Kirchhoff law.
      struct Cache
      {
        Math::Matrix<Real> E;  ///< Green-Lagrange strain @f$ \mathbf{E} @f$
        Math::Matrix<Real> S;  ///< Second Piola-Kirchhoff stress @f$ \mathbf{S} @f$
        Real trE;              ///< @f$ \operatorname{tr}(\mathbf{E}) @f$
      };

      /**
       * @brief Constructs a Saint-Venant-Kirchhoff law with the given Lamé parameters.
       * @param lambda First Lamé parameter @f$ \lambda @f$
       * @param mu Second Lamé parameter (shear modulus) @f$ \mu @f$
       */
      SaintVenantKirchhoff(Real lambda, Real mu)
        : m_lambda(lambda), m_mu(mu)
      {}

      SaintVenantKirchhoff(const SaintVenantKirchhoff&) = default;
      SaintVenantKirchhoff(SaintVenantKirchhoff&&) = default;

      /// @brief Gets the first Lamé parameter.
      Real getLambda() const { return m_lambda; }

      /// @brief Gets the shear modulus.
      Real getMu() const { return m_mu; }

      void initializeCache(Cache& cache, const KinematicState& state) const
      {
        const auto& C = state.getRightCauchyGreenTensor();
        const size_t d = state.getDimension();
        const auto I = Math::Matrix<Real>::Identity(d, d);

        // E = 0.5 (C - I)
        cache.E = 0.5 * (C - I);
        cache.trE = cache.E.trace();

        // S = lambda tr(E) I + 2 mu E
        cache.S = m_lambda * cache.trE * I + 2.0 * m_mu * cache.E;
      }

      Real strainEnergyDensity(const Cache& cache, const KinematicState&) const
      {
        return 0.5 * m_lambda * cache.trE * cache.trE
             + m_mu * (cache.E.array() * cache.E.array()).sum();
      }

      void firstPiolaKirchhoffStress(
          Math::Matrix<Real>& P,
          const Cache& cache,
          const KinematicState& state) const
      {
        // P = F S
        P = state.getDeformationGradient() * cache.S;
      }

      void materialTangentAction(
          Math::Matrix<Real>& dP,
          const Cache& cache,
          const KinematicState& state,
          const Math::Matrix<Real>& dF) const
      {
        const auto& F = state.getDeformationGradient();
        const size_t d = state.getDimension();
        const auto I = Math::Matrix<Real>::Identity(d, d);

        // dE = 0.5 (dF^T F + F^T dF)
        const Math::Matrix<Real> dE = 0.5 * (dF.transpose() * F + F.transpose() * dF);

        // dS = lambda tr(dE) I + 2 mu dE
        const Math::Matrix<Real> dS = m_lambda * dE.trace() * I + 2.0 * m_mu * dE;

        // dP = dF S + F dS
        dP = dF * cache.S + F * dS;
      }

    private:
      Real m_lambda;
      Real m_mu;
  };
}

#endif
