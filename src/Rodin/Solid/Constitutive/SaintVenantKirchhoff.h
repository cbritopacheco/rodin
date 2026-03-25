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
#include "Rodin/Math/SpatialMatrix.h"

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
        Math::SpatialMatrix<Real> E;  ///< Green-Lagrange strain @f$ \mathbf{E} @f$
        Math::SpatialMatrix<Real> S;  ///< Second Piola-Kirchhoff stress @f$ \mathbf{S} @f$
        Real trE;              ///< @f$ \operatorname{tr}(\mathbf{E}) @f$
      };

      /**
       * @brief Constructs a Saint-Venant-Kirchhoff law with the given Lamé parameters.
       * @param lambda First Lamé parameter @f$ \lambda @f$
       * @param mu Second Lamé parameter (shear modulus) @f$ \mu @f$
       */
      SaintVenantKirchhoff(Real lameFirstParameter, Real shearModulus)
        : m_lambda(lameFirstParameter), m_mu(shearModulus)
      {}

      SaintVenantKirchhoff(const SaintVenantKirchhoff&) = default;
      SaintVenantKirchhoff(SaintVenantKirchhoff&&) = default;

      /// @brief Gets the first Lamé parameter.
      Real getLameFirstParameter() const { return m_lambda; }

      /// @brief Gets the shear modulus.
      Real getShearModulus() const { return m_mu; }

      void setCache(Cache& cache, const ConstitutivePoint& cp) const
      {
        const auto& state = cp.getKinematicState();
        const auto& C = state.getRightCauchyGreenTensor();
        const size_t d = state.getDimension();
        Math::SpatialMatrix<Real> I;
        I.resize(static_cast<std::uint8_t>(d), static_cast<std::uint8_t>(d));
        I.setIdentity();

        // E = 0.5 (C - I)
        cache.E = 0.5 * C + (-0.5) * I;
        cache.trE = cache.E.trace();

        // S = lambda tr(E) I + 2 mu E
        cache.S = m_lambda * cache.trE * I + 2.0 * m_mu * cache.E;
      }

      Real getStrainEnergyDensity(const Cache& cache, const ConstitutivePoint&) const
      {
        return 0.5 * m_lambda * cache.trE * cache.trE
             + m_mu * cache.E.dot(cache.E);
      }

      void getFirstPiolaKirchhoffStress(
          Math::SpatialMatrix<Real>& P,
          const Cache& cache,
          const ConstitutivePoint& cp) const
      {
        // P = F S
        P = cp.getKinematicState().getDeformationGradient() * cache.S;
      }

      void getMaterialTangent(
          Math::SpatialMatrix<Real>& dP,
          const Cache& cache,
          const ConstitutivePoint& cp,
          const Math::SpatialMatrix<Real>& dF) const
      {
        const auto& state = cp.getKinematicState();
        const auto& F = state.getDeformationGradient();
        const size_t d = state.getDimension();
        Math::SpatialMatrix<Real> I;
        I.resize(static_cast<std::uint8_t>(d), static_cast<std::uint8_t>(d));
        I.setIdentity();

        // dE = 0.5 (dF^T F + F^T dF)
        const Math::SpatialMatrix<Real> dE = 0.5 * (dF.transpose() * F + F.transpose() * dF);

        // dS = lambda tr(dE) I + 2 mu dE
        const Math::SpatialMatrix<Real> dS = m_lambda * dE.trace() * I + 2.0 * m_mu * dE;

        // dP = dF S + F dS
        dP = dF * cache.S + F * dS;
      }

    private:
      Real m_lambda;
      Real m_mu;
  };
}

#endif
