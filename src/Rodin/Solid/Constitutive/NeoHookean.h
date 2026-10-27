/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file NeoHookean.h
 * @brief Compressible Neo-Hookean hyperelastic constitutive law.
 *
 * Implements the stored energy density:
 * @f[
 *   W(\mathbf{F}) = \frac{\mu}{2}(I_1 - d) - \mu \ln J
 *                 + \frac{\lambda}{2}(\ln J)^2
 * @f]
 * where @f$ I_1 = \operatorname{tr}(\mathbf{F}^T \mathbf{F}) @f$,
 * @f$ J = \det(\mathbf{F}) @f$, and @f$ d @f$ is the spatial dimension.
 *
 * The first Piola-Kirchhoff stress is:
 * @f[
 *   \mathbf{P} = \mu(\mathbf{F} - \mathbf{F}^{-T})
 *              + \lambda \ln(J) \, \mathbf{F}^{-T}
 * @f]
 *
 * The material tangent action on @f$ \delta\mathbf{F} @f$ is:
 * @f[
 *   D\mathbf{P}[\delta\mathbf{F}]
 *     = \mu \, \delta\mathbf{F}
 *     + \lambda (\mathbf{F}^{-T} : \delta\mathbf{F}) \, \mathbf{F}^{-T}
 *     + (\mu - \lambda \ln J) \, \mathbf{F}^{-T} \delta\mathbf{F}^T \mathbf{F}^{-T}
 * @f]
 */
#ifndef RODIN_SOLID_CONSTITUTIVE_NEOHOOKEAN_H
#define RODIN_SOLID_CONSTITUTIVE_NEOHOOKEAN_H

#include <cmath>

#include "Rodin/Types.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SpatialMatrix.h"

#include "Rodin/Solid/Kinematics/KinematicState.h"

#include "HyperElasticLaw.h"

namespace Rodin::Solid
{
  /**
   * @brief Compressible Neo-Hookean hyperelastic law.
   *
   * @see <a href="class_rodin_1_1_solid_1_1_hyper_elastic_law.html">HyperElasticLaw</a>
   */
  class NeoHookean final : public HyperElasticLaw<NeoHookean>
  {
    public:
      /**
       * @brief Constructs a Neo-Hookean law with the given Lamé parameters.
       * @param lambda First Lamé parameter @f$ \lambda @f$
       * @param mu Second Lamé parameter (shear modulus) @f$ \mu @f$
       */
      NeoHookean(Real lambda, Real mu)
        : m_lambda(lambda),
          m_mu(mu)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      NeoHookean(const NeoHookean& other) = default;

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      NeoHookean(NeoHookean&& other) = default;

      /**
       * @brief Gets the first Lamé parameter.
       * @returns The first Lamé parameter.
       */
      Real getLameFirstParameter() const
      {
        return m_lambda;
      }

      /**
       * @brief Gets the shear modulus.
       * @returns The shear modulus.
       */
      Real getShearModulus() const
      {
        return m_mu;
      }

      /**
       * @brief Binds owned evaluation data, reusing unchanged constitutive inputs.
       * @param cp Constitutive point supplying kinematics and law-specific inputs.
       * @returns This law after successful binding.
       */
      NeoHookean& setConstitutivePoint(const ConstitutivePoint& cp)
      {
        if (m_state && hasSameKinematics(*m_state, cp.getKinematicState()))
        {
          return *this;
        }
        m_state.reset();
        computeCache(m_cache, cp);
        m_state.emplace(cp.getKinematicState());
        return *this;
      }

      /**
       * @brief Returns the owned kinematics of the prepared point.
       * @returns Prepared kinematic state.
       */
      const KinematicState& getKinematicState() const
      {
        this->require(m_state.has_value(), __func__);
        return *m_state;
      }

      /**
       * @brief Returns the strain-energy density at the prepared point.
       * @returns Stored strain-energy density.
       */
      Real getStrainEnergyDensity() const
      {
        const auto& cache = m_cache;
        const auto& state = getKinematicState();
        const size_t d = state.getDimension();
        return 0.5 * m_mu * (cache.I1 - static_cast<Real>(d)) - m_mu * cache.logJ +
          0.5 * m_lambda * cache.logJ * cache.logJ;
      }

      /**
       * @brief Computes stress at the prepared point.
       * @param P Output first Piola-Kirchhoff stress.
       */
      void getFirstPiolaKirchhoffStress(Math::SpatialMatrix<Real>& P) const
      {
        const auto& cache = m_cache;
        const auto& state = getKinematicState();
        const auto& F = state.getDeformationGradient();
        const auto& FinvT = state.getDeformationGradientInverseTranspose();

        // P = mu * F + (lambda ln(J) - mu) * F^{-T}
        P = m_mu * F + (m_lambda * cache.logJ - m_mu) * FinvT;
      }

      /**
       * @brief Computes the tangent action at the prepared point.
       * @param dP Output stress increment.
       * @param dF Deformation-gradient increment.
       */
      void getMaterialTangent(
        Math::SpatialMatrix<Real>& dP, const Math::SpatialMatrix<Real>& dF) const
      {
        const auto& cache = m_cache;
        const auto& state = getKinematicState();
        const auto& FinvT = state.getDeformationGradientInverseTranspose();

        // Frobenius inner product F^{-T} : dF
        const Real FinvT_dF = FinvT.dot(dF);

        // dP = mu dF + lambda (F^{-T} : dF) F^{-T}
        //    + (mu - lambda ln J) F^{-T} dF^T F^{-T}
        dP = m_mu * dF + m_lambda * FinvT_dF * FinvT +
          (m_mu - m_lambda * cache.logJ) * FinvT * dF.transpose() * FinvT;
      }

    private:
      /// Precomputed cache for the Neo-Hookean law.
      struct Cache
      {
        Real I1;    ///< First invariant @f$ I_1 = \operatorname{tr}(\mathbf{C}) @f$
        Real J;     ///< Jacobian
        Real logJ;  ///< @f$ \ln(J) @f$
      };

      /**
       * @brief Populates the invariant cache at a constitutive point.
       * @param cache Storage for the constitutive invariant cache.
       * @param cp Constitutive evaluation point.
       */
      void computeCache(Cache& cache, const ConstitutivePoint& cp) const
      {
        const auto& state = cp.getKinematicState();
        const auto& C = state.getRightCauchyGreenTensor();
        cache.I1 = C.trace();
        cache.J = state.getJacobian();
        cache.logJ = state.getLogJacobian();
      }

      /// @brief Private intermediate values for the prepared input.
      Cache m_cache{};

      /// @brief Owned snapshot; empty until binding succeeds.
      Optional<KinematicState> m_state;

      Real m_lambda;
      Real m_mu;
  };
}

#endif
