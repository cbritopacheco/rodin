/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file HolzapfelOgden.h
 * @brief Compressible Holzapfel-Ogden hyperelastic law with one fiber family.
 */
#ifndef RODIN_SOLID_CONSTITUTIVE_HOLZAPFELOGDEN_H
#define RODIN_SOLID_CONSTITUTIVE_HOLZAPFELOGDEN_H

#include <cmath>

#include "Rodin/Types.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialVector.h"

#include "Rodin/Solid/Kinematics/KinematicState.h"
#include "Rodin/Solid/Local/FiberKinematics.h"
#include "Rodin/Solid/Local/ConstitutivePoint.h"

#include "HyperElasticLaw.h"

namespace Rodin::Solid
{
  /**
   * @brief Compressible Holzapfel-Ogden law with one preferred direction.
   *
   * The stored energy density is
   * @f[
   *   W = \mu_1(\bar I_1 - 3) + \mu_2(\bar I_2 - 3)
   *     + C_0 e^{C_1(\bar I_1 - 3)^2}
   *     + C_2 e^{C_3(\bar I_4 - 1)^2}
   *     + \kappa(J - 1 - \log J),
   * @f]
   * where @f$\bar I_1 = I_1 I_3^{-1/3}@f$,
   * @f$\bar I_2 = I_2 I_3^{-2/3}@f$,
   * @f$\bar I_4 = I_4 I_3^{-1/3}@f$, and @f$I_3=J^2@f$.
   */
  class HolzapfelOgden final : public HyperElasticLaw<HolzapfelOgden>
  {
    public:
      /// @brief Material parameters for the Holzapfel-Ogden law.
      struct Parameters
      {
        /// @brief Isochoric first-invariant coefficient.
        Real mu1 = 0.0;

        /// @brief Isochoric second-invariant coefficient.
        Real mu2 = 0.0;

        /// @brief Exponential first-invariant scale.
        Real C0 = 0.0;

        /// @brief Exponential first-invariant stiffness.
        Real C1 = 0.0;

        /// @brief Fiber exponential scale.
        Real C2 = 0.0;

        /// @brief Fiber exponential stiffness.
        Real C3 = 0.0;

        /// @brief Volumetric penalty coefficient.
        Real kappa = 0.0;
      };

      /// @brief Cached invariant values and fiber direction.
      struct Cache
      {
        /// @brief First invariant of the right Cauchy-Green tensor.
        Real I1 = 0.0;

        /// @brief Second invariant of the right Cauchy-Green tensor.
        Real I2 = 0.0;

        /// @brief Third invariant of the right Cauchy-Green tensor.
        Real I3 = 1.0;

        /// @brief Fiber invariant.
        Real I4 = 1.0;

        /// @brief Jacobian determinant.
        Real J = 1.0;

        /// @brief Cached @f$I_3^{-1/3}@f$.
        Real I3m13 = 1.0;

        /// @brief Cached @f$I_3^{-2/3}@f$.
        Real I3m23 = 1.0;

        /// @brief Isochoric first invariant.
        Real I1bar = 0.0;

        /// @brief Isochoric second invariant.
        Real I2bar = 0.0;

        /// @brief Isochoric fiber invariant.
        Real I4bar = 1.0;

        /// @brief Reference fiber direction.
        Math::SpatialVector<Real> direction;
      };

      /**
       * @brief Constructs the law from a parameter object.
       * @param params Material parameters configuring the law.
       */
      explicit HolzapfelOgden(const Parameters& params)
        : m_params(params)
      {}

      /**
       * @brief Constructs the law from scalar material parameters.
       * @param mu1 Isochoric first-invariant coefficient.
       * @param mu2 Isochoric second-invariant coefficient.
       * @param C0 Exponential first-invariant scale.
       * @param C1 Exponential first-invariant stiffness.
       * @param C2 Fiber exponential scale.
       * @param C3 Fiber exponential stiffness.
       * @param kappa Volumetric penalty coefficient.
       */
      HolzapfelOgden(Real mu1, Real mu2, Real C0, Real C1, Real C2, Real C3, Real kappa)
        : m_params{mu1, mu2, C0, C1, C2, C3, kappa}
      {}

      /**
       * @brief Returns the material parameters.
       * @returns The material parameters.
       */
      const Parameters& getParameters() const
      {
        return m_params;
      }

      /**
       * @brief Populates invariant caches at a constitutive point.
       * @param cache Storage for the constitutive invariant cache.
       * @param cp Constitutive evaluation point.
       */
      void getCache(Cache& cache, const ConstitutivePoint& cp) const
      {
        getCache(cache, cp.getKinematicState(), FiberKinematics::direction(cp));
      }

      /**
       * @brief Returns the stored strain-energy density.
       * @param cache Constitutive invariant cache populated by getCache().
       * @returns The stored strain-energy density.
       * @param cp Constitutive point required by the material-law interface; unused because the energy is obtained from the cache.
       */
      Real getStrainEnergyDensity(
        const Cache& cache, [[maybe_unused]] const ConstitutivePoint& cp) const
      {
        return m_params.mu1 * (cache.I1bar - 3.0) + m_params.mu2 * (cache.I2bar - 3.0) +
          m_params.C0 * std::exp(m_params.C1 * square(cache.I1bar - 3.0)) +
          m_params.C2 * std::exp(m_params.C3 * square(cache.I4bar - 1.0)) +
          m_params.kappa * (cache.J - 1.0 - std::log(cache.J));
      }

      /**
       * @brief Computes the first Piola-Kirchhoff stress.
       * @param P Storage for the first Piola-Kirchhoff stress.
       * @param cache Constitutive invariant cache populated by getCache().
       * @param cp Constitutive evaluation point.
       */
      void getFirstPiolaKirchhoffStress(Math::SpatialMatrix<Real>& P, const Cache& cache,
        const ConstitutivePoint& cp) const
      {
        computeFirstPiolaKirchhoffStress(P, cache, cp.getKinematicState());
      }

      /**
       * @brief Computes the analytic material tangent action at fixed reference fiber direction.
       * @param dP Storage for the resulting stress increment.
       * @param cache Constitutive invariant cache populated by getCache().
       * @param cp Constitutive evaluation point.
       * @param dF Increment of the deformation gradient.
       */
      void getMaterialTangent(Math::SpatialMatrix<Real>& dP, const Cache& cache,
        const ConstitutivePoint& cp, const Math::SpatialMatrix<Real>& dF) const
      {
        const auto& state = cp.getKinematicState();
        const auto& F = state.getDeformationGradient();
        const auto& C = state.getRightCauchyGreenTensor();
        const auto Cinv = C.inverse();
        const auto derivative = getEnergyDerivatives(cache);
        const Math::SpatialMatrix<Real> dC = F.transpose() * dF + dF.transpose() * F;
        const Real trace = (Cinv * dC).trace();
        const Real dI1 = dC.trace();
        const Real dI2 = cache.I1 * dI1 - (C * dC).trace();
        const Real dI4 = cache.direction.dot(dC * cache.direction);
        const Real dI1bar = cache.I3m13 * dI1 - cache.I1bar * trace / 3;
        const Real dI4bar = cache.I3m13 * dI4 - cache.I4bar * trace / 3;
        const Real dW1 =
          -derivative.first * trace / 3 + derivative.firstCurvature * dI1bar;
        const Real dW2 = -2 * derivative.second * trace / 3;
        const Real dW4 =
          -derivative.fiber * trace / 3 + derivative.fiberCurvature * dI4bar;
        const Real dVolume =
          -(dI1 * derivative.first + cache.I1 * dW1 + 2 * dI2 * derivative.second +
            2 * cache.I2 * dW2 + dI4 * derivative.fiber + cache.I4 * dW4) /
            3 +
          m_params.kappa * cache.J * trace / 4;
        const auto identity =
          Math::SpatialMatrix<Real>::Identity(state.getDimension(), state.getDimension());
        const auto fiber = FiberKinematics::dyad(cache.direction);
        const Math::SpatialMatrix<Real> secondInvariant = cache.I1 * identity - C;
        const Math::SpatialMatrix<Real> stress = derivative.first * identity +
          derivative.second * secondInvariant + derivative.volume * Cinv +
          derivative.fiber * fiber;
        const Math::SpatialMatrix<Real> increment = dW1 * identity +
          dW2 * secondInvariant + derivative.second * (dI1 * identity - dC) +
          dVolume * Cinv - derivative.volume * (Cinv * dC * Cinv) + dW4 * fiber;
        dP = 2 * (dF * stress + F * increment);
      }

    private:
      /// @brief Invariant derivatives of the energy and scaled volume derivative.
      struct EnergyDerivatives
      {
        /// @brief Energy derivative with respect to the first invariant.
        Real first;

        /// @brief Energy derivative with respect to the second invariant.
        Real second;

        /// @brief Energy derivative with respect to the fiber invariant.
        Real fiber;

        /// @brief Third-invariant derivative multiplied by that invariant.
        Real volume;

        /// @brief First derivative's change with the normalized first invariant.
        Real firstCurvature;

        /// @brief Fiber derivative's change with the normalized fiber invariant.
        Real fiberCurvature;
      };

      /**
       * @brief Evaluates invariant derivatives shared by stress and its analytic tangent.
       * @param cache Prepared invariants.
       * @returns Derivatives with respect to @f$I_1,I_2,I_4@f$ and @f$I_3W_{,3}@f$.
       *
       * The scaled volume derivative is
       * @f$I_3W_{,3}=-(I_1W_{,1}+2I_2W_{,2}+I_4W_{,4})/3+
       * \kappa(J-1)/2@f$. Its differential and
       * @f$dC^{-1}=-C^{-1}(dC)C^{-1}@f$ give an action linear in @f$dF@f$,
       * without a perturbation step or another constitutive evaluation.
       */
      EnergyDerivatives getEnergyDerivatives(const Cache& cache) const
      {
        const Real firstOffset = cache.I1bar - 3;
        const Real fiberOffset = cache.I4bar - 1;
        const Real firstScale =
          2 * m_params.C0 * m_params.C1 * std::exp(m_params.C1 * square(firstOffset));
        const Real fiberScale =
          2 * m_params.C2 * m_params.C3 * std::exp(m_params.C3 * square(fiberOffset));
        EnergyDerivatives result;
        result.first = cache.I3m13 * (m_params.mu1 + firstScale * firstOffset);
        result.second = m_params.mu2 * cache.I3m23;
        result.fiber = cache.I3m13 * fiberScale * fiberOffset;
        result.volume = -(cache.I1 * result.first + 2 * cache.I2 * result.second +
                          cache.I4 * result.fiber) /
            3 +
          m_params.kappa * (cache.J - 1) / 2;
        result.firstCurvature =
          cache.I3m13 * firstScale * (1 + 2 * m_params.C1 * square(firstOffset));
        result.fiberCurvature =
          cache.I3m13 * fiberScale * (1 + 2 * m_params.C3 * square(fiberOffset));
        return result;
      }

      /**
       * @brief Computes the square of a scalar.
       * @param x Scalar evaluation argument.
       * @returns The scalar multiplied by itself.
       */
      static Real square(Real x)
      {
        return x * x;
      }

      /**
       * @brief Updates the material cache from the kinematic state.
       * @param cache Constitutive cache used or updated by the evaluation.
       * @param state Kinematic state at the material point.
       * @param direction Material direction to include in the stress contribution.
       */

      void getCache(Cache& cache, const KinematicState& state,
        const Math::SpatialVector<Real>& direction) const
      {
        const auto& C = state.getRightCauchyGreenTensor();
        cache.J = state.getJacobian();
        cache.I1 = C.trace();
        cache.I2 = 0.5 * (cache.I1 * cache.I1 - (C * C).trace());
        cache.I3 = cache.J * cache.J;
        cache.I3m13 = std::pow(cache.I3, -1.0 / 3.0);
        cache.I3m23 = cache.I3m13 * cache.I3m13;
        cache.I1bar = cache.I1 * cache.I3m13;
        cache.I2bar = cache.I2 * cache.I3m23;
        cache.direction = direction;
        cache.I4 = direction.dot(C * direction);
        cache.I4bar = cache.I4 * cache.I3m13;
      }

      /**
       * @brief Accumulates the cached first Piola-Kirchhoff stress.
       * @param P First Piola-Kirchhoff stress accumulator.
       * @param cache Constitutive cache used or updated by the evaluation.
       * @param state Kinematic state at the material point.
       */

      void computeFirstPiolaKirchhoffStress(Math::SpatialMatrix<Real>& P,
        const Cache& cache, const KinematicState& state) const
      {
        const auto& F = state.getDeformationGradient();
        const auto& C = state.getRightCauchyGreenTensor();
        const auto& Cinv = C.inverse();
        const size_t d = state.getDimension();

        Math::SpatialMatrix<Real> I(
          static_cast<std::uint8_t>(d), static_cast<std::uint8_t>(d));
        I.setIdentity();
        const Math::SpatialMatrix<Real> A = FiberKinematics::dyad(cache.direction);

        const auto derivative = getEnergyDerivatives(cache);
        const Math::SpatialMatrix<Real> dWdC = derivative.first * I +
          derivative.second * (cache.I1 * I - C) + derivative.volume * Cinv +
          derivative.fiber * A;

        P = 2.0 * F * dWdC;
      }

      Parameters m_params;
  };
}

#endif
