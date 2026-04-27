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
      struct Parameters
      {
        Real mu1 = 0.0;
        Real mu2 = 0.0;
        Real C0 = 0.0;
        Real C1 = 0.0;
        Real C2 = 0.0;
        Real C3 = 0.0;
        Real kappa = 0.0;
      };

      struct Cache
      {
        Real I1 = 0.0;
        Real I2 = 0.0;
        Real I3 = 1.0;
        Real I4 = 1.0;
        Real J = 1.0;
        Real I3m13 = 1.0;
        Real I3m23 = 1.0;
        Real I1bar = 0.0;
        Real I2bar = 0.0;
        Real I4bar = 1.0;
        Math::SpatialVector<Real> direction;
      };

      explicit HolzapfelOgden(const Parameters& params)
        : m_params(params)
      {}

      HolzapfelOgden(
          Real mu1,
          Real mu2,
          Real C0,
          Real C1,
          Real C2,
          Real C3,
          Real kappa)
        : m_params{mu1, mu2, C0, C1, C2, C3, kappa}
      {}

      const Parameters& getParameters() const
      {
        return m_params;
      }

      void setCache(Cache& cache, const ConstitutivePoint& cp) const
      {
        setCache(cache, cp.getKinematicState(), FiberKinematics::direction(cp));
      }

      Real getStrainEnergyDensity(const Cache& cache, const ConstitutivePoint&) const
      {
        return m_params.mu1 * (cache.I1bar - 3.0)
             + m_params.mu2 * (cache.I2bar - 3.0)
             + m_params.C0 * std::exp(m_params.C1 * square(cache.I1bar - 3.0))
             + m_params.C2 * std::exp(m_params.C3 * square(cache.I4bar - 1.0))
             + m_params.kappa * (cache.J - 1.0 - std::log(cache.J));
      }

      void getFirstPiolaKirchhoffStress(
          Math::SpatialMatrix<Real>& P,
          const Cache& cache,
          const ConstitutivePoint& cp) const
      {
        computeFirstPiolaKirchhoffStress(P, cache, cp.getKinematicState());
      }

      void getMaterialTangent(
          Math::SpatialMatrix<Real>& dP,
          const Cache& cache,
          const ConstitutivePoint& cp,
          const Math::SpatialMatrix<Real>& dF) const
      {
        const auto& state = cp.getKinematicState();
        const auto& H = state.getDisplacementGradient();
        const Real eps = 1.e-7;

        KinematicState plus(state.getDimension());
        plus.setDisplacementGradient(H + eps * dF);
        Cache plusCache = cache;
        setCache(plusCache, plus, cache.direction);
        Math::SpatialMatrix<Real> Pplus;
        computeFirstPiolaKirchhoffStress(Pplus, plusCache, plus);

        Math::SpatialMatrix<Real> P;
        computeFirstPiolaKirchhoffStress(P, cache, state);

        dP = (1.0 / eps) * Pplus + (-1.0 / eps) * P;
      }

    private:
      static Real square(Real x)
      {
        return x * x;
      }

      void setCache(
          Cache& cache,
          const KinematicState& state,
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

      void computeFirstPiolaKirchhoffStress(
          Math::SpatialMatrix<Real>& P,
          const Cache& cache,
          const KinematicState& state) const
      {
        const auto& F = state.getDeformationGradient();
        const auto& C = state.getRightCauchyGreenTensor();
        const auto& Cinv = C.inverse();
        const size_t d = state.getDimension();

        Math::SpatialMatrix<Real> I(static_cast<std::uint8_t>(d), static_cast<std::uint8_t>(d));
        I.setIdentity();
        const Math::SpatialMatrix<Real> A =
          FiberKinematics::dyad(cache.direction);

        const Real W1 =
          2.0 * m_params.C0 * m_params.C1 * cache.I3m13
            * (cache.I1bar - 3.0)
            * std::exp(m_params.C1 * square(cache.I1bar - 3.0))
          + m_params.mu1 * cache.I3m13;
        const Real W2 = m_params.mu2 * cache.I3m23;
        const Real W4 =
          2.0 * m_params.C2 * m_params.C3 * cache.I3m13
            * (cache.I4bar - 1.0)
            * std::exp(m_params.C3 * square(cache.I4bar - 1.0));
        const Real W3 =
          -(1.0 / 3.0) * m_params.mu1 * cache.I1 * std::pow(cache.I3, -4.0 / 3.0)
          -(2.0 / 3.0) * m_params.mu2 * cache.I2 * std::pow(cache.I3, -5.0 / 3.0)
          -(2.0 / 3.0) * m_params.C0 * m_params.C1 * cache.I1
            * std::pow(cache.I3, -4.0 / 3.0) * (cache.I1bar - 3.0)
            * std::exp(m_params.C1 * square(cache.I1bar - 3.0))
          -(1.0 / 3.0) * 2.0 * m_params.C2 * m_params.C3 * cache.I4
            * std::pow(cache.I3, -4.0 / 3.0) * (cache.I4bar - 1.0)
            * std::exp(m_params.C3 * square(cache.I4bar - 1.0))
          + 0.5 * m_params.kappa * (std::pow(cache.I3, -0.5) - 1.0 / cache.I3);

        const Math::SpatialMatrix<Real> dWdC =
          W1 * I
          + W2 * (cache.I1 * I + (-1.0) * C)
          + W3 * cache.I3 * Cinv
          + W4 * A;

        P = 2.0 * F * dWdC;
      }

      Parameters m_params;
  };
}

#endif
