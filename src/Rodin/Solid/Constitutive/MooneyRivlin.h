/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file MooneyRivlin.h
 * @brief Compressible Mooney-Rivlin hyperelastic constitutive law.
 *
 * Implements the stored energy density:
 * @f[
 *   W = c_1(\bar{I}_1 - d) + c_2(\bar{I}_2 - d)
 *     + \frac{\kappa}{2}(J - 1)^2
 * @f]
 * where the modified invariants are:
 * @f[
 *   \bar{I}_1 = J^{-2/d} I_1, \quad
 *   \bar{I}_2 = J^{-4/d} I_2
 * @f]
 *
 * @note When @f$ c_2 = 0 @f$, this reduces to a Neo-Hookean-type model.
 */
#ifndef RODIN_SOLID_CONSTITUTIVE_MOONEYRIVLIN_H
#define RODIN_SOLID_CONSTITUTIVE_MOONEYRIVLIN_H

#include <cmath>

#include "Rodin/Types.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SpatialMatrix.h"

#include "Rodin/Solid/Kinematics/KinematicState.h"

#include "HyperElasticLaw.h"

namespace Rodin::Solid
{
  /**
   * @brief Compressible Mooney-Rivlin hyperelastic law.
   *
   * Parameterized by material constants @f$ c_1, c_2 @f$ and bulk modulus
   * @f$ \kappa @f$.
   *
   * @see HyperElasticLaw
   */
  class MooneyRivlin final : public HyperElasticLaw<MooneyRivlin>
  {
    public:
      /// Precomputed cache for the Mooney-Rivlin law.
      struct Cache
      {
        Real I1;      ///< @f$ I_1 = \operatorname{tr}(\mathbf{C}) @f$
        Real I2;      ///< @f$ I_2 @f$
        Real J;       ///< Jacobian
        Real Jm2d;    ///< @f$ J^{-2/d} @f$
        Real Jm4d;    ///< @f$ J^{-4/d} @f$
        Real I1bar;   ///< @f$ \bar{I}_1 = J^{-2/d} I_1 @f$
        Real I2bar;   ///< @f$ \bar{I}_2 = J^{-4/d} I_2 @f$
      };

      /**
       * @brief Constructs a Mooney-Rivlin law.
       * @param c1 First material constant
       * @param c2 Second material constant
       * @param kappa Bulk modulus @f$ \kappa @f$
       */
      MooneyRivlin(Real c1, Real c2, Real bulkModulus)
        : m_c1(c1), m_c2(c2), m_kappa(bulkModulus)
      {}

      MooneyRivlin(const MooneyRivlin&) = default;
      MooneyRivlin(MooneyRivlin&&) = default;

      /// @brief Gets @f$ c_1 @f$.
      Real getMaterialConstantC1() const { return m_c1; }

      /// @brief Gets @f$ c_2 @f$.
      Real getMaterialConstantC2() const { return m_c2; }

      /// @brief Gets the bulk modulus @f$ \kappa @f$.
      Real getBulkModulus() const { return m_kappa; }

      void setCache(Cache& cache, const ConstitutivePoint& cp) const
      {
        const auto& state = cp.getKinematicState();
        const auto& C = state.getRightCauchyGreenTensor();
        const size_t d = state.getDimension();
        const Real dd = static_cast<Real>(d);

        cache.J = state.getJacobian();
        cache.I1 = C.trace();
        cache.I2 = 0.5 * (cache.I1 * cache.I1 - (C * C).trace());

        cache.Jm2d = std::pow(cache.J, -2.0 / dd);
        cache.Jm4d = std::pow(cache.J, -4.0 / dd);

        cache.I1bar = cache.Jm2d * cache.I1;
        cache.I2bar = cache.Jm4d * cache.I2;
      }

      Real getStrainEnergyDensity(const Cache& cache, const ConstitutivePoint& cp) const
      {
        const Real dd = static_cast<Real>(cp.getKinematicState().getDimension());
        return m_c1 * (cache.I1bar - dd)
             + m_c2 * (cache.I2bar - dd)
             + 0.5 * m_kappa * (cache.J - 1.0) * (cache.J - 1.0);
      }

      void getFirstPiolaKirchhoffStress(
          Math::SpatialMatrix<Real>& P,
          const Cache& cache,
          const ConstitutivePoint& cp) const
      {
        const auto& state = cp.getKinematicState();
        const auto& F = state.getDeformationGradient();
        const auto& FinvT = state.getDeformationGradientInverseTranspose();
        const auto& C = state.getRightCauchyGreenTensor();
        const size_t d = state.getDimension();
        const Real dd = static_cast<Real>(d);

        // dI1bar/dF = Jm2d (2F - (2/d) I1 F^{-T})
        // SpatialMatrix does not provide direct matrix subtraction operators.
        // Original form:
        // dI1bar/dF = Jm2d (2F - (2/d) I1 F^{-T})
        const Math::SpatialMatrix<Real> dI1bar_dF =
          2.0 * cache.Jm2d * F + (-(2.0 / dd) * cache.Jm2d * cache.I1) * FinvT;

        // dI2bar/dF = Jm4d (2 I1 F - 2 F C - (4/d) I2 F^{-T})
        // Written as linear combination to stay within SpatialMatrix operators.
        const Math::SpatialMatrix<Real> dI2bar_dF =
          2.0 * cache.Jm4d * cache.I1 * F
          + (-2.0 * cache.Jm4d) * (F * C)
          + (-(4.0 / dd) * cache.Jm4d * cache.I2) * FinvT;

        // dJvol/dF = kappa (J - 1) J F^{-T}
        const Math::SpatialMatrix<Real> dJvol_dF =
          m_kappa * (cache.J - 1.0) * cache.J * FinvT;

        P = m_c1 * dI1bar_dF + m_c2 * dI2bar_dF + dJvol_dF;
      }

      void getMaterialTangent(
          Math::SpatialMatrix<Real>& dP,
          const Cache& cache,
          const ConstitutivePoint& cp,
          const Math::SpatialMatrix<Real>& dF) const
      {
        const auto& state = cp.getKinematicState();
        const auto& F = state.getDeformationGradient();
        const auto& FinvT = state.getDeformationGradientInverseTranspose();
        const auto& C = state.getRightCauchyGreenTensor();
        const size_t d = state.getDimension();
        const Real dd = static_cast<Real>(d);

        // Compute Frobenius products needed
        const Real FinvT_dF = FinvT.dot(dF);
        const Real F_dF = F.dot(dF);

        // Derivative of F^{-T} in direction dF: d(F^{-T}) = -F^{-T} dF^T F^{-T}
        const Math::SpatialMatrix<Real> dFinvT = -1.0 * (FinvT * dF.transpose() * FinvT);

        // ---- I1bar contribution ----
        // d(Jm2d)/ddF term + Jm2d d(2F - 2/d I1 F^{-T})/ddF
        const Real dJm2d = -(2.0 / dd) * cache.Jm2d * FinvT_dF;
        const Real dI1 = 2.0 * F_dF;
        // Original form:
        // d(Jm2d)/dF term + Jm2d d(2F - (2/d) I1 F^{-T})/dF
        const Math::SpatialMatrix<Real> d_dI1bar_dF =
          dJm2d * (2.0 * F + (-(2.0 / dd) * cache.I1) * FinvT)
          + 2.0 * cache.Jm2d * dF
          + (-(2.0 / dd) * cache.Jm2d * dI1) * FinvT
          + (-(2.0 / dd) * cache.Jm2d * cache.I1) * dFinvT;

        // ---- I2bar contribution ----
        // I2 = 0.5(I1^2 - tr(C^2)) so dI2 = I1 dI1 - tr(C dC) where dC = dF^T F + F^T dF
        const Math::SpatialMatrix<Real> dC = dF.transpose() * F + F.transpose() * dF;
        const Real dI2 = cache.I1 * dI1 - C.dot(dC);

        const Real dJm4d = -(4.0 / dd) * cache.Jm4d * FinvT_dF;
        // Original form:
        // d(Jm4d)/dF term + Jm4d d(2 I1 F - 2 F C - (4/d) I2 F^{-T})/dF
        const Math::SpatialMatrix<Real> d_dI2bar_dF =
          dJm4d * (2.0 * cache.I1 * F + (-2.0) * (F * C) + (-(4.0 / dd) * cache.I2) * FinvT)
          + 2.0 * cache.Jm4d * dI1 * F
          + 2.0 * cache.Jm4d * cache.I1 * dF
          + (-2.0 * cache.Jm4d) * (dF * C)
          + (-2.0 * cache.Jm4d) * (F * dC)
          + (-(4.0 / dd) * cache.Jm4d * dI2) * FinvT
          + (-(4.0 / dd) * cache.Jm4d * cache.I2) * dFinvT;

        // ---- Volumetric contribution ----
        const Real dJ = cache.J * FinvT_dF;
        const Math::SpatialMatrix<Real> d_dJvol_dF =
          m_kappa * ((dJ * (2.0 * cache.J - 1.0)) * FinvT
                   + (cache.J - 1.0) * cache.J * dFinvT);

        dP = m_c1 * d_dI1bar_dF + m_c2 * d_dI2bar_dF + d_dJvol_dF;
      }

    private:
      Real m_c1;
      Real m_c2;
      Real m_kappa;
  };
}

#endif
