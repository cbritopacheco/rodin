/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Yeoh.h
 * @brief Compressible Yeoh hyperelastic constitutive law.
 *
 * Implements the stored energy density:
 * @f[
 *   W = c_1(\bar{I}_1 - d) + c_2(\bar{I}_1 - d)^2 + c_3(\bar{I}_1 - d)^3
 *     + \frac{\kappa}{2}(J - 1)^2
 * @f]
 * where the modified first invariant is:
 * @f[
 *   \bar{I}_1 = J^{-2/d} I_1
 * @f]
 * and @f$ d @f$ is the spatial dimension (classically @f$ d = 3 @f$).
 *
 * The Yeoh model is a reduced polynomial model depending only on
 * @f$ \bar{I}_1 @f$. Its cubic form captures the strain-stiffening response
 * typical of soft biological tissue (e.g. arterial wall) and filled rubbers.
 *
 * @note When @f$ c_2 = c_3 = 0 @f$, this reduces to a Neo-Hookean-type model
 *       with shear modulus @f$ \mu = 2 c_1 @f$.
 */
#ifndef RODIN_SOLID_CONSTITUTIVE_YEOH_H
#define RODIN_SOLID_CONSTITUTIVE_YEOH_H

#include <cmath>

#include "Rodin/Types.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SpatialMatrix.h"

#include "Rodin/Solid/Kinematics/KinematicState.h"

#include "HyperElasticLaw.h"

namespace Rodin::Solid
{
  /**
   * @brief Compressible Yeoh hyperelastic law.
   *
   * Parameterized by material constants @f$ c_1, c_2, c_3 @f$ and bulk
   * modulus @f$ \kappa @f$. The small-strain shear modulus is
   * @f$ \mu = 2 c_1 @f$; @f$ c_2, c_3 @f$ control the strain-stiffening.
   *
   * @see HyperElasticLaw, MooneyRivlin, NeoHookean
   */
  class Yeoh final : public HyperElasticLaw<Yeoh>
  {
    public:
      /// Precomputed cache for the Yeoh law.
      struct Cache
      {
        Real I1;      ///< @f$ I_1 = \operatorname{tr}(\mathbf{C}) @f$
        Real J;       ///< Jacobian
        Real Jm2d;    ///< @f$ J^{-2/d} @f$
        Real I1bar;   ///< @f$ \bar{I}_1 = J^{-2/d} I_1 @f$
      };

      /**
       * @brief Constructs a Yeoh law.
       * @param c1 First material constant (small-strain: @f$ \mu = 2 c_1 @f$)
       * @param c2 Second material constant (quadratic stiffening)
       * @param c3 Third material constant (cubic stiffening)
       * @param bulkModulus Bulk modulus @f$ \kappa @f$
       */
      Yeoh(Real c1, Real c2, Real c3, Real bulkModulus)
        : m_c1(c1), m_c2(c2), m_c3(c3), m_kappa(bulkModulus)
      {}

      Yeoh(const Yeoh&) = default;
      Yeoh(Yeoh&&) = default;

      /// @brief Gets @f$ c_1 @f$.
      Real getMaterialConstantC1() const { return m_c1; }

      /// @brief Gets @f$ c_2 @f$.
      Real getMaterialConstantC2() const { return m_c2; }

      /// @brief Gets @f$ c_3 @f$.
      Real getMaterialConstantC3() const { return m_c3; }

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
        cache.Jm2d = std::pow(cache.J, -2.0 / dd);
        cache.I1bar = cache.Jm2d * cache.I1;
      }

      Real getStrainEnergyDensity(const Cache& cache, const ConstitutivePoint& cp) const
      {
        const Real dd = static_cast<Real>(cp.getKinematicState().getDimension());
        const Real x = cache.I1bar - dd;
        return m_c1 * x
             + m_c2 * x * x
             + m_c3 * x * x * x
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
        const size_t d = state.getDimension();
        const Real dd = static_cast<Real>(d);

        // Isochoric scalar response: dW/dI1bar evaluated at x = I1bar - d.
        //   dW1 = c1 + 2 c2 x + 3 c3 x^2
        const Real x = cache.I1bar - dd;
        const Real dW1 = m_c1 + 2.0 * m_c2 * x + 3.0 * m_c3 * x * x;

        // dI1bar/dF = Jm2d (2F - (2/d) I1 F^{-T})
        // SpatialMatrix does not provide direct matrix subtraction operators;
        // written as a linear combination (cf. MooneyRivlin).
        const Math::SpatialMatrix<Real> dI1bar_dF =
          2.0 * cache.Jm2d * F + (-(2.0 / dd) * cache.Jm2d * cache.I1) * FinvT;

        // dJvol/dF = kappa (J - 1) J F^{-T}
        const Math::SpatialMatrix<Real> dJvol_dF =
          m_kappa * (cache.J - 1.0) * cache.J * FinvT;

        P = dW1 * dI1bar_dF + dJvol_dF;
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
        const size_t d = state.getDimension();
        const Real dd = static_cast<Real>(d);

        // Frobenius products needed
        const Real FinvT_dF = FinvT.dot(dF);
        const Real F_dF = F.dot(dF);

        // Derivative of F^{-T} in direction dF: d(F^{-T}) = -F^{-T} dF^T F^{-T}
        const Math::SpatialMatrix<Real> dFinvT = -1.0 * (FinvT * dF.transpose() * FinvT);

        // ---- Scalar variations ----
        const Real dJm2d = -(2.0 / dd) * cache.Jm2d * FinvT_dF;
        const Real dI1 = 2.0 * F_dF;
        // dI1bar = d(Jm2d) I1 + Jm2d dI1   (= dI1bar/dF : dF)
        const Real dI1bar = dJm2d * cache.I1 + cache.Jm2d * dI1;

        // First and second isochoric derivatives at x = I1bar - d:
        //   dW1  = c1 + 2 c2 x + 3 c3 x^2
        //   ddW1 = (2 c2 + 6 c3 x) dI1bar
        const Real x = cache.I1bar - dd;
        const Real dW1 = m_c1 + 2.0 * m_c2 * x + 3.0 * m_c3 * x * x;
        const Real ddW1 = (2.0 * m_c2 + 6.0 * m_c3 * x) * dI1bar;

        // dI1bar/dF = Jm2d (2F - (2/d) I1 F^{-T})
        const Math::SpatialMatrix<Real> dI1bar_dF =
          2.0 * cache.Jm2d * F + (-(2.0 / dd) * cache.Jm2d * cache.I1) * FinvT;

        // Directional derivative of dI1bar/dF (cf. MooneyRivlin):
        // d(Jm2d)/dF term + Jm2d d(2F - (2/d) I1 F^{-T})/dF
        const Math::SpatialMatrix<Real> d_dI1bar_dF =
          dJm2d * (2.0 * F + (-(2.0 / dd) * cache.I1) * FinvT)
          + 2.0 * cache.Jm2d * dF
          + (-(2.0 / dd) * cache.Jm2d * dI1) * FinvT
          + (-(2.0 / dd) * cache.Jm2d * cache.I1) * dFinvT;

        // ---- Volumetric contribution ----
        const Real dJ = cache.J * FinvT_dF;
        const Math::SpatialMatrix<Real> d_dJvol_dF =
          m_kappa * ((dJ * (2.0 * cache.J - 1.0)) * FinvT
                   + (cache.J - 1.0) * cache.J * dFinvT);

        // Chain rule: P_iso = dW1(I1bar) dI1bar/dF, so
        //   dP_iso = ddW1 dI1bar/dF + dW1 d(dI1bar/dF)
        dP = ddW1 * dI1bar_dF + dW1 * d_dI1bar_dF + d_dJvol_dF;
      }

    private:
      Real m_c1;
      Real m_c2;
      Real m_c3;
      Real m_kappa;
  };
}

#endif
