/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file HolzapfelOgden.h
 * @brief Compressible anisotropic Holzapfel-Ogden hyperelastic constitutive law.
 *
 * Stored energy density:
 * @f[
 *   W = \mu_1(\bar I_1 - 3)
 *     + \mu_2(\bar I_2 - 3)
 *     + C_0 e^{C_1(\bar I_1 - 3)^2}
 *     + C_2 e^{C_3(\bar I_4 - 1)^2}
 *     + \kappa (J - 1 - \ln J)
 * @f]
 * where
 * @f[
 *   \bar I_1 = I_1 I_3^{-1/3},\quad
 *   \bar I_2 = I_2 I_3^{-2/3},\quad
 *   \bar I_4 = I_4 I_3^{-1/3},\quad
 *   I_3 = J^2.
 * @f]
 */
#ifndef RODIN_SOLID_CONSTITUTIVE_HOLZAPFELOGDEN_H
#define RODIN_SOLID_CONSTITUTIVE_HOLZAPFELOGDEN_H

#include <cassert>
#include <cmath>

#include "Rodin/Types.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialVector.h"

#include "Rodin/Solid/Kinematics/KinematicState.h"
#include "Rodin/Solid/Local/ConstitutivePoint.h"

#include "HyperElasticLaw.h"

namespace Rodin::Solid
{
  class HolzapfelOgden final : public HyperElasticLaw<HolzapfelOgden>
  {
    public:
      struct Cache
      {
        Real I1;
        Real I2;
        Real I3;
        Real I4;

        Real I3m1d3;
        Real I3m2d3;
        Real I3m4d3;
        Real I3m5d3;
        Real I3m7d3;
        Real I3m8d3;
        Real I3m10d3;

        Real I3p1d3;
        Real I3p2d3;

        Real J;
        Real logJ;

        Real E1;
        Real E4;
      };

      HolzapfelOgden(
          Real kappa,
          Real mu1,
          Real mu2,
          Real c0,
          Real c1,
          Real c2,
          Real c3)
        : m_kappa(kappa),
          m_mu1(mu1),
          m_mu2(mu2),
          m_C0(c0),
          m_C1(c1),
          m_C2(c2),
          m_C3(c3)
      {}

      HolzapfelOgden(const HolzapfelOgden&) = default;
      HolzapfelOgden(HolzapfelOgden&&) = default;

      void setCache(Cache& cache, const ConstitutivePoint& cp) const
      {
        assert(cp.has<Tags::FiberDirection>());
        const auto& state = cp.getKinematicState();
        const auto& C = state.getRightCauchyGreenTensor();
        const auto& a0 = cp.get<Tags::FiberDirection>();

        cache.I1 = C.trace();
        cache.I2 = 0.5 * (cache.I1 * cache.I1 - (C * C).trace());
        cache.I3 = state.getJacobian() * state.getJacobian();
        cache.I4 = a0.dot(C * a0);

        cache.I3m1d3 = std::pow(cache.I3, -1.0 / 3.0);
        cache.I3m2d3 = std::pow(cache.I3, -2.0 / 3.0);
        cache.I3m4d3 = std::pow(cache.I3, -4.0 / 3.0);
        cache.I3m5d3 = std::pow(cache.I3, -5.0 / 3.0);
        cache.I3m7d3 = std::pow(cache.I3, -7.0 / 3.0);
        cache.I3m8d3 = std::pow(cache.I3, -8.0 / 3.0);
        cache.I3m10d3 = std::pow(cache.I3, -10.0 / 3.0);

        cache.I3p1d3 = std::pow(cache.I3, 1.0 / 3.0);
        cache.I3p2d3 = std::pow(cache.I3, 2.0 / 3.0);

        cache.J = state.getJacobian();
        cache.logJ = state.getLogJacobian();

        cache.E1 = std::exp(m_C1 * std::pow(cache.I1 * cache.I3m1d3 - 3.0, 2));
        cache.E4 = std::exp(m_C3 * std::pow(cache.I4 * cache.I3m1d3 - 1.0, 2));
      }

      Real getStrainEnergyDensity(const Cache& cache, const ConstitutivePoint&) const
      {
        return m_mu1 * (cache.I1 * cache.I3m1d3 - 3.0)
             + m_mu2 * (cache.I2 * cache.I3m2d3 - 3.0)
             + m_C0 * (cache.E1 - 1.0)
             + m_C2 * (cache.E4 - 1.0)
             + m_kappa * (cache.J - 1.0 - cache.logJ);
      }

      void getFirstPiolaKirchhoffStress(
          Math::SpatialMatrix<Real>& P,
          const Cache& cache,
          const ConstitutivePoint& cp) const
      {
        const auto& state = cp.getKinematicState();
        const auto& F = state.getDeformationGradient();
        const auto& C = state.getRightCauchyGreenTensor();
        const auto& FinvT = state.getDeformationGradientInverseTranspose();
        const auto& a0 = cp.get<Tags::FiberDirection>();

        const auto w1 = dWdI1(cache);
        const auto w2 = dWdI2(cache);
        const auto w3 = dWdI3(cache);
        const auto w4 = dWdI4(cache);

        const auto dI1_dF = 2.0 * F;
        const auto dI2_dF = 2.0 * (cache.I1 * F + (-1.0) * (F * C));
        const auto dI3_dF = 2.0 * cache.I3 * FinvT;

        const Math::SpatialVector<Real> Fa0 = F * a0;
        Math::SpatialMatrix<Real> dI4_dF(F.rows(), F.cols());
        for (size_t i = 0; i < F.rows(); ++i)
          for (size_t j = 0; j < F.cols(); ++j)
            dI4_dF(i, j) = 2.0 * Fa0[i] * a0[j];

        P = w1 * dI1_dF + w2 * dI2_dF + w3 * dI3_dF + w4 * dI4_dF;
      }

      void getMaterialTangent(
          Math::SpatialMatrix<Real>& dP,
          const Cache& cache,
          const ConstitutivePoint& cp,
          const Math::SpatialMatrix<Real>& dF) const
      {
        const auto& state = cp.getKinematicState();
        const auto& F = state.getDeformationGradient();
        const auto& C = state.getRightCauchyGreenTensor();
        const auto& FinvT = state.getDeformationGradientInverseTranspose();
        const auto& a0 = cp.get<Tags::FiberDirection>();

        const auto dC = dF.transpose() * F + F.transpose() * dF;
        const auto dFinvT = -1.0 * (FinvT * dF.transpose() * FinvT);

        const auto A1 = 2.0 * F;
        const auto A2 = 2.0 * (cache.I1 * F + (-1.0) * (F * C));
        const auto A3 = 2.0 * cache.I3 * FinvT;
        const Math::SpatialVector<Real> Fa0 = F * a0;
        Math::SpatialMatrix<Real> A4(F.rows(), F.cols());
        for (size_t i = 0; i < F.rows(); ++i)
          for (size_t j = 0; j < F.cols(); ++j)
            A4(i, j) = 2.0 * Fa0[i] * a0[j];

        const Real dI1 = A1.dot(dF);
        const Real dI2 = A2.dot(dF);
        const Real dI3 = A3.dot(dF);
        const Real dI4 = A4.dot(dF);

        const auto dA1 = 2.0 * dF;
        const auto dA2 = 2.0 * (dI1 * F + cache.I1 * dF + (-1.0) * (dF * C) + (-1.0) * (F * dC));
        const auto dA3 = 2.0 * (dI3 * FinvT + cache.I3 * dFinvT);
        const Math::SpatialVector<Real> dFa0 = dF * a0;
        Math::SpatialMatrix<Real> dA4(F.rows(), F.cols());
        for (size_t i = 0; i < F.rows(); ++i)
          for (size_t j = 0; j < F.cols(); ++j)
            dA4(i, j) = 2.0 * dFa0[i] * a0[j];

        const auto w1 = dWdI1(cache);
        const auto w2 = dWdI2(cache);
        const auto w3 = dWdI3(cache);
        const auto w4 = dWdI4(cache);

        const auto w11 = d2WdI1I1(cache);
        const auto w22 = d2WdI2I2(cache);
        const auto w33 = d2WdI3I3(cache);
        const auto w44 = d2WdI4I4(cache);
        const auto w12 = d2WdI1I2(cache);
        const auto w13 = d2WdI1I3(cache);
        const auto w14 = d2WdI1I4(cache);
        const auto w23 = d2WdI2I3(cache);
        const auto w24 = d2WdI2I4(cache);
        const auto w34 = d2WdI3I4(cache);

        dP = w1 * dA1 + w2 * dA2 + w3 * dA3 + w4 * dA4
           + (w11 * dI1 + w12 * dI2 + w13 * dI3 + w14 * dI4) * A1
           + (w12 * dI1 + w22 * dI2 + w23 * dI3 + w24 * dI4) * A2
           + (w13 * dI1 + w23 * dI2 + w33 * dI3 + w34 * dI4) * A3
           + (w14 * dI1 + w24 * dI2 + w34 * dI3 + w44 * dI4) * A4;
      }

    private:
      Real dWdI1(const Cache& c) const
      {
        return 2.0 * m_C0 * m_C1 * c.I3m2d3 * (c.I1 - 3.0 * c.I3p1d3) * c.E1
             + m_mu1 * c.I3m1d3;
      }

      Real dWdI2(const Cache& c) const
      {
        return m_mu2 * c.I3m2d3;
      }

      Real dWdI3(const Cache& c) const
      {
        return -2.0 / 3.0 * m_C0 * m_C1 * c.I3m5d3 * c.I1 * (c.I1 - 3.0 * c.I3p1d3) * c.E1
             + 2.0 / 3.0 * m_C2 * m_C3 * c.I3m5d3 * c.I4 * (c.I3p1d3 - c.I4) * c.E4
             - 1.0 / 3.0 * m_mu1 * c.I1 * c.I3m4d3
             - 2.0 / 3.0 * m_mu2 * c.I2 * c.I3m5d3
             + 0.5 * m_kappa * (std::pow(c.I3, -0.5) - 1.0 / c.I3);
      }

      Real dWdI4(const Cache& c) const
      {
        return -2.0 * m_C2 * m_C3 * c.I3m2d3 * (c.I3p1d3 - c.I4) * c.E4;
      }

      Real d2WdI1I1(const Cache& c) const
      {
        return 2.0 * m_C0 * m_C1 * c.I3m4d3
             * (2.0 * m_C1 * std::pow(c.I1 - 3.0 * c.I3p1d3, 2) + c.I3p2d3)
             * c.E1;
      }

      Real d2WdI2I2(const Cache&) const
      {
        return 0.0;
      }

      Real d2WdI3I3(const Cache& c) const
      {
        return 2.0 / 9.0 * c.I3m10d3
             * (m_C0 * m_C1 * c.I1
                * (2.0 * m_C1 * c.I1 * std::pow(c.I1 - 3.0 * c.I3p1d3, 2)
                   + 5.0 * c.I1 * c.I3p2d3 - 12.0 * c.I3) * c.E1
                + m_C2 * m_C3 * c.I4
                * (-4.0 * c.I3 + (5.0 + 2.0 * m_C3) * c.I4 * c.I3p2d3
                   - 4.0 * m_C3 * c.I3p1d3 * c.I4 * c.I4
                   + 2.0 * m_C3 * c.I4 * c.I4 * c.I4) * c.E4)
             + m_mu1 * 4.0 / 9.0 * c.I1 * c.I3m7d3
             + m_mu2 * 10.0 / 9.0 * c.I2 * c.I3m8d3
             + m_kappa * (0.5 * std::pow(c.I3, -2.0) - 0.25 * std::pow(c.I3, -1.5));
      }

      Real d2WdI4I4(const Cache& c) const
      {
        return 2.0 * m_C2 * m_C3 * c.I3m4d3
             * ((1.0 + 2.0 * m_C3) * c.I3p2d3 - 4.0 * m_C3 * c.I3p1d3 * c.I4
                + 2.0 * m_C3 * c.I4 * c.I4) * c.E4;
      }

      Real d2WdI1I2(const Cache&) const
      {
        return 0.0;
      }

      Real d2WdI1I3(const Cache& c) const
      {
        return -2.0 / 3.0 * m_C0 * m_C1 * c.I3m7d3
             * (2.0 * m_C1 * c.I1 * std::pow(c.I1 - 3.0 * c.I3p1d3, 2)
                + 2.0 * c.I1 * c.I3p2d3 - 3.0 * c.I3) * c.E1
             - m_mu1 * 1.0 / 3.0 * c.I3m4d3;
      }

      Real d2WdI1I4(const Cache&) const
      {
        return 0.0;
      }

      Real d2WdI2I3(const Cache& c) const
      {
        return -2.0 / 3.0 * m_mu2 * c.I3m5d3;
      }

      Real d2WdI2I4(const Cache&) const
      {
        return 0.0;
      }

      Real d2WdI3I4(const Cache& c) const
      {
        return -2.0 / 3.0 * m_C2 * m_C3 * c.I3m7d3
             * (-c.I3 + 2.0 * (1.0 + m_C3) * c.I3p2d3 * c.I4
                - 4.0 * m_C3 * c.I3p1d3 * c.I4 * c.I4
                + 2.0 * m_C3 * c.I4 * c.I4 * c.I4) * c.E4;
      }

      Real m_kappa;
      Real m_mu1;
      Real m_mu2;
      Real m_C0;
      Real m_C1;
      Real m_C2;
      Real m_C3;
  };
}

#endif
