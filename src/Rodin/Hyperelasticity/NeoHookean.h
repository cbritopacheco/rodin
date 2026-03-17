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
 * Implements the stored energy:
 * @f[
 *   W(\mathbf{F}) = \frac{\mu}{2}(I_1 - d) - \mu \ln J
 *                 + \frac{\lambda}{2}(\ln J)^2
 * @f]
 * where @f$ I_1 = \mathrm{tr}(\mathbf{F}^T \mathbf{F}) @f$,
 *       @f$ J = \det \mathbf{F} @f$, and @f$ d @f$ is the spatial dimension.
 */
#ifndef RODIN_HYPERELASTICITY_NEOHOOKEAN_H
#define RODIN_HYPERELASTICITY_NEOHOOKEAN_H

#include <cmath>
#include <cassert>

#include <Eigen/Dense>

#include "Rodin/Types.h"

#include "ConstitutiveLaw.h"

namespace Rodin::Hyperelasticity
{
  /**
   * @brief Compressible Neo-Hookean material.
   *
   * ## Constitutive relations
   *
   * First Piola-Kirchhoff stress:
   * @f[
   *   \mathbf{P} = \mu (\mathbf{F} - \mathbf{F}^{-T})
   *              + \lambda (\ln J)\, \mathbf{F}^{-T}
   * @f]
   *
   * Material tangent @f$ \mathbb{A}_{iJkL} =
   *   \partial P_{iJ}/\partial F_{kL} @f$:
   * @f[
   *   \mathbb{A}_{iJkL} = \mu\,\delta_{ik}\delta_{JL}
   *     + (\mu - \lambda \ln J)\, F^{-1}_{Lk}\, F^{-1}_{Ji}
   *     + \lambda\, F^{-1}_{Ji}\, F^{-1}_{Lk}
   * @f]
   * (see Bonet & Wood, Nonlinear Continuum Mechanics for FEA, 2nd ed.)
   *
   * @tparam (none) – concrete, non-template class
   */
  class NeoHookean final : public ConstitutiveLawBase<NeoHookean>
  {
    public:
      using Parent = ConstitutiveLawBase<NeoHookean>;

      /**
       * @brief Constructs a compressible Neo-Hookean material.
       * @param mu    Shear modulus (second Lamé parameter)
       * @param lambda First Lamé parameter
       */
      NeoHookean(Real mu, Real lambda)
        : m_mu(mu), m_lambda(lambda)
      {}

      NeoHookean(const NeoHookean&) = default;

      NeoHookean* copy() const noexcept override
      {
        return new NeoHookean(*this);
      }

      Real getMu()     const { return m_mu; }
      Real getLambda() const { return m_lambda; }

      /**
       * @brief Computes the first Piola-Kirchhoff stress.
       */
      void stressImpl(const Eigen::Ref<const Eigen::MatrixXd>& F,
                      Eigen::Ref<Eigen::MatrixXd> P) const
      {
        const auto d = F.rows();
        assert(F.cols() == d);
        assert(P.rows() == d && P.cols() == d);

        const Real J = F.determinant();
        assert(J > 0.0);

        const Eigen::MatrixXd Finv = F.inverse();
        const Eigen::MatrixXd FinvT = Finv.transpose();
        const Real lnJ = std::log(J);

        P = m_mu * (F - FinvT) + m_lambda * lnJ * FinvT;
      }

      /**
       * @brief Computes the material tangent modulus.
       */
      void tangentImpl(const Eigen::Ref<const Eigen::MatrixXd>& F,
                       Eigen::Ref<Eigen::MatrixXd> C) const
      {
        const auto d = F.rows();
        assert(F.cols() == d);
        const auto dd = d * d;
        assert(C.rows() == dd && C.cols() == dd);

        const Real J = F.determinant();
        assert(J > 0.0);

        const Eigen::MatrixXd Finv = F.inverse();
        const Real lnJ = std::log(J);

        C.setZero();
        for (Eigen::Index i = 0; i < d; ++i)
          for (Eigen::Index J_ = 0; J_ < d; ++J_)
            for (Eigen::Index k = 0; k < d; ++k)
              for (Eigen::Index L = 0; L < d; ++L)
              {
                Real val = 0.0;
                // mu * delta_ik * delta_JL
                if (i == k && J_ == L)
                  val += m_mu;
                // (mu - lambda*lnJ) * Finv(L,k) * Finv(J,i)  [from FinvT derivative]
                val += (m_mu - m_lambda * lnJ) * Finv(L, k) * Finv(J_, i);
                // lambda * Finv(J,i) * Finv(L,k)  [from lnJ derivative]
                val += m_lambda * Finv(J_, i) * Finv(L, k);
                C(i * d + J_, k * d + L) = val;
              }
      }

      /**
       * @brief Computes the strain energy density.
       */
      Real energyImpl(const Eigen::Ref<const Eigen::MatrixXd>& F) const
      {
        const auto d = F.rows();
        const Real I1 = (F.transpose() * F).trace();
        const Real J = F.determinant();
        assert(J > 0.0);
        const Real lnJ = std::log(J);
        return 0.5 * m_mu * (I1 - d) - m_mu * lnJ
             + 0.5 * m_lambda * lnJ * lnJ;
      }

    private:
      Real m_mu;
      Real m_lambda;
  };
}

#endif
