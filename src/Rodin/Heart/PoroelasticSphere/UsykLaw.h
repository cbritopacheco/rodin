/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file UsykLaw.h
 * @brief Usyk orthotropic exponential skeletal energy reduced to the radial
 * deformation of the spherical wall (Ruz et al., coronary FSI, Sec. 4.3).
 */
#ifndef RODIN_HEART_POROELASTICSPHERE_USYKLAW_H
#define RODIN_HEART_POROELASTICSPHERE_USYKLAW_H

#include <cmath>

#include "Rodin/Heart/PoroelasticSphere/PassiveLaw.h"

namespace Rodin::Heart
{
  /**
   * @brief Reduced Usyk passive energy on @f$ \mathbf F = \mathrm{diag}(\lambda_r,
   * \lambda_\theta, \lambda_\theta) @f$.
   *
   * @f$ \Psi_M = \tfrac{c_0}{2}(e^{Q} - 1) @f$,
   * @f$ Q = b_{rr}E_{RR}^2 + 2 b_{\theta\theta}E_{\theta\theta}^2 + 2 b_{r\theta}E_{RR}E_{\theta\theta} @f$,
   * @f$ E_{RR} = \tfrac12(\lambda_r^2-1) @f$, @f$ E_{\theta\theta} = \tfrac12(\lambda_\theta^2-1) @f$.
   *
   * The radial Piola components follow the convention of @ref RadialPiolaStress:
   * @f$ P_{M,R} = \partial_{\lambda_r}\Psi_M = \tfrac{c_0}{2}e^{Q}\lambda_r[b_{rr}(\lambda_r^2-1) + b_{r\theta}(\lambda_\theta^2-1)] @f$,
   * @f$ P_{M,\theta} = \tfrac12\partial_{\lambda_\theta}\Psi_M = \tfrac{c_0}{4}e^{Q}\lambda_\theta[2b_{\theta\theta}(\lambda_\theta^2-1) + b_{r\theta}(\lambda_r^2-1)] @f$.
   *
   * @note The manuscript's Eq. (92) carries @f$ c_0/2 @f$ in @f$ P_{M,\theta} @f$;
   * the factor consistent with @f$ \Psi_M @f$ (and verified by finite
   * differences in the unit tests) is @f$ c_0/4 @f$, used here.
   *
   * @tparam Scalar Scalar numeric type.
   */
  template <class Scalar>
  class UsykRadialLaw
  {
    public:
      /// @brief Usyk parameters (stress scale and dimensionless exponents).
      struct Parameters
      {
        Scalar c0 = Scalar(0);  ///< Stress scale @f$ c_0 @f$.
        Scalar brr = Scalar(0); ///< Transmural (radial) exponent @f$ b_{rr} @f$.
        Scalar btt = Scalar(0); ///< In-plane exponent @f$ b_{\theta\theta} @f$.
        Scalar brt = Scalar(0); ///< Radial-circumferential coupling @f$ b_{r\theta} @f$.
      };

      explicit UsykRadialLaw(const Parameters& params = Parameters())
        : m_params(params)
      {}

      /// @brief Skeletal energy @f$ \Psi_M(\lambda_r, \lambda_\theta) @f$.
      Scalar energy(const Scalar lr, const Scalar lt) const
      {
        return Scalar(0.5) * m_params.c0 * (std::exp(exponent(lr, lt)) - Scalar(1));
      }

      /// @brief Radial Piola components and their tangent.
      RadialPiolaStress<Scalar> evaluate(const Scalar lr, const Scalar lt) const
      {
        const Scalar A = lr * lr - Scalar(1);
        const Scalar B = lt * lt - Scalar(1);
        const Scalar brr = m_params.brr, btt = m_params.btt, brt = m_params.brt;

        const Scalar Q = exponent(lr, lt);
        const Scalar Qr = lr * (brr * A + brt * B);
        const Scalar Qt = lt * (Scalar(2) * btt * B + brt * A);
        const Scalar Qrr = brr * A + brt * B + Scalar(2) * brr * lr * lr;
        const Scalar Qrt = Scalar(2) * brt * lr * lt;
        const Scalar Qtt = Scalar(2) * btt * B + brt * A + Scalar(4) * btt * lt * lt;

        const Scalar cr = Scalar(0.5) * m_params.c0 * std::exp(Q);
        const Scalar ct = Scalar(0.5) * cr;

        RadialPiolaStress<Scalar> out;
        out.PRR = cr * Qr;
        out.Ptt = ct * Qt;
        out.dPRR_dlr = cr * (Qr * Qr + Qrr);
        out.dPRR_dlt = cr * (Qr * Qt + Qrt);
        out.dPtt_dlr = ct * (Qt * Qr + Qrt);
        out.dPtt_dlt = ct * (Qt * Qt + Qtt);
        return out;
      }

      const Parameters& getParameters() const noexcept
      {
        return m_params;
      }

    private:
      Scalar exponent(const Scalar lr, const Scalar lt) const
      {
        const Scalar A = lr * lr - Scalar(1);
        const Scalar B = lt * lt - Scalar(1);
        return Scalar(0.25) *
          (m_params.brr * A * A + Scalar(2) * m_params.btt * B * B +
            Scalar(2) * m_params.brt * A * B);
      }

      Parameters m_params;
  };

  /**
   * @brief Passive stress operator plugging @ref UsykRadialLaw into the
   * 0D poroelastic sphere (same functor signature as
   * @ref PoroelasticSpherePassiveLaw).
   */
  template <class Scalar>
  struct UsykPassiveLaw
  {
      RadialPiolaStress<Scalar> operator()(
        const UsykRadialLaw<Scalar>& law, const Scalar lr, const Scalar lt) const
      {
        return law.evaluate(lr, lt);
      }
  };
}

#endif
