/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file PassiveLaw.h
 * @brief Radial first Piola-Kirchhoff components of the passive skeleton
 * of the spherical wall.
 *
 * Two passive responses are provided: the reduced exponential energy of Usyk
 * type restricted to the spherical stretches (default), and the CCMLC2014
 * reduced-invariant energy evaluated on the same kinematics.
 */
#ifndef RODIN_HEART_POROELASTICSPHERE_PASSIVELAW_H
#define RODIN_HEART_POROELASTICSPHERE_PASSIVELAW_H

#include <cmath>

#include "Rodin/Heart/CCMLC2014/PassiveEnergy.h"

namespace Rodin::Heart
{
  /**
   * @brief Radial first Piola-Kirchhoff components and their tangent.
   *
   * For @f$ \mathbf F = \mathrm{diag}(\lambda_r, \lambda_\theta, \lambda_\theta) @f$
   * and @f$ \widehat W(\lambda_r, \lambda_\theta) = \Psi_M(\mathbf F) @f$:
   * @f$ P_{RR} = \partial_{\lambda_r} \widehat W @f$,
   * @f$ P_{\theta\theta} = \tfrac12 \partial_{\lambda_\theta} \widehat W @f$.
   *
   * @tparam Scalar Scalar numeric type.
   */
  template <class Scalar>
  struct RadialPiolaStress
  {
      Scalar PRR = Scalar(0);   ///< @f$ P_{RR} @f$.
      Scalar Ptt = Scalar(0);   ///< @f$ P_{\theta\theta} = P_{\varphi\varphi} @f$.
      Scalar dPRR_dlr = Scalar(0); ///< @f$ \partial P_{RR} / \partial \lambda_r @f$.
      Scalar dPRR_dlt = Scalar(0); ///< @f$ \partial P_{RR} / \partial \lambda_\theta @f$.
      Scalar dPtt_dlr =
        Scalar(0); ///< @f$ \partial P_{\theta\theta} / \partial \lambda_r @f$.
      Scalar dPtt_dlt =
        Scalar(0); ///< @f$ \partial P_{\theta\theta} / \partial \lambda_\theta @f$.
  };

  /**
   * @brief Reduced exponential passive energy of Usyk type restricted to
   * @f$ \mathbf F = \mathrm{diag}(\lambda_r, \lambda_\theta, \lambda_\theta) @f$,
   * @f[
   *   \widehat\Psi_M = \frac{c_0}{2}\big(e^Q - 1\big), \quad
   *   Q = \tfrac14\big[b_{rr}(\lambda_r^2 - 1)^2 + 2 b_{\theta\theta}(\lambda_\theta^2 - 1)^2
   *     + 2 b_{r\theta}(\lambda_r^2 - 1)(\lambda_\theta^2 - 1)\big],
   * @f]
   * with @f$ P_{M,R} = \partial_{\lambda_r}\widehat\Psi_M @f$ and
   * @f$ P_{M,\theta} = \tfrac12 \partial_{\lambda_\theta}\widehat\Psi_M @f$.
   *
   * @tparam Scalar Scalar numeric type.
   */
  template <class Scalar>
  class UsykSphericalEnergy
  {
    public:
      /// @brief Material parameters.
      struct Parameters
      {
        Scalar c0 = Scalar(100);    ///< Stress scale @f$ c_0 @f$ (Pa).
        Scalar brr = Scalar(1.24);  ///< Radial exponent @f$ b_{rr} @f$.
        Scalar btt = Scalar(8.4);   ///< Tangential exponent @f$ b_{\theta\theta} @f$.
        Scalar brt = Scalar(0);     ///< Coupling exponent @f$ b_{r\theta} @f$.
      };

      /// @brief Constructs the law with a parameter set.
      explicit UsykSphericalEnergy(const Parameters& params = Parameters())
        : m_params(params)
      {}

      /// @brief Returns the parameter set.
      const Parameters& getParameters() const noexcept
      {
        return m_params;
      }

      /**
       * @brief Radial Piola components and their tangent.
       * @param[in] lr Radial stretch @f$ \lambda_r @f$.
       * @param[in] lt Hoop stretch @f$ \lambda_\theta @f$.
       */
      RadialPiolaStress<Scalar> radialStress(const Scalar lr, const Scalar lt) const
      {
        const auto& p = m_params;
        const Scalar er = lr * lr - Scalar(1);
        const Scalar et = lt * lt - Scalar(1);
        const Scalar Q = Scalar(0.25) *
          (p.brr * er * er + Scalar(2) * p.btt * et * et + Scalar(2) * p.brt * er * et);
        const Scalar e = std::exp(Q);

        const Scalar Qr = lr * (p.brr * er + p.brt * et);
        const Scalar Qt = lt * (Scalar(2) * p.btt * et + p.brt * er);
        const Scalar Qrr = p.brr * er + p.brt * et + Scalar(2) * p.brr * lr * lr;
        const Scalar Qrt = Scalar(2) * p.brt * lr * lt;
        const Scalar Qtt =
          Scalar(2) * p.btt * et + p.brt * er + Scalar(4) * p.btt * lt * lt;

        const Scalar aR = Scalar(0.5) * p.c0 * e;
        const Scalar aT = Scalar(0.25) * p.c0 * e;

        RadialPiolaStress<Scalar> out;
        out.PRR = aR * Qr;
        out.Ptt = aT * Qt;
        out.dPRR_dlr = aR * (Qr * Qr + Qrr);
        out.dPRR_dlt = aR * (Qt * Qr + Qrt);
        out.dPtt_dlr = aT * (Qr * Qt + Qrt);
        out.dPtt_dlt = aT * (Qt * Qt + Qtt);
        return out;
      }

    private:
      Parameters m_params;
  };

  /**
   * @brief Default passive stress law of the poroelastic sphere: forwards to
   * UsykSphericalEnergy::radialStress.
   *
   * @tparam Scalar Scalar numeric type.
   */
  template <class Scalar>
  struct UsykSpherePassiveLaw
  {
      /**
       * @brief Evaluate the radial Piola components and their tangent.
       * @param[in] law Usyk energy.
       * @param[in] lr Radial stretch.
       * @param[in] lt Hoop stretch.
       */
      template <class PassiveEnergyLaw>
      RadialPiolaStress<Scalar> operator()(
        const PassiveEnergyLaw& law, const Scalar lr, const Scalar lt) const
      {
        return law.radialStress(lr, lt);
      }
  };

  /**
   * @brief CCMLC2014 reduced-invariant passive stress on the radial kinematics.
   *
   * Evaluates the CCMLC2014 reduced-invariant energy @f$ W(J_1, J_2, J_4) @f$
   * on the radial deformation @f$ \mathbf F = \mathrm{diag}(\lambda_r,
   * \lambda_\theta, \lambda_\theta) @f$, @f$ J = \lambda_r \lambda_\theta^2 @f$,
   * with the isochoric invariants
   * @f$ J_1 = I_1 J^{-2/3} = \lambda_r^{4/3}\lambda_\theta^{-4/3} + 2\lambda_r^{-2/3}\lambda_\theta^{2/3} @f$,
   * @f$ J_2 = I_2 J^{-4/3} = 2\lambda_r^{2/3}\lambda_\theta^{-2/3} + \lambda_r^{-4/3}\lambda_\theta^{4/3} @f$
   * and the in-plane fiber invariant @f$ J_4 = \lambda_\theta^2 @f$. At
   * @f$ J = 1 @f$ these are the CCMLC2014 reduced invariants
   * @f$ J_1 = 2C + C^{-2} @f$, @f$ J_2 = C^2 + 2C^{-1} @f$, @f$ J_4 = C @f$
   * with @f$ C = \lambda_\theta^2 @f$; the volumetric response is carried by
   * the incompressibility multiplier, not by the law.
   *
   * @tparam Scalar Scalar numeric type.
   */
  template <class Scalar>
  struct PoroelasticSpherePassiveLaw
  {
      /**
     * @brief Evaluate the radial Piola components and their tangent.
     *
     * @tparam PassiveEnergyLaw Passive-energy evaluator type.
     * @param[in] law Passive-energy law.
     * @param[in] lr Radial stretch @f$ \lambda_r @f$.
     * @param[in] lt Hoop stretch @f$ \lambda_\theta @f$.
     * @returns Piola components and their derivatives.
     */
      template <class PassiveEnergyLaw>
      RadialPiolaStress<Scalar> operator()(
        const PassiveEnergyLaw& law, const Scalar lr, const Scalar lt) const
      {
        const Scalar third = Scalar(1) / Scalar(3);
        auto p = [](Scalar x, Scalar e) { return std::pow(x, e); };

        // Invariants J1 = a^{4/3} b^{-4/3} + 2 a^{-2/3} b^{2/3},
        // J2 = 2 a^{2/3} b^{-2/3} + a^{-4/3} b^{4/3}, J4 = b^2, (a, b) = (lr, lt).
        const Scalar J1 = p(lr, 4 * third) * p(lt, -4 * third) +
          Scalar(2) * p(lr, -2 * third) * p(lt, 2 * third);
        const Scalar J2 = Scalar(2) * p(lr, 2 * third) * p(lt, -2 * third) +
          p(lr, -4 * third) * p(lt, 4 * third);
        const Scalar J4 = lt * lt;

        // First derivatives.
        const Scalar J1_r = (Scalar(4) * third) *
          (p(lr, third) * p(lt, -4 * third) - p(lr, -5 * third) * p(lt, 2 * third));
        const Scalar J1_t = (Scalar(4) * third) *
          (-p(lr, 4 * third) * p(lt, -7 * third) + p(lr, -2 * third) * p(lt, -third));
        const Scalar J2_r = (Scalar(4) * third) *
          (p(lr, -third) * p(lt, -2 * third) - p(lr, -7 * third) * p(lt, 4 * third));
        const Scalar J2_t = (Scalar(4) * third) *
          (-p(lr, 2 * third) * p(lt, -5 * third) + p(lr, -4 * third) * p(lt, third));
        const Scalar J4_t = Scalar(2) * lt;

        // Second derivatives.
        const Scalar ninth = third * third;
        const Scalar J1_rr = ninth *
          (Scalar(4) * p(lr, -2 * third) * p(lt, -4 * third) +
            Scalar(20) * p(lr, -8 * third) * p(lt, 2 * third));
        const Scalar J1_rt = -ninth *
          (Scalar(16) * p(lr, third) * p(lt, -7 * third) +
            Scalar(8) * p(lr, -5 * third) * p(lt, -third));
        const Scalar J1_tt = ninth *
          (Scalar(28) * p(lr, 4 * third) * p(lt, -10 * third) -
            Scalar(4) * p(lr, -2 * third) * p(lt, -4 * third));
        const Scalar J2_rr = ninth *
          (-Scalar(4) * p(lr, -4 * third) * p(lt, -2 * third) +
            Scalar(28) * p(lr, -10 * third) * p(lt, 4 * third));
        const Scalar J2_rt = -ninth *
          (Scalar(8) * p(lr, -third) * p(lt, -5 * third) +
            Scalar(16) * p(lr, -7 * third) * p(lt, third));
        const Scalar J2_tt = ninth *
          (Scalar(20) * p(lr, 2 * third) * p(lt, -8 * third) +
            Scalar(4) * p(lr, -4 * third) * p(lt, -2 * third));
        const Scalar J4_tt = Scalar(2);

        const ReducedInvariants<Scalar> I{J1, J2, J4};
        const auto D = law.evaluate(I);

        const Scalar W1 = D.grad.dW_dJ1;
        const Scalar W2 = D.grad.dW_dJ2;
        const Scalar W4 = D.grad.dW_dJ4;
        const Scalar W11 = D.hess.d2W_dJ1dJ1;
        const Scalar W12 = D.hess.d2W_dJ1dJ2;
        const Scalar W14 = D.hess.d2W_dJ1dJ4;
        const Scalar W22 = D.hess.d2W_dJ2dJ2;
        const Scalar W24 = D.hess.d2W_dJ2dJ4;
        const Scalar W44 = D.hess.d2W_dJ4dJ4;

        // dW_a / dlr and dW_a / dlt by the chain rule through (J1, J2, J4).
        const Scalar W1_r = W11 * J1_r + W12 * J2_r;
        const Scalar W2_r = W12 * J1_r + W22 * J2_r;
        const Scalar W4_r = W14 * J1_r + W24 * J2_r;
        const Scalar W1_t = W11 * J1_t + W12 * J2_t + W14 * J4_t;
        const Scalar W2_t = W12 * J1_t + W22 * J2_t + W24 * J4_t;
        const Scalar W4_t = W14 * J1_t + W24 * J2_t + W44 * J4_t;

        RadialPiolaStress<Scalar> out;
        out.PRR = W1 * J1_r + W2 * J2_r;
        out.Ptt = Scalar(0.5) * (W1 * J1_t + W2 * J2_t + W4 * J4_t);

        out.dPRR_dlr = W1_r * J1_r + W1 * J1_rr + W2_r * J2_r + W2 * J2_rr;
        out.dPRR_dlt = W1_t * J1_r + W1 * J1_rt + W2_t * J2_r + W2 * J2_rt;
        out.dPtt_dlr = Scalar(0.5) *
          (W1_r * J1_t + W1 * J1_rt + W2_r * J2_t + W2 * J2_rt + W4_r * J4_t);
        out.dPtt_dlt = Scalar(0.5) *
          (W1_t * J1_t + W1 * J1_tt + W2_t * J2_t + W2 * J2_tt + W4_t * J4_t +
            W4 * J4_tt);
        return out;
      }
  };
}

#endif
