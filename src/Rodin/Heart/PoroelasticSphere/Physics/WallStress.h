/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file WallStress.h
 * @brief Skeletal stress, incompressibility multiplier and pressure law of the
 * reduced poroelastic shell.
 *
 * The skeletal nominal stresses collect the passive, viscous and active parts,
 * @f[
 *   P_{M,R} = P^{pas}_{M,R} + \eta \lambda_r^2 \dot\lambda_r, \qquad
 *   P_{M,\theta} = P^{pas}_{M,\theta} + \eta \lambda_\theta^2 \dot\lambda_\theta
 *     + \tfrac12 \sigma_{1D} \lambda_\theta,
 * @f]
 * (viscous @f$ \mathbf S^v = \eta \dot{\mathbf E} @f$; active in-plane average
 * of the CCMLC2014 fiber stress). The radial-stress closure
 * @f$ \bar\sigma_{rr} = -p_v/2 @f$ fixes the multiplier and the shell
 * equilibrium gives the pressure law:
 * @f[
 *   \bar\lambda = -\frac{p_v}{2} - \frac{P_{M,R}}{\lambda_\theta^2}, \qquad
 *   \bar P_{\theta\theta} = P_{M,\theta} + \bar\lambda \lambda_r \lambda_\theta, \qquad
 *   p_v = \frac{2 h_0}{R_0 \lambda_\theta^2} \bar P_{\theta\theta}.
 * @f]
 */
#ifndef RODIN_HEART_POROELASTICSPHERE_PHYSICS_WALLSTRESS_H
#define RODIN_HEART_POROELASTICSPHERE_PHYSICS_WALLSTRESS_H

#include <cmath>

namespace Rodin::Heart::PoroelasticSphere::Physics
{
  /**
   * @brief Evaluates the skeletal stresses, @f$ \bar\lambda @f$ and the shell
   * pressure law with their derivatives with respect to
   * @f$ (y, \bar\Phi, p_v, e_c) @f$.
   *
   * @tparam PassiveLaw Functor type returning the radial Piola components.
   * @tparam Input Model input parameter type.
   */
  template <class PassiveLaw, class Input>
  class WallStressEvaluator
  {
    public:
      /**
       * @brief Construct with model input parameters.
       * @param[in] input Model parameters.
       */
      explicit WallStressEvaluator(const Input& input)
        : m_input(input)
      {}

      /**
       * @brief Evaluate stresses, multiplier and wall pressure.
       *
       * Requires the kinematics and the active stress in @p data.
       *
       * @tparam EvalData Evaluation data structure type.
       * @param[in,out] data Evaluation data.
       */
      template <class EvalData>
      void evaluate(EvalData& data) const
      {
        using Scalar = decltype(data.y);
        const Scalar eta = m_input.eta;
        const Scalar lt = data.lt;
        const Scalar lr = data.lr;
        const Scalar lt2 = lt * lt;
        const Scalar lt3 = lt2 * lt;

        const Scalar sigma = data.active.activeStressOneDimensional;
        const Scalar sigma_y = data.active.partialActiveStressWrtDisplacement;
        const Scalar sigma_ec = data.active.partialActiveStressWrtFiberDeformation;

        PassiveLaw passiveLaw;
        const auto pas = passiveLaw(m_input.passiveEnergy, lr, lt);

        // Radial skeletal stress (passive + viscous).
        data.PMR = pas.PRR + eta * lr * lr * data.lrDot;
        data.PMR_y = pas.dPRR_dlr * data.lr_y + pas.dPRR_dlt * data.lt_y +
          eta * (Scalar(2) * lr * data.lr_y * data.lrDot + lr * lr * data.lrDot_y);
        data.PMR_phi = pas.dPRR_dlr * data.lr_phi +
          eta * (Scalar(2) * lr * data.lr_phi * data.lrDot + lr * lr * data.lrDot_phi);

        // Tangential skeletal stress (passive + viscous + active).
        data.PMt = pas.Ptt + eta * lt2 * data.ltDot + Scalar(0.5) * sigma * lt;
        data.PMt_y = pas.dPtt_dlr * data.lr_y + pas.dPtt_dlt * data.lt_y +
          eta * (Scalar(2) * lt * data.lt_y * data.ltDot + lt2 * data.ltDot_y) +
          Scalar(0.5) * (sigma_y * lt + sigma * data.lt_y);
        data.PMt_phi = pas.dPtt_dlr * data.lr_phi;
        data.PMt_ec = Scalar(0.5) * sigma_ec * lt;

        // Multiplier from the radial-stress closure sigma_rr = -p_v / 2.
        data.lambdaBar = -Scalar(0.5) * data.pv - data.PMR / lt2;
        data.dLambdaBar_dy =
          -(data.PMR_y / lt2 - Scalar(2) * data.PMR * data.lt_y / lt3);
        data.dLambdaBar_dphi = -data.PMR_phi / lt2;
        data.dLambdaBar_dpv = -Scalar(0.5);
        data.dLambdaBar_dec = Scalar(0);

        // Total tangential nominal stress and shell pressure law.
        const Scalar Ptt = data.PMt + data.lambdaBar * lr * lt;
        const Scalar Ptt_y = data.PMt_y + data.dLambdaBar_dy * lr * lt +
          data.lambdaBar * (data.lr_y * lt + lr * data.lt_y);
        const Scalar Ptt_phi = data.PMt_phi + data.dLambdaBar_dphi * lr * lt +
          data.lambdaBar * data.lr_phi * lt;
        const Scalar Ptt_pv = data.dLambdaBar_dpv * lr * lt;
        const Scalar Ptt_ec = data.PMt_ec;
        data.Ptt = Ptt;

        const Scalar c = Scalar(2) * m_input.d0 / m_input.R0;
        data.wallPressure = c * Ptt / lt2;
        data.dWallPressure_dy = c * (Ptt_y / lt2 - Scalar(2) * Ptt * data.lt_y / lt3);
        data.dWallPressure_dphi = c * Ptt_phi / lt2;
        data.dWallPressure_dpv = c * Ptt_pv / lt2;
        data.dWallPressure_dec = c * Ptt_ec / lt2;
      }

    private:
      const Input& m_input;
  };
}

#endif
