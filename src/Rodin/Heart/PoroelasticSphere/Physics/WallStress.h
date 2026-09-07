/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file WallStress.h
 * @brief Transmural quadrature of the two tangent projections of the
 * momentum balance of the poroelastic sphere.
 *
 * With the raw first Piola-Kirchhoff components
 * @f$ P_{RR}, P_{\theta\theta} @f$ (passive + viscous + active) evaluated
 * along the exact kinematics, the projection on @f$ v_1 @f$ gives the
 * pressure-volume law
 * @f[
 *   P_v\, r_{in}^2 = \int_{R_{in}}^{R_{out}}
 *     \Big[ P_{RR} v_1' + 2 P_{\theta\theta} \frac{v_1}{R} \Big] R^2 \,\mathrm dR,
 * @f]
 * in which the incompressibility multiplier cancels identically, and the
 * projection on @f$ v_2 @f$ determines its wall average
 * @f[
 *   \bar\lambda = -\frac{3}{R_{out}^3 - R_{in}^3} \int_{R_{in}}^{R_{out}}
 *     \Big[ P_{RR} v_2' + 2 P_{\theta\theta} \frac{v_2}{R} \Big] R^2 \,\mathrm dR .
 * @f]
 * The active stress is the in-plane average of the CCMLC2014 fiber stress,
 * @f$ \mathbf\Sigma^{act} = \tfrac12 \sigma_{1D} (\mathbf I - \mathbf e_R \otimes \mathbf e_R) @f$,
 * i.e. @f$ P^{act}_{\theta\theta} = \tfrac12 \sigma_{1D} \lambda_\theta @f$,
 * and the viscous stress is @f$ \mathbf\Sigma^{v} = \eta \dot{\mathbf e} @f$,
 * i.e. @f$ P^{v}_{RR} = \eta \lambda_r^2 \dot\lambda_r @f$,
 * @f$ P^{v}_{\theta\theta} = \eta \lambda_\theta^2 \dot\lambda_\theta @f$;
 * both reduce to the CCMLC2014 terms in the membrane limit.
 */
#ifndef RODIN_HEART_POROELASTICSPHERE_PHYSICS_WALLSTRESS_H
#define RODIN_HEART_POROELASTICSPHERE_PHYSICS_WALLSTRESS_H

#include <cassert>
#include <cmath>
#include <numbers>
#include <vector>

#include "Rodin/Types.h"
#include "Rodin/Heart/PoroelasticSphere/Physics/Kinematics.h"

namespace Rodin::Heart::PoroelasticSphere::Physics
{
  /**
   * @brief Evaluates the wall pressure-volume law and the averaged
   * multiplier by Gauss-Legendre quadrature across the wall.
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
       * @param[in] input Model parameters (passive energy law, geometry).
       */
      explicit WallStressEvaluator(const Input& input)
        : m_input(input), m_kinematics(input)
      {
        gaussLegendreUnit(input.wallQuadraturePoints, m_x, m_w);
      }

      /**
       * @brief Evaluate the wall pressure and averaged multiplier with their
       * derivatives.
       *
       * Reads the kinematic fields and @p data.active, writes
       * @p data.wallPressure, @p data.lambdaBar and their derivatives with
       * respect to @f$ y, \Phi, e_c @f$.
       *
       * @tparam EvalData Evaluation data structure type.
       * @param[in,out] data Evaluation data with kinematics and active stress
       *   already computed.
       */
      template <class EvalData>
      void evaluate(EvalData& data) const
      {
        using Scalar = decltype(data.y);
        const Scalar Rin = m_input.R0;
        const Scalar Rout = m_input.R0 + m_input.d0;
        const Scalar eta = m_input.eta;
        const Scalar sigma = data.active.activeStressOneDimensional;
        const Scalar sigma_y = data.active.partialActiveStressWrtDisplacement;
        const Scalar sigma_phi = data.active.partialActiveStressWrtPorosity;
        const Scalar sigma_ec = data.active.partialActiveStressWrtFiberDeformation;

        Scalar I1 = 0, I1_y = 0, I1_phi = 0, I1_ec = 0;
        Scalar I2 = 0, I2_y = 0, I2_phi = 0, I2_ec = 0;

        PassiveLaw passiveLaw;
        for (size_t k = 0; k < m_x.size(); ++k)
        {
          const Scalar R = Rin + (Rout - Rin) * m_x[k];
          const Scalar wq = (Rout - Rin) * m_w[k] * R * R;
          const auto q = m_kinematics.at(data, R);
          const auto pas = passiveLaw(m_input.passiveEnergy, q.lr, q.lt);

          // P_RR = passive + eta lr^2 lrDot.
          const Scalar PRR = pas.PRR + eta * q.lr * q.lr * q.lrDot;
          const Scalar PRR_y = pas.dPRR_dlr * q.lr_y + pas.dPRR_dlt * q.lt_y
            + eta * (Scalar(2) * q.lr * q.lr_y * q.lrDot + q.lr * q.lr * q.lrDot_y);
          const Scalar PRR_phi = pas.dPRR_dlr * q.lr_phi + pas.dPRR_dlt * q.lt_phi
            + eta * (Scalar(2) * q.lr * q.lr_phi * q.lrDot + q.lr * q.lr * q.lrDot_phi);

          // P_tt = passive + eta lt^2 ltDot + sigma lt / 2.
          const Scalar Ptt = pas.Ptt + eta * q.lt * q.lt * q.ltDot + Scalar(0.5) * sigma * q.lt;
          const Scalar Ptt_y = pas.dPtt_dlr * q.lr_y + pas.dPtt_dlt * q.lt_y
            + eta * (Scalar(2) * q.lt * q.lt_y * q.ltDot + q.lt * q.lt * q.ltDot_y)
            + Scalar(0.5) * (sigma_y * q.lt + sigma * q.lt_y);
          const Scalar Ptt_phi = pas.dPtt_dlr * q.lr_phi + pas.dPtt_dlt * q.lt_phi
            + eta * (Scalar(2) * q.lt * q.lt_phi * q.ltDot + q.lt * q.lt * q.ltDot_phi)
            + Scalar(0.5) * (sigma_phi * q.lt + sigma * q.lt_phi);
          const Scalar Ptt_ec = Scalar(0.5) * sigma_ec * q.lt;

          I1 += wq * (PRR * q.dv1 + Scalar(2) * Ptt * q.v1 / R);
          I1_y += wq * (PRR_y * q.dv1 + PRR * q.dv1_y
            + Scalar(2) * (Ptt_y * q.v1 + Ptt * q.v1_y) / R);
          I1_phi += wq * (PRR_phi * q.dv1 + PRR * q.dv1_phi
            + Scalar(2) * (Ptt_phi * q.v1 + Ptt * q.v1_phi) / R);
          I1_ec += wq * Scalar(2) * Ptt_ec * q.v1 / R;

          I2 += wq * (PRR * q.dv2 + Scalar(2) * Ptt * q.v2 / R);
          I2_y += wq * (PRR_y * q.dv2 + PRR * q.dv2_y
            + Scalar(2) * (Ptt_y * q.v2 + Ptt * q.v2_y) / R);
          I2_phi += wq * (PRR_phi * q.dv2 + PRR * q.dv2_phi
            + Scalar(2) * (Ptt_phi * q.v2 + Ptt * q.v2_phi) / R);
          I2_ec += wq * Scalar(2) * Ptt_ec * q.v2 / R;
        }

        const Scalar rin = data.rin;
        const Scalar rin2 = rin * rin;
        data.wallPressure = I1 / rin2;
        data.dWallPressure_dy = I1_y / rin2 - Scalar(2) * I1 / (rin2 * rin);
        data.dWallPressure_dphi = I1_phi / rin2;
        data.dWallPressure_dec = I1_ec / rin2;

        const Scalar c = -Scalar(3) / (Rout * Rout * Rout - Rin * Rin * Rin);
        data.lambdaBar = c * I2;
        data.dLambdaBar_dy = c * I2_y;
        data.dLambdaBar_dphi = c * I2_phi;
        data.dLambdaBar_dec = c * I2_ec;
      }

    private:
      /**
       * @brief @f$ n @f$-point Gauss-Legendre nodes and weights on @f$ [0, 1] @f$
       * (zeros of @f$ P_n @f$ by Newton iteration on the three-term recurrence).
       */
      static void gaussLegendreUnit(
          size_t n, std::vector<Real>& x, std::vector<Real>& w)
      {
        assert(n >= 1);
        x.assign(n, Real(0));
        w.assign(n, Real(0));
        const size_t m = (n + 1) / 2;
        for (size_t k = 0; k < m; ++k)
        {
          Real z = std::cos(std::numbers::pi_v<Real> * (Real(k) + Real(0.75))
              / (Real(n) + Real(0.5)));
          Real dp = Real(1);
          for (size_t it = 0; it < 100; ++it)
          {
            Real p0 = Real(1);
            Real p1 = z;
            for (size_t j = 2; j <= n; ++j)
            {
              const Real p = ((Real(2) * j - Real(1)) * z * p1 - (Real(j) - Real(1)) * p0) / Real(j);
              p0 = p1;
              p1 = p;
            }
            dp = Real(n) * (z * p1 - p0) / (z * z - Real(1));
            const Real dz = p1 / dp;
            z -= dz;
            if (std::abs(dz) < Real(1e-15))
              break;
          }
          // Map [-1, 1] -> [0, 1]; symmetric node n - 1 - k.
          const Real wk = Real(1) / ((Real(1) - z * z) * dp * dp);
          x[k] = Real(0.5) * (Real(1) - z);
          x[n - 1 - k] = Real(0.5) * (Real(1) + z);
          w[k] = wk;
          w[n - 1 - k] = wk;
        }
      }

      const Input& m_input;
      WallKinematics<Input> m_kinematics;
      std::vector<Real> m_x;
      std::vector<Real> m_w;
  };
}

#endif
