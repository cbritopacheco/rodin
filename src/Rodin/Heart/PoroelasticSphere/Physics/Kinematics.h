/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Kinematics.h
 * @brief Representative shell kinematics of the reduced poroelastic sphere.
 *
 * The wall is a spherical shell of mid-surface radius @f$ R_0 @f$ and
 * reference thickness @f$ h_0 @f$, described by the mid-surface displacement
 * @f$ y @f$ and the volume-averaged Lagrangian porosity @f$ \bar\Phi @f$:
 * @f[
 *   \bar\lambda_\theta = 1 + \frac{y}{R_0}, \qquad
 *   \bar\lambda_r = \frac{1 - \phi_0 + \bar\Phi}{\bar\lambda_\theta^2}, \qquad
 *   h = h_0 \bar\lambda_r,
 * @f]
 * the second relation being the solid-volume constraint
 * @f$ \bar\lambda_r \bar\lambda_\theta^2 = 1 - \phi_0 + \bar\Phi @f$.
 */
#ifndef RODIN_HEART_POROELASTICSPHERE_PHYSICS_KINEMATICS_H
#define RODIN_HEART_POROELASTICSPHERE_PHYSICS_KINEMATICS_H

#include <cassert>
#include <cmath>
#include <numbers>

namespace Rodin::Heart::PoroelasticSphere::Physics
{
  /**
   * @brief Evaluates the representative stretches, their rates and their
   * total derivatives with respect to @f$ (y, \bar\Phi) @f$.
   *
   * Rates are discrete (BE/BDF2), so @f$ \partial \dot q / \partial q = a_0 @f$.
   *
   * @tparam Input Model input parameter type.
   */
  template <class Input>
  class WallKinematics
  {
    public:
      /**
       * @brief Construct with a reference to model input parameters.
       * @param[in] input Model parameters (geometry, reference porosity).
       */
      explicit WallKinematics(const Input& input)
        : m_input(input)
      {}

      /**
       * @brief Evaluate the shell kinematics.
       *
       * Reads @p data.y, @p data.phi, @p data.yDot, @p data.phiDot, @p data.a0
       * and writes the stretches, rates, cavity radius, wall volume and the
       * hoop Green-Lagrange strain of the contractile unit.
       *
       * @tparam EvalData Evaluation data structure type.
       * @param[in,out] data Evaluation data.
       */
      template <class EvalData>
      void evaluate(EvalData& data) const
      {
        using Scalar = decltype(data.y);
        const Scalar R0 = m_input.R0;
        const Scalar h0 = m_input.d0;
        const Scalar pi = std::numbers::pi_v<Scalar>;

        data.radius = R0 + data.y;
        assert(data.radius > Scalar(0));
        data.J = Scalar(1) - m_input.phi0 + data.phi;
        assert(data.J > Scalar(0));

        const Scalar Rin = R0 - Scalar(0.5) * h0;
        const Scalar Rout = R0 + Scalar(0.5) * h0;
        data.Vw0 = Scalar(4) / Scalar(3) * pi * (Rout * Rout * Rout - Rin * Rin * Rin);
        data.Aw0 = Scalar(4) * pi * R0 * R0;

        const Scalar lt = Scalar(1) + data.y / R0;
        const Scalar lt2 = lt * lt;
        const Scalar lt3 = lt2 * lt;
        const Scalar lt4 = lt2 * lt2;
        const Scalar J = data.J;

        data.lt = lt;
        data.lt_y = Scalar(1) / R0;
        data.lr = J / lt2;
        data.lr_y = -Scalar(2) * J / lt3 * data.lt_y;
        data.lr_phi = Scalar(1) / lt2;
        data.thickness = h0 * data.lr;

        data.ltDot = data.yDot / R0;
        data.ltDot_y = data.a0 / R0;
        data.lrDot = data.phiDot / lt2 - Scalar(2) * J * data.ltDot / lt3;
        data.lrDot_y = -Scalar(2) * data.phiDot * data.lt_y / lt3 -
          Scalar(2) * J * data.ltDot_y / lt3 +
          Scalar(6) * J * data.ltDot * data.lt_y / lt4;
        data.lrDot_phi = data.a0 / lt2 - Scalar(2) * data.ltDot / lt3;

        // Hoop Green-Lagrange strain of the single contractile unit.
        data.lambdaMid = lt;
        data.strain1D = Scalar(0.5) * (lt2 - Scalar(1));
        data.diffGreen = lt * data.lt_y;
      }

    private:
      const Input& m_input;
  };
}

#endif
