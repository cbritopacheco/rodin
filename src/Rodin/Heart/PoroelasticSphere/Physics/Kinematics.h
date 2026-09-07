/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Kinematics.h
 * @brief Exact radial kinematics of the uniformly swollen incompressible
 * spherical wall.
 *
 * With a transmurally uniform Lagrangian porosity @f$ \Phi @f$ the solid
 * incompressibility constraint @f$ J = 1 - \phi_0 + \Phi @f$ integrates
 * exactly,
 * @f[
 *   r^3(R) = r_{in}^3 + J \,(R^3 - R_{in}^3), \qquad
 *   \lambda_\theta = r / R, \qquad \lambda_r = J / \lambda_\theta^2,
 * @f]
 * and variations of @f$ (r_{in}, \Phi) @f$ generate the tangent fields
 * @f$ v_1 = r_{in}^2 / r^2 @f$ (isochoric) and
 * @f$ v_2 = (R^3 - R_{in}^3) / (3 r^2) @f$ (unit swelling).
 */
#ifndef RODIN_HEART_POROELASTICSPHERE_PHYSICS_KINEMATICS_H
#define RODIN_HEART_POROELASTICSPHERE_PHYSICS_KINEMATICS_H

#include <cassert>
#include <cmath>
#include <numbers>
#include <utility>

namespace Rodin::Heart::PoroelasticSphere::Physics
{
  /**
   * @brief Kinematic quantities at one material radius @f$ R @f$, with their
   * total derivatives with respect to the endocardial displacement @f$ y @f$
   * and the porosity @f$ \Phi @f$ (the rates @f$ \dot y, \dot\Phi @f$ are
   * discrete, @f$ \partial \dot q / \partial q = a_0 @f$).
   *
   * @tparam Scalar Scalar numeric type.
   */
  template <class Scalar>
  struct RadialPoint
  {
    Scalar R = 0;  ///< Material radius.
    Scalar r = 0;  ///< Deformed radius.
    Scalar lt = 1; ///< Hoop stretch @f$ \lambda_\theta @f$.
    Scalar lr = 1; ///< Radial stretch @f$ \lambda_r @f$.
    Scalar v1 = 0; ///< Tangent field @f$ v_1 @f$.
    Scalar v2 = 0; ///< Tangent field @f$ v_2 @f$.
    Scalar dv1 = 0; ///< @f$ v_1' = \mathrm d v_1 / \mathrm d R @f$.
    Scalar dv2 = 0; ///< @f$ v_2' = \mathrm d v_2 / \mathrm d R @f$.
    Scalar ltDot = 0; ///< @f$ \dot\lambda_\theta @f$.
    Scalar lrDot = 0; ///< @f$ \dot\lambda_r @f$.

    Scalar lt_y = 0, lt_phi = 0;   ///< Derivatives of @f$ \lambda_\theta @f$.
    Scalar lr_y = 0, lr_phi = 0;   ///< Derivatives of @f$ \lambda_r @f$.
    Scalar v1_y = 0, v1_phi = 0;   ///< Derivatives of @f$ v_1 @f$.
    Scalar v2_y = 0, v2_phi = 0;   ///< Derivatives of @f$ v_2 @f$.
    Scalar dv1_y = 0, dv1_phi = 0; ///< Derivatives of @f$ v_1' @f$.
    Scalar dv2_y = 0, dv2_phi = 0; ///< Derivatives of @f$ v_2' @f$.
    Scalar ltDot_y = 0, ltDot_phi = 0; ///< Total derivatives of @f$ \dot\lambda_\theta @f$.
    Scalar lrDot_y = 0, lrDot_phi = 0; ///< Total derivatives of @f$ \dot\lambda_r @f$.
  };

  /**
   * @brief Evaluates the exact radial kinematics of the reduced wall.
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
       * @brief Evaluate the kinematics at one material radius.
       *
       * @tparam EvalData Evaluation data structure type.
       * @param[in] data Evaluation data; @p rin, @p J, @p yDot, @p phiDot,
       *   @p a0 must already be set.
       * @param[in] R Material radius in @f$ [R_{in}, R_{out}] @f$.
       * @returns Kinematic quantities and derivatives at @p R.
       */
      template <class EvalData>
      RadialPoint<decltype(std::declval<EvalData>().y)> at(
          const EvalData& data, const decltype(std::declval<EvalData>().y) R) const
      {
        using Scalar = decltype(data.y);
        const Scalar Rin = m_input.R0;
        const Scalar rin = data.rin;
        const Scalar J = data.J;

        RadialPoint<Scalar> q;
        q.R = R;
        const Scalar shell = R * R * R - Rin * Rin * Rin;
        q.r = std::cbrt(rin * rin * rin + J * shell);
        assert(q.r > Scalar(0));
        const Scalar r = q.r;
        const Scalar r2 = r * r;

        q.v1 = rin * rin / r2;
        q.v2 = shell / (Scalar(3) * r2);
        q.lt = r / R;
        q.lr = J / (q.lt * q.lt);

        // dr/dy = v1, dr/dphi = v2.
        q.lt_y = q.v1 / R;
        q.lt_phi = q.v2 / R;
        const Scalar lt3 = q.lt * q.lt * q.lt;
        q.lr_y = -Scalar(2) * J / lt3 * q.lt_y;
        q.lr_phi = Scalar(1) / (q.lt * q.lt) - Scalar(2) * J / lt3 * q.lt_phi;

        q.v1_y = Scalar(2) * rin / r2 - Scalar(2) * q.v1 * q.v1 / r;
        q.v1_phi = -Scalar(2) * q.v1 * q.v2 / r;
        q.v2_y = -Scalar(2) * q.v1 * q.v2 / r;
        q.v2_phi = -Scalar(2) * q.v2 * q.v2 / r;

        // v1' = -2 lr v1 / r, v2' = R^2 / r^2 - 2 lr v2 / r.
        q.dv1 = -Scalar(2) * q.lr * q.v1 / r;
        q.dv2 = R * R / r2 - Scalar(2) * q.lr * q.v2 / r;
        q.dv1_y = -Scalar(2) * (q.lr_y * q.v1 / r + q.lr * q.v1_y / r - q.lr * q.v1 * q.v1 / r2);
        q.dv1_phi = -Scalar(2) * (q.lr_phi * q.v1 / r + q.lr * q.v1_phi / r - q.lr * q.v1 * q.v2 / r2);
        q.dv2_y = -Scalar(2) * R * R * q.v1 / (r2 * r)
          - Scalar(2) * (q.lr_y * q.v2 / r + q.lr * q.v2_y / r - q.lr * q.v2 * q.v1 / r2);
        q.dv2_phi = -Scalar(2) * R * R * q.v2 / (r2 * r)
          - Scalar(2) * (q.lr_phi * q.v2 / r + q.lr * q.v2_phi / r - q.lr * q.v2 * q.v2 / r2);

        // Rates: rDot = v1 yDot + v2 phiDot, ltDot = rDot / R,
        // lrDot = phiDot / lt^2 - 2 J ltDot / lt^3.
        const Scalar rDot = q.v1 * data.yDot + q.v2 * data.phiDot;
        const Scalar rDot_y = q.v1_y * data.yDot + q.v1 * data.a0 + q.v2_y * data.phiDot;
        const Scalar rDot_phi = q.v1_phi * data.yDot + q.v2_phi * data.phiDot + q.v2 * data.a0;
        q.ltDot = rDot / R;
        q.ltDot_y = rDot_y / R;
        q.ltDot_phi = rDot_phi / R;

        const Scalar lt2 = q.lt * q.lt;
        const Scalar lt4 = lt2 * lt2;
        q.lrDot = data.phiDot / lt2 - Scalar(2) * J * q.ltDot / lt3;
        q.lrDot_y = -Scalar(2) * data.phiDot * q.lt_y / lt3
          - Scalar(2) * J * q.ltDot_y / lt3 + Scalar(6) * J * q.ltDot * q.lt_y / lt4;
        q.lrDot_phi = data.a0 / lt2 - Scalar(2) * data.phiDot * q.lt_phi / lt3
          - Scalar(2) * q.ltDot / lt3 - Scalar(2) * J * q.ltDot_phi / lt3
          + Scalar(6) * J * q.ltDot * q.lt_phi / lt4;
        return q;
      }

      /**
       * @brief Evaluate the global kinematics and the mid-wall strain used
       * by the (single) contractile unit.
       *
       * Writes @p data.rin, @p data.J, @p data.Vw0, @p data.lambdaMid,
       * @p data.strain1D, @p data.diffGreen, @p data.diffGreenWrtPorosity.
       *
       * @tparam EvalData Evaluation data structure type.
       * @param[in,out] data Evaluation data; candidate unknowns and rates
       *   must already be set.
       */
      template <class EvalData>
      void evaluate(EvalData& data) const
      {
        using Scalar = decltype(data.y);
        const Scalar Rin = m_input.R0;
        const Scalar Rout = m_input.R0 + m_input.d0;

        data.rin = Rin + data.y;
        assert(data.rin > Scalar(0));
        data.J = Scalar(1) - m_input.phi0 + data.phi;
        assert(data.J > Scalar(0));
        data.Vw0 = Scalar(4) / Scalar(3) * std::numbers::pi_v<Scalar>
          * (Rout * Rout * Rout - Rin * Rin * Rin);

        const auto m = at(data, Scalar(0.5) * (Rin + Rout));
        data.lambdaMid = m.lt;
        data.strain1D = Scalar(0.5) * (m.lt * m.lt - Scalar(1));
        data.diffGreen = m.lt * m.lt_y;
        data.diffGreenWrtPorosity = m.lt * m.lt_phi;
      }

    private:
      const Input& m_input;
  };
}

#endif
