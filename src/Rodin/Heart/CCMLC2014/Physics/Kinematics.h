/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Kinematics.h
 * @brief Midpoint kinematics evaluator for the CCMLC2014 reduced 0D cavity model.
 *
 * Computes midpoint time-discretized quantities and reduced kinematic scalars
 * (right Cauchy-Green measure, Green-Lagrange strain) at the midpoint
 * configuration, following Caruel et al. (2014), §4.
 */
#ifndef RODIN_HEART_CCMLC2014_PHYSICS_KINEMATICS_H
#define RODIN_HEART_CCMLC2014_PHYSICS_KINEMATICS_H

#include <cassert>

namespace Rodin::Heart::CCMLC2014::Physics
{
  /**
   * @brief Evaluates midpoint kinematics for the 0D reduced cavity model.
   *
   * Given the candidate unknowns at @f$ t_{n+1} @f$ and the previous values
   * at @f$ t_n @f$, computes:
   * - Midpoint displacement @f$ y_{\text{mid}} @f$
   * - Midpoint pressures @f$ p_{v,\text{mid}},\; p_{ar,\text{mid}},\; p_{d,\text{mid}} @f$
   * - Midpoint radial velocity @f$ \dot{y} \approx (y^{n+1} - y^n)/\Delta t @f$
   * - Reduced right Cauchy-Green scalar @f$ C = (1 + y_{\text{mid}} / R_0)^2 @f$
   * - Green-Lagrange strain @f$ e = \tfrac{1}{2}(C - 1) @f$
   * - Strain derivative @f$ \partial e / \partial y = \sqrt{C} / R_0 @f$
   *
   * @tparam Input Model input parameter type.
   */
  template <class Input>
  class MidpointKinematics
  {
    public:
      /**
       * @brief Construct with a reference to model input parameters.
       * @param[in] input Model parameters (geometry, material constants).
       */
      explicit MidpointKinematics(const Input& input)
        : m_input(input)
      {}

      /**
       * @brief Evaluate midpoint kinematics and populate derived fields in @p data.
       *
       * @tparam EvalData Evaluation data structure type.
       * @param[in,out] data Evaluation data; candidate unknowns and previous
       *   state must already be set.
       */
      template <class EvalData>
      void evaluate(EvalData& data) const
      {
        using Scalar = decltype(data.y);
        data.yMid = Scalar(0.5) * (data.y + data.yPrev);
        data.pvMid = Scalar(0.5) * (data.pv + data.pvPrev);
        data.parMid = Scalar(0.5) * (data.par + data.parPrev);
        data.pdMid = Scalar(0.5) * (data.pd + data.pdPrev);
        data.vel = (data.y - data.yPrev) / data.dt;

        data.sqrtC = Scalar(1) + data.yMid / m_input.R0;
        assert(data.sqrtC > Scalar(0));
        data.C = data.sqrtC * data.sqrtC;

        data.strain1D = Scalar(0.5) * (data.C - Scalar(1));
        data.diffGreen = data.sqrtC / m_input.R0;
      }

    private:
      const Input& m_input;
  };
}

#endif
