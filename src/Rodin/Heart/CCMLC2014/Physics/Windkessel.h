/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Windkessel.h
 * @brief Windkessel outflow evaluator for the arterial branch.
 *
 * Implements the three-element Windkessel outflow term
 * @f$ Q_{ar} = K_{ar} (p_v - p_{ar}) @f$ when the aortic valve is open
 * (Caruel et al. 2014, §4, arterial outflow model).
 */
#ifndef RODIN_HEART_CCMLC2014_PHYSICS_WINDKESSEL_H
#define RODIN_HEART_CCMLC2014_PHYSICS_WINDKESSEL_H

namespace Rodin::Heart::CCMLC2014::Physics
{
  /**
   * @brief Evaluates the Windkessel outflow and its pressure derivatives.
   *
   * The outflow is the positive part of the aortic flow rate
   * @f$ Q_{ar}^+ = \max(0, K_{ar}(p_{v,\text{mid}} - p_{ar,\text{mid}})) @f$.
   *
   * @tparam Input Model input parameter type.
   */
  template <class Input>
  class WindkesselOutflowEvaluator
  {
    public:
      /**
       * @brief Construct with model input parameters.
       * @param[in] input Model parameters (aortic valve conductance).
       */
      explicit WindkesselOutflowEvaluator(const Input& input)
        : m_input(input)
      {}

      /**
       * @brief Evaluate outflow and pressure derivatives.
       *
       * Reads midpoint pressures from @p data and writes
       * @p data.windkesselOutflow and its pressure derivatives.
       *
       * @tparam EvalData Evaluation data structure type.
       * @param[in,out] data Evaluation data with midpoint pressures set.
       */
      template <class EvalData>
      void evaluate(EvalData& data) const
      {
        using Scalar = decltype(data.y);
        const Scalar aorticFlowRate =
          m_input.Kar * (data.pvMid - data.parMid);

        data.windkesselOutflow =
          (aorticFlowRate > Scalar(0)) ? aorticFlowRate : Scalar(0);

        const Scalar midpointDerivative =
          (aorticFlowRate > Scalar(0))
            ? (Scalar(0.5) * m_input.Kar)
            : Scalar(0);

        data.dWindkesselOutflow_dPv = midpointDerivative;
        data.dWindkesselOutflow_dPar = -midpointDerivative;
      }

    private:
      const Input& m_input;
  };
}

#endif
