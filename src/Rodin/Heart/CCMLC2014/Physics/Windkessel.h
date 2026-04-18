/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Windkessel.h
 * @brief Outflow law toward the arterial Windkessel branch.
 */
#ifndef RODIN_HEART_CCMLC2014_PHYSICS_WINDKESSEL_H
#define RODIN_HEART_CCMLC2014_PHYSICS_WINDKESSEL_H

namespace Rodin::Heart::CCMLC2014::Physics
{
  /**
   * @brief Compute Windkessel outflow and midpoint-pressure derivatives.
   */
  template <class Input, class EvalData>
  inline void computeWindkesselOutflow(const Input& in, EvalData& evalData)
  {
    using Scalar = decltype(evalData.y);
    const Scalar aorticFlowRate = in.Kar * (evalData.pvMid - evalData.parMid);
    evalData.windkesselOutflow = (aorticFlowRate > Scalar(0)) ? aorticFlowRate : Scalar(0);

    const Scalar midpointOutflowDerivative =
      (aorticFlowRate > Scalar(0)) ? (Scalar(0.5) * in.Kar) : Scalar(0);
    evalData.dWindkesselOutflow_dPv = midpointOutflowDerivative;
    evalData.dWindkesselOutflow_dPar = -midpointOutflowDerivative;
  }
}

#endif
