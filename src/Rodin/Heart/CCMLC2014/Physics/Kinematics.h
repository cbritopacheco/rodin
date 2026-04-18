/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Kinematics.h
 * @brief Midpoint kinematics for the CCMLC2014 reduced 0D cavity model.
 */
#ifndef RODIN_HEART_CCMLC2014_PHYSICS_KINEMATICS_H
#define RODIN_HEART_CCMLC2014_PHYSICS_KINEMATICS_H

#include <cassert>

namespace Rodin::Heart::CCMLC2014::Physics
{
  /**
   * @brief Compute midpoint pressures/displacement and derived reduced kinematics.
   */
  template <class Input, class EvalData>
  inline void computeMidpointKinematics(const Input& in, EvalData& evalData)
  {
    using Scalar = decltype(evalData.y);
    evalData.yMid = Scalar(0.5) * (evalData.y + evalData.yPrev);
    evalData.pvMid = Scalar(0.5) * (evalData.pv + evalData.pvPrev);
    evalData.parMid = Scalar(0.5) * (evalData.par + evalData.parPrev);
    evalData.pdMid = Scalar(0.5) * (evalData.pd + evalData.pdPrev);
    evalData.vel = (evalData.y - evalData.yPrev) / evalData.dt;

    evalData.sqrtC = Scalar(1) + evalData.yMid / in.R0;
    assert(evalData.sqrtC > Scalar(0));
    evalData.C = evalData.sqrtC * evalData.sqrtC;

    evalData.strain1D = Scalar(0.5) * (evalData.C - Scalar(1));
    evalData.diffGreen = evalData.sqrtC / in.R0;
  }
}

#endif
