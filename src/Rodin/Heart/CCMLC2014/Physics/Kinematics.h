/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_HEART_CCMLC2014_PHYSICS_KINEMATICS_H
#define RODIN_HEART_CCMLC2014_PHYSICS_KINEMATICS_H

#include <cassert>

namespace Rodin::Heart::CCMLC2014::Physics
{
  template <class Input, class EvalData>
  inline void computeMidpointKinematics(const Input& in, EvalData& d)
  {
    using Scalar = decltype(d.y);
    d.yMid = Scalar(0.5) * (d.y + d.yPrev);
    d.pvMid = Scalar(0.5) * (d.pv + d.pvPrev);
    d.parMid = Scalar(0.5) * (d.par + d.parPrev);
    d.pdMid = Scalar(0.5) * (d.pd + d.pdPrev);
    d.vel = (d.y - d.yPrev) / d.dt;

    d.sqrtC = Scalar(1) + d.yMid / in.R0;
    assert(d.sqrtC > Scalar(0));
    d.C = d.sqrtC * d.sqrtC;

    d.strain1D = Scalar(0.5) * (d.C - Scalar(1));
    d.diffGreen = d.sqrtC / in.R0;
  }
}

#endif
