/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_HEART_CCMLC2014_PHYSICS_WINDKESSEL_H
#define RODIN_HEART_CCMLC2014_PHYSICS_WINDKESSEL_H

namespace Rodin::Heart::CCMLC2014::Physics
{
  template <class Input, class EvalData>
  inline void computeWindkesselOutflow(const Input& in, EvalData& d)
  {
    using Scalar = decltype(d.y);
    const Scalar flowAr = in.Kar * (d.pvMid - d.parMid);
    d.windkesselOutflow = (flowAr > Scalar(0)) ? flowAr : Scalar(0);

    const Scalar difFlow = (flowAr > Scalar(0)) ? (Scalar(0.5) * in.Kar) : Scalar(0);
    d.dWindkesselOutflow_dPv = difFlow;
    d.dWindkesselOutflow_dPar = -difFlow;
  }
}

#endif
