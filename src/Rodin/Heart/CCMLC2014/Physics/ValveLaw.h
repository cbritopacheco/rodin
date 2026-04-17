/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_HEART_CCMLC2014_PHYSICS_VALVELAW_H
#define RODIN_HEART_CCMLC2014_PHYSICS_VALVELAW_H

namespace Rodin::Heart::CCMLC2014::Physics
{
  template <class Input, class EvalData>
  inline void computeCavityFluxes(const Input& in, EvalData& d)
  {
    using Scalar = decltype(d.y);
    const bool mitralOpenCur = d.pv <= d.pAtCur;
    const bool bothClosedCur = d.pAtCur <= d.pv && d.pv <= d.par;

    if (mitralOpenCur)
    {
      d.cavityFluxCur = in.Kat * (d.pv - d.pAtCur);
      d.dCavityFluxCur_dPv = in.Kat;
      d.dCavityFluxCur_dPar = Scalar(0);
    }
    else if (bothClosedCur)
    {
      d.cavityFluxCur = in.Kp * (d.pv - d.pAtCur);
      d.dCavityFluxCur_dPv = in.Kp;
      d.dCavityFluxCur_dPar = Scalar(0);
    }
    else
    {
      d.cavityFluxCur = in.Kar * (d.pv - d.par) + in.Kp * (d.par - d.pAtCur);
      d.dCavityFluxCur_dPv = in.Kar;
      d.dCavityFluxCur_dPar = -in.Kar + in.Kp;
    }

    const bool mitralOpenPrev = d.pvPrev <= d.pAtPrev;
    const bool bothClosedPrev = d.pAtPrev <= d.pvPrev && d.pvPrev <= d.parPrev;

    if (mitralOpenPrev)
    {
      d.cavityFluxPrev = in.Kat * (d.pvPrev - d.pAtPrev);
    }
    else if (bothClosedPrev)
    {
      d.cavityFluxPrev = in.Kp * (d.pvPrev - d.pAtPrev);
    }
    else
    {
      d.cavityFluxPrev = in.Kar * (d.pvPrev - d.parPrev) + in.Kp * (d.parPrev - d.pAtPrev);
    }
  }
}

#endif
