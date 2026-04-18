/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ValveLaw.h
 * @brief Valve branch laws used in the 0D cavity mass balance.
 */
#ifndef RODIN_HEART_CCMLC2014_PHYSICS_VALVELAW_H
#define RODIN_HEART_CCMLC2014_PHYSICS_VALVELAW_H

namespace Rodin::Heart::CCMLC2014::Physics
{
  /**
   * @brief Compute cavity inflow/outflow terms and their pressure derivatives.
   */
  template <class Input, class EvalData>
  inline void computeCavityFluxes(const Input& in, EvalData& evalData)
  {
    using Scalar = decltype(evalData.y);
    const bool mitralValveOpenAtCurrentStep = evalData.pv <= evalData.pAtCur;
    const bool bothValvesClosedAtCurrentStep =
      evalData.pAtCur <= evalData.pv && evalData.pv <= evalData.par;

    if (mitralValveOpenAtCurrentStep)
    {
      evalData.cavityFluxCur = in.Kat * (evalData.pv - evalData.pAtCur);
      evalData.dCavityFluxCur_dPv = in.Kat;
      evalData.dCavityFluxCur_dPar = Scalar(0);
    }
    else if (bothValvesClosedAtCurrentStep)
    {
      evalData.cavityFluxCur = in.Kp * (evalData.pv - evalData.pAtCur);
      evalData.dCavityFluxCur_dPv = in.Kp;
      evalData.dCavityFluxCur_dPar = Scalar(0);
    }
    else
    {
      evalData.cavityFluxCur = in.Kar * (evalData.pv - evalData.par) + in.Kp * (evalData.par - evalData.pAtCur);
      evalData.dCavityFluxCur_dPv = in.Kar;
      evalData.dCavityFluxCur_dPar = -in.Kar + in.Kp;
    }

    const bool mitralValveOpenAtPreviousStep = evalData.pvPrev <= evalData.pAtPrev;
    const bool bothValvesClosedAtPreviousStep =
      evalData.pAtPrev <= evalData.pvPrev && evalData.pvPrev <= evalData.parPrev;

    if (mitralValveOpenAtPreviousStep)
    {
      evalData.cavityFluxPrev = in.Kat * (evalData.pvPrev - evalData.pAtPrev);
    }
    else if (bothValvesClosedAtPreviousStep)
    {
      evalData.cavityFluxPrev = in.Kp * (evalData.pvPrev - evalData.pAtPrev);
    }
    else
    {
      evalData.cavityFluxPrev =
        in.Kar * (evalData.pvPrev - evalData.parPrev) + in.Kp * (evalData.parPrev - evalData.pAtPrev);
    }
  }
}

#endif
