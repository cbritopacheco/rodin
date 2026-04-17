/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_HEART_CCMLC2014_PHYSICS_PASSIVESTRESS_H
#define RODIN_HEART_CCMLC2014_PHYSICS_PASSIVESTRESS_H

#include <cmath>

namespace Rodin::Heart::CCMLC2014::Physics
{
  template <class PassiveLaw, class Input, class EvalData>
  inline void computePassiveContribution(const Input& in, EvalData& d)
  {
    using Scalar = decltype(d.y);
    const Scalar dC_dyMid = Scalar(2) * d.sqrtC / in.R0;

    Scalar sigmaPassive = 0.0;
    Scalar dsigmaPassive_dyMid = 0.0;
    PassiveLaw passiveLaw;
    passiveLaw(in.passiveEnergy, d.C, dC_dyMid, sigmaPassive, dsigmaPassive_dyMid);

    d.stressPassive = sigmaPassive;
    d.diffStressPassive = Scalar(0.5) * dsigmaPassive_dyMid;
  }

  template <class Input, class EvalData>
  inline void computeViscousContribution(const Input& in, EvalData& d)
  {
    using Scalar = decltype(d.y);
    const Scalar R0 = in.R0;
    const Scalar nu = in.eta;
    const Scalar vel = d.vel;
    const Scalar sqrtC = d.sqrtC;
    const Scalar C = d.C;

    const Scalar diffGreen_rr = Scalar(-2) / R0 * std::pow(sqrtC, -5);
    const Scalar diffGreen_pp = d.diffGreen;

    const Scalar stress_rr = nu * diffGreen_rr * vel;
    const Scalar stress_pp = nu * diffGreen_pp * vel;
    const Scalar stress_tt = stress_pp;

    d.stressViscous =
      stress_pp + stress_tt - Scalar(2) * std::pow(C, -3) * stress_rr;

    const Scalar d_dotC_rr_dy =
      Scalar(10) / (R0 * R0) * nu * std::pow(sqrtC, -6) * vel
      + Scalar(2) * nu * diffGreen_rr / d.dt;

    const Scalar d_dotC_pp_dy =
      nu / (R0 * R0) * vel
      + Scalar(2) * nu * diffGreen_pp / d.dt;

    const Scalar diffStress =
      d_dotC_pp_dy + d_dotC_pp_dy
      - Scalar(2) * (std::pow(C, -3) * d_dotC_rr_dy
          - Scalar(6) / R0 * std::pow(sqrtC, -7) * stress_rr);

    d.diffStressViscous = diffStress;
  }
}

#endif
