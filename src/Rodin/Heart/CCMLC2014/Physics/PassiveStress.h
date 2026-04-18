/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file PassiveStress.h
 * @brief Passive and viscous stress contributions for CCMLC2014 reduced dynamics.
 */
#ifndef RODIN_HEART_CCMLC2014_PHYSICS_PASSIVESTRESS_H
#define RODIN_HEART_CCMLC2014_PHYSICS_PASSIVESTRESS_H

#include <cmath>

namespace Rodin::Heart::CCMLC2014::Physics
{
  /**
   * @brief Evaluate passive stress and tangent contribution from the selected passive law.
   */
  template <class PassiveLaw, class Input, class EvalData>
  inline void computePassiveContribution(const Input& in, EvalData& evalData)
  {
    using Scalar = decltype(evalData.y);
    const Scalar dC_dyMid = Scalar(2) * evalData.sqrtC / in.R0;

    Scalar passiveStress = 0.0;
    Scalar passiveStressDerivativeWrtMidpointDisplacement = 0.0;
    PassiveLaw passiveLaw;
    passiveLaw(
        in.passiveEnergy,
        evalData.C,
        dC_dyMid,
        passiveStress,
        passiveStressDerivativeWrtMidpointDisplacement);

    evalData.stressPassive = passiveStress;
    evalData.diffStressPassive = Scalar(0.5) * passiveStressDerivativeWrtMidpointDisplacement;
  }

  /**
   * @brief Evaluate viscous stress and tangent contribution.
   */
  template <class Input, class EvalData>
  inline void computeViscousContribution(const Input& in, EvalData& evalData)
  {
    using Scalar = decltype(evalData.y);
    const Scalar referenceRadius = in.R0;
    const Scalar viscosity = in.eta;
    const Scalar radialVelocity = evalData.vel;
    const Scalar reducedStretch = evalData.sqrtC;
    const Scalar reducedCauchyGreen = evalData.C;

    const Scalar dGreenRRWrtDisplacement = Scalar(-2) / referenceRadius * std::pow(reducedStretch, -5);
    const Scalar dGreenThetaThetaWrtDisplacement = evalData.diffGreen;

    const Scalar viscousStressRR = viscosity * dGreenRRWrtDisplacement * radialVelocity;
    const Scalar viscousStressThetaTheta = viscosity * dGreenThetaThetaWrtDisplacement * radialVelocity;
    const Scalar viscousStressPhiPhi = viscousStressThetaTheta;

    evalData.stressViscous =
      viscousStressThetaTheta + viscousStressPhiPhi
      - Scalar(2) * std::pow(reducedCauchyGreen, -3) * viscousStressRR;

    const Scalar dDotCRRWrtDisplacement =
      Scalar(10) / (referenceRadius * referenceRadius) * viscosity * std::pow(reducedStretch, -6) * radialVelocity
      + Scalar(2) * viscosity * dGreenRRWrtDisplacement / evalData.dt;

    const Scalar dDotCThetaThetaWrtDisplacement =
      viscosity / (referenceRadius * referenceRadius) * radialVelocity
      + Scalar(2) * viscosity * dGreenThetaThetaWrtDisplacement / evalData.dt;

    const Scalar viscousStressDerivativeWrtDisplacement =
      dDotCThetaThetaWrtDisplacement + dDotCThetaThetaWrtDisplacement
      - Scalar(2) * (std::pow(reducedCauchyGreen, -3) * dDotCRRWrtDisplacement
          - Scalar(6) / referenceRadius * std::pow(reducedStretch, -7) * viscousStressRR);

    evalData.diffStressViscous = viscousStressDerivativeWrtDisplacement;
  }
}

#endif
