/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file PassiveLaw.h
 * @brief Conversion from reduced passive energy derivatives to 0D passive stress terms.
 */
#ifndef RODIN_HEART_CCMLC2014_PASSIVELAW_H
#define RODIN_HEART_CCMLC2014_PASSIVELAW_H

#include <cmath>

#include "Rodin/Heart/CCMLC2014/PassiveEnergy.h"

namespace Rodin::Heart
{
  /**
   * @brief Default passive stress law for the CCMLC2014 reduced model.
   *
   * @tparam Scalar Scalar numeric type.
   */
  template <class Scalar>
  struct CCMLC2014PassiveLaw
  {
    /**
     * @brief Evaluate passive stress and tangent with respect to radial displacement.
     *
     * @tparam PassiveEnergyLaw Passive-energy evaluator type.
     * @param[in] law Passive-energy law.
     * @param[in] C Reduced right Cauchy-Green scalar.
     * @param[in] dC_dy Derivative of @f$ C @f$ wrt displacement.
     * @param[out] sigmaPassive Passive stress.
     * @param[out] dsigmaPassive_dy Derivative of passive stress wrt displacement.
     */
    template <class PassiveEnergyLaw>
    void operator()(
        const PassiveEnergyLaw& law,
        const Scalar C,
        const Scalar dC_dy,
        Scalar& sigmaPassive,
        Scalar& dsigmaPassive_dy) const
    {
      const Scalar J1 = Scalar(2) * C + std::pow(C, -2);
      const Scalar J2 = C * C + Scalar(2) * std::pow(C, -1);
      const Scalar J4 = C;

      const Scalar dJ1_dC = Scalar(2) - Scalar(2) * std::pow(C, -3);
      const Scalar dJ2_dC = Scalar(2) * C - Scalar(2) * std::pow(C, -2);
      const Scalar dJ4_dC = Scalar(1);

      const ReducedInvariants<Scalar> I{J1, J2, J4};
      const auto D = law.evaluate(I);

      const Scalar dW1 = D.grad.dW_dJ1;
      const Scalar dW2 = D.grad.dW_dJ2;
      const Scalar dW4 = D.grad.dW_dJ4;

      sigmaPassive =
          Scalar(4) * (Scalar(1) - std::pow(C, -3)) * (dW1 + C * dW2)
        + Scalar(2) * dW4;

      const Scalar ddW1_dC =
          D.hess.d2W_dJ1dJ1 * dJ1_dC
        + D.hess.d2W_dJ1dJ2 * dJ2_dC
        + D.hess.d2W_dJ1dJ4 * dJ4_dC;

      const Scalar ddW2_dC =
          D.hess.d2W_dJ1dJ2 * dJ1_dC
        + D.hess.d2W_dJ2dJ2 * dJ2_dC
        + D.hess.d2W_dJ2dJ4 * dJ4_dC;

      const Scalar ddW4_dC =
          D.hess.d2W_dJ1dJ4 * dJ1_dC
        + D.hess.d2W_dJ2dJ4 * dJ2_dC
        + D.hess.d2W_dJ4dJ4 * dJ4_dC;

      const Scalar dsigmaPassive_dC =
          Scalar(12) * std::pow(C, -4) * (dW1 + C * dW2)
        + Scalar(4) * (Scalar(1) - std::pow(C, -3)) * (ddW1_dC + dW2 + C * ddW2_dC)
        + Scalar(2) * ddW4_dC;

      dsigmaPassive_dy = dsigmaPassive_dC * dC_dy;
    }
  };
}

#endif
