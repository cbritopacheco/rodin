/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_HEART_CCMLC2014_PASSIVELAW_H
#define RODIN_HEART_CCMLC2014_PASSIVELAW_H

#include <cmath>

namespace Rodin::Heart
{
  template <class Scalar>
  struct CCMLC2014PassiveLaw
  {
    template <class Input>
    void operator()(
        const Input& input,
        const Scalar C,
        const Scalar dC_dy,
        Scalar& sigmaPassive,
        Scalar& dsigmaPassive_dy) const
    {
      const Scalar J1 = Scalar(2) * C + std::pow(C, -2);
      const Scalar J2 = C * C + Scalar(2) * std::pow(C, -1);
      const Scalar J4 = C;

      const Scalar dW1 = input.dWe_dJ1(J1, J2, J4);
      const Scalar dW2 = input.dWe_dJ2(J1, J2, J4);
      const Scalar dW4 = input.dWe_dJ4(J1, J2, J4);

      sigmaPassive =
          Scalar(4) * (Scalar(1) - std::pow(C, -3)) * (dW1 + C * dW2)
        + Scalar(2) * dW4;

      dsigmaPassive_dy =
        input.dSigmaPassive_dC(C, J1, J2, J4, dW1, dW2, dW4) * dC_dy;
    }
  };
}

#endif
