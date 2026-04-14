/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_HEART_CCMLC2014_PASSIVELAW_H
#define RODIN_HEART_CCMLC2014_PASSIVELAW_H

#include <cmath>

#include "Rodin/Types.h"

namespace Rodin::Heart::CCMLC2014Laws
{
  struct PassiveResponse
  {
    Real dW1 = 0.0;
    Real dW2 = 0.0;
    Real dW4 = 0.0;
    Real dSigmaPassive_dC = 0.0;
  };

  /**
   * @brief Default Holzapfel-type passive response provider.
   *
   * This law is intentionally lightweight and only exposes the terms required
   * by the 0D spherical model equations.
   */
  struct HolzapfelOgdenLaw
  {
    struct Input
    {
      Real mu1 = 0.0;
      Real mu2 = 0.0;
      Real C0 = 0.0;
      Real C1 = 0.0;
      Real C2 = 0.0;
      Real C3 = 0.0;
    };

    Input input;

    PassiveResponse operator()(Real C, Real j1, Real /*j2*/, Real j4) const
    {
      const Real dj1 = j1 - 3.0;
      const Real dj4 = j4 - 1.0;
      const Real expj1 = std::exp(input.C1 * dj1 * dj1);
      const Real expj4 = std::exp(input.C3 * dj4 * dj4);

      PassiveResponse out;
      out.dW1 = input.mu1 + 2.0 * input.C0 * input.C1 * dj1 * expj1;
      out.dW2 = input.mu2;
      out.dW4 = 2.0 * input.C2 * input.C3 * dj4 * expj4;

      // Frozen-derivative approximation in terms of dW's (consistent with
      // previous model behavior and required Jacobian term).
      out.dSigmaPassive_dC =
        12.0 * std::pow(C, -4) * (out.dW1 + C * out.dW2)
        + 4.0 * (1.0 - std::pow(C, -3)) * out.dW2;
      return out;
    }
  };
}

#endif
