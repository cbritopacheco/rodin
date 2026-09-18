/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Metrics.h"

namespace KelvinBall
{
  Values Metrics::evaluateChamber(Mesh& fluid) const
  {
    prepare(fluid);
    VelocitySpace velocitySpace(fluid, 3);
    PressureSpace pressureSpace(fluid);
    RotatedNitscheIntegrator coupling(fluid, MasterCuts, 0.01, 0.25);
    GridFunction uT0(velocitySpace), uT1(velocitySpace), uT2(velocitySpace),
      uR0(velocitySpace), uR1(velocitySpace), uR2(velocitySpace);
    GridFunction pT0(pressureSpace), pT1(pressureSpace), pT2(pressureSpace),
      pR0(pressureSpace), pR1(pressureSpace), pR2(pressureSpace);
    return evaluateChamber(velocitySpace, pressureSpace, coupling, uT0, uT1, uT2, uR0,
      uR1, uR2, pT0, pT1, pT2, pR0, pR1, pR2);
  }

  Values Metrics::evaluateSewed(Mesh& fluid) const
  {
    prepare(fluid);
    VelocitySpace velocitySpace(fluid, 3);
    PressureSpace pressureSpace(fluid);
    GridFunction uT0(velocitySpace), uT1(velocitySpace), uT2(velocitySpace),
      uR0(velocitySpace), uR1(velocitySpace), uR2(velocitySpace);
    GridFunction pT0(pressureSpace), pT1(pressureSpace), pT2(pressureSpace),
      pR0(pressureSpace), pR1(pressureSpace), pR2(pressureSpace);
    solveFullFamily(velocitySpace, pressureSpace, VectorFunction{1, 0, 0},
      VectorFunction{0, 1, 0}, VectorFunction{0, 0, 1}, uT0, uT1, uT2, pT0, pT1, pT2);
    solveFullFamily(velocitySpace, pressureSpace, VectorFunction{0, -F::z, F::y},
      VectorFunction{F::z, 0, -F::x}, VectorFunction{-F::y, F::x, 0}, uR0, uR1, uR2, pR0,
      pR1, pR2);
    return resistance(velocitySpace, uT0, uT1, uT2, uR0, uR1, uR2, 1);
  }
}
