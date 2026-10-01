/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_WNGIRDIRECTIONALNEWTON_H
#define RODIN_ADAPTATION_WNGIRDIRECTIONALNEWTON_H

#include <algorithm>
#include <cmath>
#include "Rodin/Types.h"

namespace Rodin::Adaptation::Detail
{
  /// Outer quadratic-model minimizer; action=-DE[v], curvature=D2E[v,v].
  /// Invalid/nonpositive curvature or an unusably small proposal keeps alpha=1.
  inline Real wngirDirectionalNewtonStep(Real action, Real curvature,
                                         Real alphaMin, Real alphaMax)
  {
    if (!(action > Real(0)) || !(curvature > Real(0)) ||
        !std::isfinite(action) || !std::isfinite(curvature))
      return Real(1);
    const Real proposal = action / curvature;
    if (!(proposal >= alphaMin) || !std::isfinite(proposal))
      return Real(1);
    return std::min(proposal, alphaMax);
  }
}
#endif
