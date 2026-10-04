/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_WNGIR_DIRECTIONALNEWTON_H
#define RODIN_ADAPTATION_WNGIR_DIRECTIONALNEWTON_H

#include <algorithm>
#include <cmath>
#include "Rodin/Types.h"

namespace Rodin::Adaptation
{
  /// @brief Quadratic-model minimizer, bounded by physical motion rather than raw alpha.
  /// Positive fitting curvature supplies the fallback when robust curvature is nonpositive.
  inline Real wngirDirectionalNewtonStep(Real action, Real curvature,
    Real fittingCurvature, Real directionNorm, Real maximumStep)
  {
    if (!(action > Real(0)) || !std::isfinite(action) || !(directionNorm > Real(0)) ||
      !std::isfinite(directionNorm) || !(maximumStep > Real(0)) ||
      !std::isfinite(maximumStep))
      return Real(0);
    const Real positiveCurvature =
      curvature > Real(0) && std::isfinite(curvature) ? curvature : fittingCurvature;
    if (!(positiveCurvature > Real(0)) || !std::isfinite(positiveCurvature))
      return Real(0);
    return std::min(action / positiveCurvature, maximumStep / directionNorm);
  }
}
#endif
