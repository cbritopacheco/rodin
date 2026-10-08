/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_WNGIR_DIRECTIONALNEWTON_H
#define RODIN_ADAPTATION_WNGIR_DIRECTIONALNEWTON_H

#include <algorithm>
#include <cmath>
#include "Rodin/Types.h"

namespace Rodin::Adaptation::WNGIR
{
  /**
   * @brief Quadratic-model minimizer with an optional physical-motion bound.
   * Positive fitting curvature supplies the fallback when robust curvature is nonpositive.
   * A zero maximumStep leaves the directional Newton scale unrestricted.
   * @param action Positive fitting-force action on the predictor.
   * @param curvature Robust-energy directional curvature.
   * @param fittingCurvature Positive squared-fit curvature used as a fallback.
   * @param directionNorm Physical displacement norm of the predictor.
   * @param maximumStep Physical-motion bound; zero leaves scaling unrestricted.
   * @returns Directional scale, or zero when the model inputs are invalid.
   */
  inline Real getDirectionalNewtonStep(Real action, Real curvature,
    Real fittingCurvature, Real directionNorm, Real maximumStep)
  {
    if (!(action > Real(0)) || !std::isfinite(action) || !(directionNorm > Real(0)) ||
      !std::isfinite(directionNorm) || maximumStep < Real(0) ||
      !std::isfinite(maximumStep))
      return Real(0);
    const Real positiveCurvature =
      curvature > Real(0) && std::isfinite(curvature) ? curvature : fittingCurvature;
    if (!(positiveCurvature > Real(0)) || !std::isfinite(positiveCurvature))
      return Real(0);
    const Real scale = action / positiveCurvature;
    return maximumStep > Real(0) ? std::min(scale, maximumStep / directionNorm) : scale;
  }
}
#endif
