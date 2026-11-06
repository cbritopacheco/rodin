/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_SWIFT_DIRECTIONALNEWTON_H
#define RODIN_ADAPTATION_SWIFT_DIRECTIONALNEWTON_H

#include <algorithm>
#include <cmath>
#include "Rodin/Types.h"

namespace Rodin::Adaptation::SWIFT
{
  /**
   * @brief Quadratic-model minimizer with an optional physical-motion bound.
   *
   * Nonpositive curvature selects the physical-motion bound when it is positive.
   * With no bound, positive curvature is required for a finite minimizing step.
   * @param action Positive fitting-force action on the predictor.
   * @param curvature Robust-energy directional curvature.
   * @param directionNorm Physical displacement norm of the predictor.
   * @param maximumStep Physical-motion bound; zero leaves scaling unrestricted.
   * @returns Directional scale, or zero for invalid inputs or an unbounded model.
   */
  inline Real getDirectionalNewtonStep(
    Real action, Real curvature, Real directionNorm, Real maximumStep)
  {
    if (!(action > Real(0)) || !std::isfinite(action) || !(directionNorm > Real(0)) ||
      !std::isfinite(directionNorm) || !std::isfinite(curvature) ||
      maximumStep < Real(0) || !std::isfinite(maximumStep))
      return Real(0);
    if (curvature <= Real(0))
      return maximumStep > Real(0) ? maximumStep / directionNorm : Real(0);
    const Real scale = action / curvature;
    return maximumStep > Real(0) ? std::min(scale, maximumStep / directionNorm) : scale;
  }
}
#endif
