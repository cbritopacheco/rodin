/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_LSRPARAMETERS_H
#define RODIN_ADAPTATION_LSRPARAMETERS_H

#include "Rodin/Types.h"

#include "LSRIntegrators.h"

namespace Rodin::Adaptation
{
  enum class LSRInitialGuess
  {
    Zero,
    Current,
    Hilbert
  };

  enum class LSRTangent
  {
    GaussNewton,
    Newton,
    PSDProjectedNewton
  };

  enum class LSRHilbertMetric
  {
    Harmonic,
    Elasticity,
    ShapeHessian
  };

  enum class LSRInitialGuessScaling
  {
    Unnormalized,
    EnergyNormalized,
    BandNormalized
  };

  struct LSRParameters
  {
    Real rhoS = 1;
    Real deltaW = 0;
    Real hRef = 0;
    Real normalizer = 0; ///< <= 0 means compute 1 / (M_w h_ref^2).

    Real shapeWeight = 1e-1;

    /// Smooth Q-shape barrier weight (`qBarrierWeight = 0` disables).
    Real qBarrierWeight = 0;
    /// Activation threshold; default +inf disables.
    Real qBarrierAct = std::numeric_limits<Real>::infinity();
    /// Asymptote; default +inf disables.
    Real qBarrierMax = std::numeric_limits<Real>::infinity();

    // jMinRatio: dimensionless singularity floor for the line-search
    // admissibility check + the defensive singular-cell fallback.
    // jSafeRatio: the line-search safety threshold; the line search
    // requires min_K j_K^u > lineSearchSafetyMargin * jSafeRatio.
    Real jMinRatio = 1e-8;
    Real jSafeRatio = 1e-3;
    Real lineSearchSafetyMargin = 10;

    Real alphaInit = 1;
    Real alphaReduction = 0.5;
    Real alphaMin = 1e-6;

    Real absoluteTolerance = 1e-8;
    Real relativeTolerance = 1e-7;
    Real stepTolerance = 0;
    std::size_t maxNewtonIterations = 40;
    std::size_t stallPatience = 3;

    // Best-effort quality safeguard. Line search rejects any trial step
    // whose worst-case cell would have intrinsic shape quality
    // Q_shape = ||A_K^u||_F^2 / (d * (sigma_K det A_K^u)^(2/d))
    // exceeding `qShapeMax`. The default infinity disables the cap.
    // A finite cap turns Newton into a "best-effort" iteration: it
    // proceeds as long as the mesh stays below the quality threshold,
    // and halts (returning the best-residual iterate) when no further
    // step satisfies it.
    Real qShapeMax = std::numeric_limits<Real>::infinity();

    LSRInitialGuess initialGuess = LSRInitialGuess::Hilbert;
    LSRHilbertMetric initialGuessMetric = LSRHilbertMetric::Harmonic;
    LSRInitialGuessScaling initialGuessScaling =
      LSRInitialGuessScaling::Unnormalized;
    Real initialGuessElasticityLambda = 0;
    Real initialGuessElasticityMu = 0.5;
    LSRTangent tangent = LSRTangent::PSDProjectedNewton;
  };

  inline LSRIntegratorParameters
  makeLSRIntegratorParameters(const LSRParameters& params)
  {
    LSRIntegratorParameters out;
    out.rhoS = params.rhoS;
    out.deltaW = params.deltaW;
    out.hRef = params.hRef;
    out.normalizer = params.normalizer;
    return out;
  }

  inline LSRIntegratorParameters makeLSRParameters(const LSRParameters& params)
  {
    return makeLSRIntegratorParameters(params);
  }
}

#endif
