/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_SDFRPARAMETERS_H
#define RODIN_ADAPTATION_SDFRPARAMETERS_H

#include "Rodin/Types.h"

#include "SDFRIntegrators.h"

namespace Rodin::Adaptation
{
  enum class SDFRInitialGuess
  {
    Zero,
    Current,
    Hilbert
  };

  enum class SDFRTangent
  {
    GaussNewton,
    Newton,
    PSDProjectedNewton
  };

  enum class SDFRHilbertMetric
  {
    Harmonic,
    Elasticity,
    ShapeHessian
  };

  enum class SDFRInitialGuessScaling
  {
    Unnormalized,
    EnergyNormalized,
    BandNormalized
  };

  struct SDFRParameters
  {
    Real rhoS = 1;
    Real deltaW = 0;
    Real hRef = 0;
    Real normalizer = 0; ///< <= 0 means compute 1 / (M_w h_ref^2).

    Real shapeWeight = 1e-1;
    Real floorWeight = 1e-2;

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

    Real harmonicEps = 0;
    Real tikhonovEps = 0;

    SDFRInitialGuess initialGuess = SDFRInitialGuess::Hilbert;
    SDFRHilbertMetric initialGuessMetric = SDFRHilbertMetric::Harmonic;
    SDFRInitialGuessScaling initialGuessScaling =
      SDFRInitialGuessScaling::Unnormalized;
    Real initialGuessElasticityLambda = 0;
    Real initialGuessElasticityMu = 0.5;
    SDFRTangent tangent = SDFRTangent::PSDProjectedNewton;
  };

  inline SDFRIntegratorParameters
  makeSDFRIntegratorParameters(const SDFRParameters& params)
  {
    SDFRIntegratorParameters out;
    out.rhoS = params.rhoS;
    out.deltaW = params.deltaW;
    out.hRef = params.hRef;
    out.normalizer = params.normalizer;
    return out;
  }

  inline SDFRIntegratorParameters makeSDRParameters(const SDFRParameters& params)
  {
    return makeSDFRIntegratorParameters(params);
  }
}

#endif
