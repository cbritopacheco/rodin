/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_LSRPARAMETERS_H
#define RODIN_ADAPTATION_LSRPARAMETERS_H

#include <functional>

#include "Rodin/Geometry/Types.h"
#include "Rodin/Types.h"

#include "LSRIntegrators.h"
#include "LSRReport.h"

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
    std::size_t quadratureOrder = 0; ///< 0 selects the default FE-based rule.
    LSRIntegratorParameters::FieldEvaluation fieldEvaluation =
      LSRIntegratorParameters::FieldEvaluation::PhysicalPoint;

    /// Weight for the relative distortion energy Q_rel(I + grad u) - 1.
    Real shapeWeight = 1e-1;
    Real h1RegularizationWeight = 0;

    /// Interface-distance term
    ///   E_Gamma = 0.5 * interfaceWeight * interfaceNormalizer
    ///             * int_{Gamma_h} (phi(X + u) / |grad phi(X + u)|)^2 ds.
    /// This term does not assume that phi is a signed distance function.
    Real interfaceWeight = 0;
    Geometry::Attribute interfaceAttribute = 0;
    Real interfaceNormalizer = 0; ///< <= 0 means compute 1 / (|Gamma_h| h_ref^2).
    Real interfaceGradientFloor = 1e-12;

    /// Band-weighted Hilbert metric stiffness:
    ///   a(u, v) = int_Omega (1 + s_H * (1 - W(psi))^2) grad u : grad v dX.
    /// With s_H = 0 the Harmonic metric is recovered. With s_H > 0 the
    /// Riesz lift used by the Hilbert initial guess is stiffer outside
    /// the data band W(psi), so the lift is localised to the band and the
    /// harmonic-extension extremum in the geometric centre of the
    /// classified interior is suppressed. Only affects LSRHilbertMetric::
    /// Harmonic; ignored for Elasticity and ShapeHessian.
    Real bandHilbertStiffness = 0;

    /// Outside-band L^2 Tikhonov damping
    ///   E_damp = 0.5 * outsideBandTikhonovWeight
    ///            * int_Omega (1 - W(psi(X)))^2 |u|^2 dX,
    ///   W(s) = exp(-s^2 / (2 deltaW^2)).
    /// Localises u to the data band by penalising L^2 mass of u where the
    /// LSR data weight W(psi) is small (interior + far exterior). Cures
    /// the harmonic-extension bulge in the geometric centre of the
    /// classified interior that the Hilbert initial guess otherwise
    /// produces. Identity-neutral (zero at u = 0).
    Real outsideBandTikhonovWeight = 0;

    /// Smooth relative-Q barrier weight (`qBarrierWeight = 0` disables).
    Real qBarrierWeight = 0;
    /// Activation threshold; default +inf disables.
    Real qBarrierAct = std::numeric_limits<Real>::infinity();
    /// Asymptote; default +inf disables.
    Real qBarrierMax = std::numeric_limits<Real>::infinity();

    /// Smooth admissibility barrier on the dimensionless Jacobian ratio
    /// j = det(I + grad u_h). Inactive when the weight is zero or the
    /// safe ratio is not larger than jMinRatio.
    Real jBarrierWeight = 0;
    Real jBarrierSafeRatio = 0;

    /// Centered volume tether 0.5 * (log j)^2 on the same dimensionless
    /// ratio. The residual is zero at j = 1, but the tangent is active.
    Real jVolumeTetherWeight = 0;

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
    Real energyDecreaseTolerance = 1e-12;

    /// Warm-start the next iteration's initial step at
    /// alphaWarmStartGrowth * (previously accepted alpha), clamped to
    /// [alphaMin, 1]. Enabled by default: on hard cases this cuts the
    /// total backtrack count by ~50% with no convergence cost; on easy
    /// cases the previously accepted alpha is 1 and warm-start is a
    /// no-op.
    bool useWarmStartAlpha = true;
    Real alphaWarmStartGrowth = 2.0;
    bool useSampledQuadraticAlphaPredictor = true;
    Real alphaPredictorSafety = 0.9;

    Real absoluteTolerance = 1e-8;
    Real relativeTolerance = 1e-7;
    Real stepTolerance = 0;
    std::size_t maxNewtonIterations = 40;
    std::size_t stallPatience = 3;

    // Optional problem-level stop criterion evaluated after an admissible
    // accepted step, with the displacement already updated to that step.
    // This lets callers stop on geometric criteria without exposing their
    // interface diagnostics inside the LSR solver.
    std::function<bool(const LSRReport&)> acceptedStateConvergenceTest;

    // Optional cap on the relative-distortion quality of the moved cell:
    //   Q_rel(F) = ||F||^2 / (d * det(F)^(2/d)),  F = I + grad u_h.
    // Identity-neutral (Q_rel = 1 at u = 0, regardless of background shape).
    // The line search rejects any step that pushes max_K Q_rel above this
    // cap. Default +infinity disables the cap. A value of ~1.5 means
    // "no cell may be distorted by more than ~50% from its background
    // shape during this Newton step".
    Real qRelMax = std::numeric_limits<Real>::infinity();

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
    out.quadratureOrder = params.quadratureOrder;
    out.fieldEvaluation = params.fieldEvaluation;
    out.interfaceWeight = params.interfaceWeight;
    out.interfaceAttribute = params.interfaceAttribute;
    out.interfaceNormalizer = params.interfaceNormalizer;
    out.interfaceGradientFloor = params.interfaceGradientFloor;
    return out;
  }

  inline LSRIntegratorParameters makeLSRParameters(const LSRParameters& params)
  {
    return makeLSRIntegratorParameters(params);
  }
}

#endif
