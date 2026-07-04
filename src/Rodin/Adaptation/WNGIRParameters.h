/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_WNGIRPARAMETERS_H
#define RODIN_ADAPTATION_WNGIRPARAMETERS_H

#include "WNGIRCommon.h"

namespace Rodin::Adaptation
{
  enum class WNGIRMetricActivation
  {
    Hard,
    Smooth
  };

  struct WNGIRParameters
  {
    Real h = 0;                 ///< reference mesh size (required).
    Real gammaM = 0;            ///< L² weight; ≤0 ⇒ 1/h.
    Real gammaH = 0;            ///< deviatoric-strain weight; ≤0 ⇒ 1/h.
    Real gammaDiv = 0;          ///< divergence weight; ≤0 ⇒ gammaH.
    Real ellM = 0;              ///< Sobolev length; ≤0 ⇒ 3h.
    Real gammaObs = 1;          ///< surface observation metric weight.
    bool residualStabilizedObservationMetric = true;
    Real initialGuessGamma = 1000;
    Real initialGuessCapH = 2;
    Real gammaJ = 1;            ///< j-barrier weight.
    Real gammaQ = 1;            ///< Q-barrier weight.
    Real jSafe = 1e-2;          ///< barrier floor on j.
    Real qMax = 10;             ///< barrier + line-search ceiling on Q.
    Real s0J = 0.25;            ///< j-barrier activation width.
    Real s0Q = 2;               ///< Q-barrier activation width.
    /// One-sided relative-distortion quality metric:
    ///   K_Q(v,z) = gammaQual ∫_{Q_rel>qStar} a_Q(v) a_Q(z) dX.
    /// This changes the Riesz metric only; no quality force is added to the
    /// right-hand side.
    /// gammaQual ≤ 0 disables the Q hinge.
    Real gammaQual = 1;
    Real qStar = Real(1.75);
    /// Optional one-sided size hinge:
    ///   K_j(v,z) = gammaSize ∫_{j<jStar} a_j(v) a_j(z) dX.
    /// Disabled by default. Inversion is handled by the near-zero j barrier
    /// and the true-geometry line search, allowing small well-shaped cells.
    Real gammaSize = 0;
    Real jStar = Real(0.3);
    WNGIRMetricActivation metricActivation = WNGIRMetricActivation::Hard;
    Real qualitySmoothDelta = Real(0.1);
    Real jBarrierSmoothDelta = Real(0.05);
    Real qBarrierSmoothDelta = Real(0.1);
    Real metricActivationEpsilon = Real(1e-8);
    Real omegaMin = 0.1;        ///< active-set threshold on ω.
    Real alphaMin = 1e-4;       ///< line-search floor.
    bool energyLineSearch = true;
    Real jMinRatio = 1e-8;      ///< hard inadmissibility floor.
    Real jLineSearchRatio = 1e-2;
    Real activeRMSTol = 0;      ///< ≤0 ⇒ 4h².
    Real activeSupTol = 0;      ///< ≤0 ⇒ 10h².
    Real activeRMSOverHTol = 0; ///< >0 enables scale-aware RMS stopping.
    Real activeSupOverHTol = 0; ///< >0 enables scale-aware sup stopping.
    bool geometryAwareTolerances = true;
    Real rmsFloor2D = Real(0.05);
    Real supFloor2D = Real(0.25);
    Real rmsFloor3D = Real(0.03);
    Real supFloor3D = Real(0.20);
    Real rmsNormalJumpFactor = Real(0.03);
    Real supNormalJumpFactor = Real(0.05);
    Real energyStagTol = 1e-4;
    Real stepTol = 0;           ///< ≤0 ⇒ 1e-4·h.
    Real acceptedStepOverHTol = Real(5e-3); ///< >0 stops best-effort when accepted step/h is small.
    Real pointLocationTolerance = 0; ///< ≤0 ⇒ 1e-10 for moved FE evaluation.
    Real cgRelativeTolerance = 1e-6; ///< relative residual tolerance for CG.
    std::size_t cgMaxIterations = 0; ///< 0 ⇒ min(2000, max(100, 2*ndofs)).
    std::size_t andersonMemory = 3;  ///< 0 disables safeguarded Anderson.
    std::size_t andersonStart = 2;   ///< first iteration where AA may be used.
    Real andersonDamping = 1;        ///< first AA damping trial.
    Real andersonMinDamping = 0.125; ///< smallest AA damping trial.
    std::size_t maxIterations = 200;
    std::size_t quadratureOrder = 0; ///< 0 ⇒ 2·(FE order).
    bool hasInterfaceAttribute = false;
    Geometry::Attribute interfaceAttribute = 0;
    bool trace = false;
    /// If true, add near-boundary admissibility barriers to the metric.
    /// The true-geometry line search remains the final admissibility check.
    bool includeAdmissibilityMetric = true;
    /// If true, add the Q_rel and optional j-size hinge Gauss--Newton terms
    /// to the metric.
    bool includeQualityMetric = true;
    /// Optimal 1-D rescale of the lifted step along itself:
    ///   β = ⟨d, v⟩_Γ / ⟨v, v⟩_Γ  (surface inner products),
    /// clamped to [1, betaMax]; line search starts at β·v instead of
    /// v. The H¹ lift systematically under-scales the skeleton trace
    /// (gain ≈ surface-weight / M-diagonal ≈ 1/20 at default γ), so
    /// without β the iteration is linearly convergent with ρ ≈ 0.95.
    /// β recovers Newton-matched magnitude while preserving the lift's
    /// smooth admissibility-aware shape. Since β only scales the same
    /// descent direction, the nonlinear line search remains the final
    /// admissibility and energy-decrease guard.
    Real betaMax = 50;
  };
}

#endif
