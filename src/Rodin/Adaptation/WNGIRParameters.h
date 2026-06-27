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
  struct WNGIRParameters
  {
    Real h = 0;                 ///< reference mesh size (required).
    Real gammaM = 0;            ///< L² weight; ≤0 ⇒ 1/h.
    Real gammaH = 0;            ///< H¹ weight; ≤0 ⇒ 1/h.
    Real ellM = 0;              ///< Sobolev length; ≤0 ⇒ 3h.
    Real gammaObs = 1;          ///< surface observation metric weight.
    bool residualStabilizedObservationMetric = true;
    Real gammaJ = 1;            ///< j-barrier weight.
    Real gammaQ = 1;            ///< Q-barrier weight.
    Real jSafe = 1e-2;          ///< barrier floor on j.
    Real qMax = 10;             ///< barrier + line-search ceiling on Q.
    Real s0J = 0.25;            ///< j-barrier activation width.
    Real s0Q = 2;               ///< Q-barrier activation width.
    /// Hinge quality regularizer (energy, on the RHS). Penalizes only
    /// cells whose relative distortion exceeds a *good-quality*
    /// threshold qStar (chosen ≪ qMax):
    ///   E_q = (gammaQual/2) ∫_Ω max(0, Q_rel − qStar)² dX.
    /// ρ'(Q)=max(0,Q−qStar) is the redistribution force; ρ''(Q)=1 the
    /// Gauss–Newton Hessian (added to the metric). The energy is FLAT
    /// for well-shaped cells (Q ≤ qStar) — zero force, so it does not
    /// fight the fit — and pushes the degrading tail back to qStar.
    /// Its fixed point is the fit-optimal configuration subject to "no
    /// element worse than qStar", reached without touching the cells
    /// already good; displaced strain flows to under-threshold
    /// neighbours (which offer no resistance), so strain redistributes
    /// implicitly. The barrier (metric + line search) remains the hard
    /// validity wall underneath. gammaQual ≤ 0 disables.
    ///
    /// Default preset (on): gammaQual=1 is an O(1) trade-off weight
    /// against the data term; qStar=1.75 is a dimensionless distortion
    /// threshold (Q_rel=1 is ideal, ≪ qMax=10 the hard wall) and is
    /// h-independent — refining the mesh does not rescale it.
    Real gammaQual = 1;
    Real qStar = Real(1.75);
    /// Size hinge (energy, on the RHS). Companion to the shape hinge:
    /// penalizes cells whose Jacobian drops below a *good* level jStar
    /// (chosen ≫ jSafe):
    ///   E_s = (gammaSize/2) ∫_Ω max(0, jStar − j)² dX.
    /// ρ'(j) = −max(0, jStar − j) (restoring force pushing j up toward
    /// jStar); ρ''(j) = 1 for j < jStar (GN Hessian → metric). Flat for
    /// j ≥ jStar (no fit-fighting). Unlike the shape hinge on Q_rel —
    /// which is scale-invariant and undefined for j ≤ 0 — the size
    /// hinge depends only on j and the cofactor action a_j, so it is
    /// evaluated for ALL cells, INCLUDING inverted ones (j ≤ 0), giving
    /// a strong restoring force that pulls inverted elements back to
    /// validity. The metric barrier + line search remain the hard wall.
    /// gammaSize ≤ 0 disables.
    ///
    /// Default preset (on): gammaSize=1 (O(1) trade-off weight);
    /// jStar=0.3 is a dimensionless volume-ratio threshold (j=1 at
    /// identity, ≫ jSafe=1e-2) and is h-independent.
    Real gammaSize = 1;
    Real jStar = Real(0.3);
    Real omegaMin = 0.1;        ///< active-set threshold on ω.
    Real alphaMin = 1e-4;       ///< line-search floor.
    bool energyLineSearch = true;
    Real jMinRatio = 1e-8;      ///< hard inadmissibility floor.
    Real jLineSearchRatio = 1e-2;
    Real activeRMSTol = 0;      ///< ≤0 ⇒ 4h².
    Real activeSupTol = 0;      ///< ≤0 ⇒ 10h².
    Real energyStagTol = 1e-4;
    Real stepTol = 0;           ///< ≤0 ⇒ 1e-4·h.
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
    /// If true, also add the nonlinear barrier first variation to the
    /// RHS. Default false: admissibility enters as metric stiffness and
    /// true-geometry line-search validation, not as a competing energy.
    bool includeAdmissibilityGradient = false;
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
