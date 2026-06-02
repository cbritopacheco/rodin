/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//
// Diagnostic + line-search helpers for `LevelSetLSRWavyCircleSweep.cpp`.
//
// Kept in a separate header to keep the example itself focused on the
// per-frame orchestration. The contents here are:
//
//   1. `AdmissibilityReport` — black-box (minJRatio,
//       inadmissibleCount) summary at the cell scale.
//   2. `evaluateAdmissibilityP1` — P1-vector admissibility check
//       computed in closed form from the existing cell-geometry cache.
//   3. `admissibilityLineSearch` — FES-INDEPENDENT backtracking line
//       search; only sees dof-vectors and an `Evaluator` callable.
//   4. `AdmissibilityDiagnostics` — X1-style per-iteration inspection
//       (min_j_ratio_before, min_j_ratio_after_full, inadm_after,
//       min_pred_margin, pred_inadm). Useful for understanding WHY a
//       step had to backtrack but never used to decide whether to accept.
//
#ifndef LEVELSETLSRWAVYCIRCLESWEEP_DIAGNOSTICS_H
#define LEVELSETLSRWAVYCIRCLESWEEP_DIAGNOSTICS_H

#include <cstddef>
#include <limits>
#include <vector>

#include "Rodin/Adaptation/CellGeomCache.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Types.h"

namespace LevelSetLSRWavyCircleSweepSupport
{
  using Rodin::Index;
  using Rodin::Real;
  namespace Math   = Rodin::Math;
  namespace Adapta = Rodin::Adaptation;

  // ---------- Black-box admissibility report --------------------------------
  struct AdmissibilityReport
  {
    // Dimensionless branch-normalized Jacobian ratio:
    //   j_K^u = sigma_K det(A_K^u) / J_K_scale.
    Real        minJRatio         = std::numeric_limits<Real>::infinity();
    std::size_t inadmissibleCount = 0;
  };

  // ---------- P1-vector admissibility evaluator -----------------------------
  //
  // For affine P1 triangles `j_K^u = sigma_K · det(F) · det(A_K) /
  // J_K_scale` is constant per cell, so a single per-cell evaluation
  // is both exact and cheap. Reads the 3-vertex displacement directly
  // through the FES's `getDOFs(0, vertex)` interface, identical to the
  // existing barrier integrator's convention.
  //
  // The line-search wrapper below is FES-independent and never calls
  // this directly — the example wires it in through a closure so a
  // higher-order or non-affine FES can substitute a different sampler
  // without touching the search loop.
  template <class FES>
  AdmissibilityReport evaluateAdmissibilityP1(
      const Math::Vector<Real>& uData,
      const std::vector<Adapta::CellGeomCache>& cellCache,
      const FES& fes,
      Real jMinRatio)
  {
    AdmissibilityReport rep;
    for (const auto& cellG : cellCache)
    {
      Math::SpatialMatrix<Real> U(2, 3);
      for (std::size_t i = 0; i < 3; ++i)
      {
        const Index v = cellG.vertices[i];
        const auto& dofs = fes.getDOFs(0, v);
        U(0, i) = uData(dofs[0]);
        U(1, i) = uData(dofs[1]);
      }
      const Math::SpatialMatrix<Real> F =
        Math::SpatialMatrix<Real>::Identity(2, 2) + U * cellG.gradN;
      const Real sigDetA =
        static_cast<Real>(cellG.sigmaK) * cellG.detAK;
      const Real j = sigDetA * F.determinant() / cellG.Jscale;
      if (j < rep.minJRatio)    rep.minJRatio         = j;
      if (j <= jMinRatio)       ++rep.inadmissibleCount;
    }
    return rep;
  }

  // ---------- FES-independent backtracking line search ----------------------
  struct LineSearchResult
  {
    Real        alphaAccepted             = Real(0);
    std::size_t backtracks                = 0;
    Real        minJRatioAtAlpha1         = std::numeric_limits<Real>::quiet_NaN();
    std::size_t inadmissibleCountAtAlpha1 = 0;
    Real        minJRatioAccepted         = std::numeric_limits<Real>::quiet_NaN();
    std::size_t inadmissibleCountAccepted = 0;
    bool        succeeded                 = false;
  };

  // The search wrapper only does
  //     uTrial = uOld + alpha * p
  //     adm    = evaluator(uTrial)
  // and never inspects the FES or element type. `Evaluator` is any
  // callable with signature `(Math::Vector<Real>) -> AdmissibilityReport`.
  template <class Evaluator>
  LineSearchResult admissibilityLineSearch(
      Math::Vector<Real>& u,
      const Math::Vector<Real>& p,
      Evaluator&& evaluator,
      Real jLineSearchRatio,
      Real alphaInit = Real(1),
      Real reduction = Real(0.5),
      Real alphaMin  = Real(1e-6))
  {
    LineSearchResult result;
    const Math::Vector<Real> uOld = u;

    Real alpha = alphaInit;
    bool firstIter = true;
    AdmissibilityReport lastAdm;

    while (alpha >= alphaMin)
    {
      const Math::Vector<Real> uTrial = uOld + alpha * p;
      const AdmissibilityReport adm = evaluator(uTrial);

      if (firstIter)
      {
        result.minJRatioAtAlpha1         = adm.minJRatio;
        result.inadmissibleCountAtAlpha1 = adm.inadmissibleCount;
        firstIter = false;
      }
      lastAdm = adm;

      if (adm.minJRatio > jLineSearchRatio
          && adm.inadmissibleCount == 0)
      {
        u = uTrial;
        result.alphaAccepted             = alpha;
        result.minJRatioAccepted         = adm.minJRatio;
        result.inadmissibleCountAccepted = adm.inadmissibleCount;
        result.succeeded                 = true;
        return result;
      }
      alpha *= reduction;
      ++result.backtracks;
    }

    // Failed: leave u unchanged.
    u = uOld;
    result.minJRatioAccepted         = lastAdm.minJRatio;
    result.inadmissibleCountAccepted = lastAdm.inadmissibleCount;
    return result;
  }

  // ---------- X1-style per-iter admissibility diagnostics -------------------
  //
  // Pure observation: never modifies the solve. The values are computed
  // from `(u_old, p_full)` where `p_full` is the un-damped Newton
  // direction. Useful to expose WHY a step had to backtrack.
  struct AdmissibilityDiagnostics
  {
    Real        minJRatioBefore     = std::numeric_limits<Real>::infinity();
    Real        minJRatioAfterFull  = std::numeric_limits<Real>::infinity();
    std::size_t inadmAfter     = 0;
    Real        minPredMargin  = std::numeric_limits<Real>::infinity();
    std::size_t predInadm      = 0;
  };

  template <class FES>
  AdmissibilityDiagnostics computeAdmissibilityDiagnostics(
      const Math::Vector<Real>& uOld,
      const Math::Vector<Real>& p,
      const std::vector<Adapta::CellGeomCache>& cellCache,
      const FES& fes,
      Real jMinRatio,
      Real tau)
  {
    AdmissibilityDiagnostics d;
    for (const auto& cellG : cellCache)
    {
      Math::SpatialMatrix<Real> U(2, 3), P(2, 3);
      for (std::size_t i = 0; i < 3; ++i)
      {
        const Index v = cellG.vertices[i];
        const auto& dofs = fes.getDOFs(0, v);
        U(0, i) = uOld(dofs[0]);
        U(1, i) = uOld(dofs[1]);
        P(0, i) = p   (dofs[0]);
        P(1, i) = p   (dofs[1]);
      }
      const Math::SpatialMatrix<Real> F =
        Math::SpatialMatrix<Real>::Identity(2, 2) + U * cellG.gradN;
      const Math::SpatialMatrix<Real> gradP = P * cellG.gradN;
      const Real sigDetA =
        static_cast<Real>(cellG.sigmaK) * cellG.detAK;
      const Real jU      = sigDetA * F.determinant() / cellG.Jscale;
      const Real jAfter  =
        sigDetA * (F + gradP).determinant() / cellG.Jscale;
      const Math::SpatialMatrix<Real> Finv = F.inverse();
      const Real trTerm  = (Finv * gradP).trace();
      const Real Dj      = jU * trTerm;
      const Real margin  = jU + Dj - jMinRatio;

      if (jU      < d.minJRatioBefore)    d.minJRatioBefore    = jU;
      if (jAfter  < d.minJRatioAfterFull) d.minJRatioAfterFull = jAfter;
      if (jAfter <= jMinRatio)            ++d.inadmAfter;
      if (margin  < d.minPredMargin) d.minPredMargin = margin;
      if (Dj < -tau * (jU - jMinRatio))   ++d.predInadm;
    }
    return d;
  }
}

#endif // LEVELSETLSRWAVYCIRCLESWEEP_DIAGNOSTICS_H
