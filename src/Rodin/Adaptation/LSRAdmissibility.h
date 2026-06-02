/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_LSRADMISSIBILITY_H
#define RODIN_ADAPTATION_LSRADMISSIBILITY_H

#include <cstddef>
#include <limits>
#include <vector>

#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Types.h"

#include "CellGeomCache.h"

namespace Rodin::Adaptation
{
  struct LSRAdmissibilityReport
  {
    // Dimensionless branch-normalized Jacobian ratio:
    //   j_K^u = sigma_K det(A_K^u) / J_K_scale.
    Real minJRatio = std::numeric_limits<Real>::infinity();
    std::size_t inadmissibleCount = 0;
    // Per-cell intrinsic shape quality, max over cells:
    //   Q_K = ||A_K^u||_F^2 / (d * (sigma_K det A_K^u)^(2/d)).
    // Equilateral cells give Q = 1; larger values indicate anisotropy
    // (a sliver has Q >> 1). Bounding maxQShape via the line search is
    // what keeps the mesh clean on under-resolved or partial fits.
    Real maxQShape = Real(0);
  };

  struct LSRLineSearchResult
  {
    Real alphaAccepted = Real(0);
    std::size_t backtracks = 0;
    Real minJRatioAtAlpha1 = std::numeric_limits<Real>::quiet_NaN();
    std::size_t inadmissibleCountAtAlpha1 = 0;
    Real maxQShapeAtAlpha1 = std::numeric_limits<Real>::quiet_NaN();
    Real minJRatioAccepted = std::numeric_limits<Real>::quiet_NaN();
    std::size_t inadmissibleCountAccepted = 0;
    Real maxQShapeAccepted = std::numeric_limits<Real>::quiet_NaN();
    bool succeeded = false;
  };

  template <class FES>
  LSRAdmissibilityReport evaluateLSRAdmissibilityP1(
      const Math::Vector<Real>& uData,
      const std::vector<CellGeomCache>& cellCache,
      const FES& fes,
      Real jMinRatio)
  {
    LSRAdmissibilityReport rep;
    const Real dExp = Real(1); // 2/d with d=2
    for (const auto& cellG : cellCache)
    {
      Math::SpatialMatrix<Real> U(2, 3);
      for (std::size_t i = 0; i < 3; ++i)
      {
        const Index vertex = cellG.vertices[i];
        const auto& dofs = fes.getDOFs(0, vertex);
        U(0, i) = uData(dofs[0]);
        U(1, i) = uData(dofs[1]);
      }

      const Math::SpatialMatrix<Real> F =
        Math::SpatialMatrix<Real>::Identity(2, 2) + U * cellG.gradN;
      const Real detF = F.determinant();
      const Real sigDetA =
        static_cast<Real>(cellG.sigmaK) * cellG.detAK;
      const Real j = sigDetA * detF / cellG.Jscale;

      if (j < rep.minJRatio)
        rep.minJRatio = j;
      if (j <= jMinRatio)
        ++rep.inadmissibleCount;

      // Q_shape = ||F A_K||_F^2 / (2 * (sigma_K det A_K det F)^(1))
      // (2D affine prototype, d=2, so the exponent 2/d = 1).
      // For a positively-oriented admissible cell the denominator is
      // > 0; if the cell is inadmissible (j <= 0) we skip Q to avoid a
      // negative denominator dominating the report.
      if (j > Real(0))
      {
        const Math::SpatialMatrix<Real> FA = F * cellG.A;
        const Real numerator = FA.squaredNorm();
        const Real denom =
          Real(2) * std::pow(sigDetA * detF, dExp);
        const Real q = numerator / denom;
        if (q > rep.maxQShape) rep.maxQShape = q;
      }
    }
    return rep;
  }

  template <class Evaluator>
  LSRLineSearchResult runLSRAdmissibilityLineSearch(
      Math::Vector<Real>& u,
      const Math::Vector<Real>& p,
      Evaluator&& evaluator,
      Real jLineSearchRatio,
      Real alphaInit = Real(1),
      Real reduction = Real(0.5),
      Real alphaMin = Real(1e-6),
      Real qShapeMax = std::numeric_limits<Real>::infinity())
  {
    LSRLineSearchResult result;
    const Math::Vector<Real> uOld = u;

    Real alpha = alphaInit;
    bool firstIter = true;
    LSRAdmissibilityReport lastAdm;

    while (alpha >= alphaMin)
    {
      const Math::Vector<Real> uTrial = uOld + alpha * p;
      const LSRAdmissibilityReport adm = evaluator(uTrial);

      if (firstIter)
      {
        result.minJRatioAtAlpha1 = adm.minJRatio;
        result.inadmissibleCountAtAlpha1 = adm.inadmissibleCount;
        result.maxQShapeAtAlpha1 = adm.maxQShape;
        firstIter = false;
      }
      lastAdm = adm;

      const bool jOK =
           adm.minJRatio > jLineSearchRatio
        && adm.inadmissibleCount == 0;
      const bool qOK = adm.maxQShape <= qShapeMax;

      if (jOK && qOK)
      {
        u = uTrial;
        result.alphaAccepted = alpha;
        result.minJRatioAccepted = adm.minJRatio;
        result.inadmissibleCountAccepted = adm.inadmissibleCount;
        result.maxQShapeAccepted = adm.maxQShape;
        result.succeeded = true;
        return result;
      }

      alpha *= reduction;
      ++result.backtracks;
    }

    u = uOld;
    result.minJRatioAccepted = lastAdm.minJRatio;
    result.inadmissibleCountAccepted = lastAdm.inadmissibleCount;
    result.maxQShapeAccepted = lastAdm.maxQShape;
    return result;
  }
}

#endif
