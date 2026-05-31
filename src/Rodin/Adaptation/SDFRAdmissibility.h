/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_SDFRADMISSIBILITY_H
#define RODIN_ADAPTATION_SDFRADMISSIBILITY_H

#include <cstddef>
#include <limits>
#include <vector>

#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Types.h"

#include "CellGeomCache.h"

namespace Rodin::Adaptation
{
  struct SDFRAdmissibilityReport
  {
    // Dimensionless branch-normalized Jacobian ratio:
    //   j_K^u = sigma_K det(A_K^u) / J_K_scale.
    Real minJRatio = std::numeric_limits<Real>::infinity();
    std::size_t inadmissibleCount = 0;
  };

  struct SDFRLineSearchResult
  {
    Real alphaAccepted = Real(0);
    std::size_t backtracks = 0;
    Real minJRatioAtAlpha1 = std::numeric_limits<Real>::quiet_NaN();
    std::size_t inadmissibleCountAtAlpha1 = 0;
    Real minJRatioAccepted = std::numeric_limits<Real>::quiet_NaN();
    std::size_t inadmissibleCountAccepted = 0;
    bool succeeded = false;
  };

  template <class FES>
  SDFRAdmissibilityReport evaluateSDFRAdmissibilityP1(
      const Math::Vector<Real>& uData,
      const std::vector<CellGeomCache>& cellCache,
      const FES& fes,
      Real jMinRatio)
  {
    SDFRAdmissibilityReport rep;
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
      const Real sigDetA =
        static_cast<Real>(cellG.sigmaK) * cellG.detAK;
      const Real j = sigDetA * F.determinant() / cellG.Jscale;

      if (j < rep.minJRatio)
        rep.minJRatio = j;
      if (j <= jMinRatio)
        ++rep.inadmissibleCount;
    }
    return rep;
  }

  template <class Evaluator>
  SDFRLineSearchResult runSDFRAdmissibilityLineSearch(
      Math::Vector<Real>& u,
      const Math::Vector<Real>& p,
      Evaluator&& evaluator,
      Real jLineSearchRatio,
      Real alphaInit = Real(1),
      Real reduction = Real(0.5),
      Real alphaMin = Real(1e-6))
  {
    SDFRLineSearchResult result;
    const Math::Vector<Real> uOld = u;

    Real alpha = alphaInit;
    bool firstIter = true;
    SDFRAdmissibilityReport lastAdm;

    while (alpha >= alphaMin)
    {
      const Math::Vector<Real> uTrial = uOld + alpha * p;
      const SDFRAdmissibilityReport adm = evaluator(uTrial);

      if (firstIter)
      {
        result.minJRatioAtAlpha1 = adm.minJRatio;
        result.inadmissibleCountAtAlpha1 = adm.inadmissibleCount;
        firstIter = false;
      }
      lastAdm = adm;

      if (adm.minJRatio > jLineSearchRatio
          && adm.inadmissibleCount == 0)
      {
        u = uTrial;
        result.alphaAccepted = alpha;
        result.minJRatioAccepted = adm.minJRatio;
        result.inadmissibleCountAccepted = adm.inadmissibleCount;
        result.succeeded = true;
        return result;
      }

      alpha *= reduction;
      ++result.backtracks;
    }

    u = uOld;
    result.minJRatioAccepted = lastAdm.minJRatio;
    result.inadmissibleCountAccepted = lastAdm.inadmissibleCount;
    return result;
  }
}

#endif
