/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_LSRADMISSIBILITY_H
#define RODIN_ADAPTATION_LSRADMISSIBILITY_H

#include <cstddef>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Types.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/Variational/IntegrationPoint.h"
#include "Rodin/Variational/Jacobian.h"

#include "CellGeomCache.h"
#include "LSRIntegrators.h"

namespace Rodin::Adaptation
{
  struct LSRAdmissibilityReport
  {
    // Dimensionless branch-normalized Jacobian ratio:
    //   j_K^u = sigma_K det(A_K^u) / J_K_scale.
    Real minJRatio = std::numeric_limits<Real>::infinity();
    std::size_t inadmissibleCount = 0;
    // Per-cell RELATIVE distortion quality, max over cells:
    //   Q_rel(F) = ||F||_F^2 / (d * (det F)^(2/d)),
    //   F = I + grad u = A_K^u (A_K)^{-1}.
    // Identity-neutral: Q_rel = 1 at u = 0 for every cell regardless of
    // background shape. Larger values indicate deformation introduced by
    // u (anisotropy / stretching relative to the background cell). This
    // is what `LSRParameters::qRelMax` caps in the line search.
    Real maxQRel = Real(0);
  };

  struct LSRLineSearchResult
  {
    Real alphaAccepted = Real(0);
    std::size_t backtracks = 0;
    Real minJRatioAccepted = std::numeric_limits<Real>::quiet_NaN();
    std::size_t inadmissibleCountAccepted = 0;
    Real maxQRelAccepted = std::numeric_limits<Real>::quiet_NaN();
    bool succeeded = false;
  };

  struct LSRAlphaPredictorReport
  {
    Real alphaMax = std::numeric_limits<Real>::infinity();
    std::size_t limitingPointCount = 0;
  };

  namespace detail
  {
    inline Real firstPositiveQuadraticRoot(
        Real a, Real b, Real c)
    {
      constexpr Real eps = Real(1e-14);
      Real root = std::numeric_limits<Real>::infinity();

      if (std::abs(a) <= eps)
      {
        if (b < Real(0))
          root = -c / b;
        return root > Real(0) ? root : std::numeric_limits<Real>::infinity();
      }

      const Real disc = b * b - Real(4) * a * c;
      if (disc < Real(0))
        return root;

      const Real s = std::sqrt(std::max(disc, Real(0)));
      const Real sign = b >= Real(0) ? Real(1) : Real(-1);
      const Real q = -Real(0.5) * (b + sign * s);
      Real r1, r2;
      if (q == Real(0))
      {
        r1 = (-b - s) / (Real(2) * a);
        r2 = (-b + s) / (Real(2) * a);
      }
      else
      {
        r1 = q / a;
        r2 = c / q;
      }

      if (r1 > Real(0))
        root = std::min(root, r1);
      if (r2 > Real(0))
        root = std::min(root, r2);
      return root;
    }
  }

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

      // Q_rel(F) = ||F||_F^2 / (d * (det F)^(2/d))
      // (2D affine prototype, d=2, so the exponent 2/d = 1).
      // For a positively-oriented admissible cell det F > 0 (and the
      // dimensionless ratio j = sigma_K det(F) det(A_K) / J_scale > 0);
      // if the cell is inadmissible we skip Q_rel to avoid a negative
      // denominator dominating the report.
      if (detF > Real(0))
      {
        const Real q = F.squaredNorm() / (Real(2) * std::pow(detF, dExp));
        if (q > rep.maxQRel) rep.maxQRel = q;
      }
    }
    return rep;
  }

  template <class Displacement>
  LSRAdmissibilityReport evaluateLSRAdmissibilitySampled(
      Displacement& u,
      const Math::Vector<Real>& uData,
      Real jMinRatio,
      std::size_t quadratureOrder = 0)
  {
    using Variational::IntegrationPoint;
    using Variational::Jacobian;

    auto& uMutable = u;
    uMutable.getData() = uData;

    LSRAdmissibilityReport rep;
    const auto& fes = uMutable.getFiniteElementSpace();
    const auto& mesh = fes.getMesh();
    const std::size_t dim = mesh.getDimension();
    const std::size_t vdim = fes.getVectorDimension();
    if (dim != vdim)
      throw std::runtime_error(
          "evaluateLSRAdmissibilitySampled: displacement dimension "
          "must equal mesh dimension.");

    auto gradU = Jacobian(uMutable);
    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const auto& cell = *cellIt;
      const auto& fe =
        fes.getFiniteElement(cell.getDimension(), cell.getIndex());
      const auto& qf =
        QF::PolytopeQuadratureFormula::get(
            quadratureOrder > 0 ? quadratureOrder : lsrQuadOrderFor(fe.getOrder()),
            cell.getGeometry());
      const auto& quadrature = cell.getQuadrature(qf);
      for (std::size_t q = 0; q < quadrature.getSize(); ++q)
      {
        const auto& pt = quadrature.getPoint(q);
        const IntegrationPoint ip(pt, &qf, q);
        Math::SpatialMatrix<Real> F =
          Math::SpatialMatrix<Real>::Identity(dim, dim)
          + gradU.getValue(ip);
        const Real j = F.determinant();

        if (j < rep.minJRatio)
          rep.minJRatio = j;
        if (j <= jMinRatio)
          ++rep.inadmissibleCount;

        // Q_rel(F) is defined directly from F = I + grad u — no need for
        // the cell Jacobian. Cap on the relative-distortion quality of
        // the moved cell wrt its background shape.
        if (j > Real(0))
        {
          const Real qRel =
            F.squaredNorm()
            / (static_cast<Real>(dim)
               * std::pow(j, Real(2) / static_cast<Real>(dim)));
          if (qRel > rep.maxQRel)
            rep.maxQRel = qRel;
        }
      }
    }
    return rep;
  }

  template <class Displacement>
  LSRAlphaPredictorReport predictSampledQuadraticAlpha(
      Displacement& u,
      const Math::Vector<Real>& uData,
      const Math::Vector<Real>& direction,
      Real jLineSearchRatio,
      Real alphaSafety,
      std::size_t quadratureOrder = 0)
  {
    using Variational::IntegrationPoint;
    using Variational::Jacobian;

    auto& uMutable = u;
    const auto& fes = uMutable.getFiniteElementSpace();
    const auto& mesh = fes.getMesh();
    const std::size_t dim = mesh.getDimension();
    if (dim != 2)
      return {};

    const std::size_t vdim = fes.getVectorDimension();
    if (dim != vdim)
      throw std::runtime_error(
          "predictSampledQuadraticAlpha: displacement dimension "
          "must equal mesh dimension.");

    struct Sample
    {
      Math::SpatialMatrix<Real> F;
    };
    std::vector<Sample> samples;
    samples.reserve(mesh.getCellCount() * 16);

    uMutable.getData() = uData;
    auto gradU = Jacobian(uMutable);
    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const auto& cell = *cellIt;
      const auto& fe =
        fes.getFiniteElement(cell.getDimension(), cell.getIndex());
      const auto& qf =
        QF::PolytopeQuadratureFormula::get(
            quadratureOrder > 0 ? quadratureOrder : lsrQuadOrderFor(fe.getOrder()),
            cell.getGeometry());
      const auto& quadrature = cell.getQuadrature(qf);
      for (std::size_t q = 0; q < quadrature.getSize(); ++q)
      {
        const auto& pt = quadrature.getPoint(q);
        const IntegrationPoint ip(pt, &qf, q);
        samples.push_back({
            Math::SpatialMatrix<Real>::Identity(dim, dim)
              + gradU.getValue(ip)});
      }
    }

    uMutable.getData() = direction;
    auto gradP = Jacobian(uMutable);
    LSRAlphaPredictorReport rep;
    std::size_t sampleIndex = 0;
    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const auto& cell = *cellIt;
      const auto& fe =
        fes.getFiniteElement(cell.getDimension(), cell.getIndex());
      const auto& qf =
        QF::PolytopeQuadratureFormula::get(
            quadratureOrder > 0 ? quadratureOrder : lsrQuadOrderFor(fe.getOrder()),
            cell.getGeometry());
      const auto& quadrature = cell.getQuadrature(qf);
      for (std::size_t q = 0; q < quadrature.getSize(); ++q, ++sampleIndex)
      {
        const auto& pt = quadrature.getPoint(q);
        const IntegrationPoint ip(pt, &qf, q);
        const auto& F = samples[sampleIndex].F;
        const Math::SpatialMatrix<Real> P = gradP.getValue(ip);

        const Real c = F.determinant() - jLineSearchRatio;
        if (c <= Real(0))
        {
          rep.alphaMax = Real(0);
          ++rep.limitingPointCount;
          continue;
        }

        const Real g1 = (F + P).determinant() - jLineSearchRatio;
        const Real gm1 = (F - P).determinant() - jLineSearchRatio;
        const Real a = Real(0.5) * (g1 + gm1 - Real(2) * c);
        const Real b = Real(0.5) * (g1 - gm1);
        const Real root = detail::firstPositiveQuadraticRoot(a, b, c);
        if (std::isfinite(root))
        {
          const Real safeRoot = alphaSafety * root;
          if (safeRoot < rep.alphaMax)
          {
            rep.alphaMax = std::max<Real>(Real(0), safeRoot);
            rep.limitingPointCount = 1;
          }
          else if (safeRoot == rep.alphaMax)
          {
            ++rep.limitingPointCount;
          }
        }
      }
    }

    uMutable.getData() = uData;
    return rep;
  }
}

#endif
