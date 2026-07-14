/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_CELLGEOMCACHE_H
#define RODIN_ADAPTATION_CELLGEOMCACHE_H

/**
 * @file
 * @brief Per-cell geometry cache for the 2D triangular affine prototype of
 *        the Jacobian admissibility barrier and the intrinsic shape
 *        quality energy.
 *
 * Scope.
 *   This header is specialised to PLANAR 2D MESHES with AFFINE TRIANGLE
 *   cells. It is NOT FES-independent: it assumes a P1 nodal basis with
 *   three vertex dofs per cell and constant gradient on the cell.
 *   Extending it to general FES / high-order polytopes requires replacing
 *   the per-cell A_K = const, gradN = const cache by a quadrature-driven
 *   evaluation of A_K^u = D(F_K + u_h o F_K) at each integration point.
 */

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include <unordered_map>

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Geometry/PolytopeTransformation.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/Types.h"

namespace Rodin::Adaptation
{
  /**
   * @brief Per-cell geometry cache for the 2D triangular affine prototype.
   *
   * For each triangle cell K the cache stores:
   *   - sigma_K: the verified-consistent sign of det A_K(xhat_q) across
   *              the cell's reference quadrature.
   *   - Jscale:  J_K_scale = (1 / |Khat|) * integral_Khat |det A_K| dxhat.
   *   - A, Ainv, AinvT, C = A A^T: the affine map and its associates.
   *   - gradN: 3x2 spatial gradients of the P1 nodal basis (constant on K).
   *
   * The matrix A is sampled from Rodin's PolytopeTransformation at the
   * reference centroid; for an affine triangle A is constant in xhat and
   * one sample is exact.
   */
  struct CellGeomCache
  {
    /// @brief Mesh cell index for this cache entry.
      Index index = 0;
      Real area = 0;                    ///< |K|
      Real detAK = 0;                   ///< signed det(A_K), constant on K
      Real Jscale = 0;                  ///< (1/|Khat|) integral |det A_K| dxhat
      int sigmaK = 1;                  ///< verified sign(det A_K(xhat_q))
      Math::SpatialMatrix<Real> A;      ///< A_K = D F_K
      Math::SpatialMatrix<Real> Ainv;   ///< A_K^{-1}
      Math::SpatialMatrix<Real> AinvT;  ///< A_K^{-T}
      Math::SpatialMatrix<Real> C;      ///< A_K A_K^T
      Math::Matrix<Real> gradN;         ///< 3x2: rows are spatial grads of P1 basis
    /// @brief Vertex indices of the triangular cell.
      std::array<Index, 3> vertices = {{0, 0, 0}};
  };

  /**
   * @brief Build a `CellGeomCache` per cell, in mesh-iteration order.
   *
   * The returned vector is indexed by LOCAL iteration order — not by mesh
   * cell index. The cache entry `.index` records the parent mesh's cell
   * index; the second returned map gives the parent-index -> local-index
   * translation for use by the rest of the pipeline.
   *
   * Throws if A_K(xhat_q) is singular at any quadrature point, or if the
   * branch sign sigma_K(xhat_q) is inconsistent across quadrature points
   * (which would mean a curved or inverted element; this 2D affine
   * prototype rejects them).
   */
  inline std::pair<std::vector<CellGeomCache>, std::unordered_map<Index, size_t>>
  precomputeCellGeometry(const Geometry::Mesh<Rodin::Context::Local>& mesh)
  {
    std::vector<CellGeomCache> cache;
    cache.reserve(mesh.getCellCount());
    std::unordered_map<Index, size_t> cellToLocal;
    cellToLocal.reserve(mesh.getCellCount());

    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const auto& cell = *cellIt;
      const auto& vertices = cell.getVertices();
      CellGeomCache c;
      c.index = cell.getIndex();
      for (size_t i = 0; i < 3; ++i)
        c.vertices[i] = vertices[i];

      const auto& qf = QF::PolytopeQuadratureFormula::get(2, cell.getGeometry());
      const auto& quadrature = cell.getQuadrature(qf);
      const size_t nqp = quadrature.getSize();
      if (nqp == 0)
        throw std::runtime_error("precomputeCellGeometry: cell quadrature is empty.");

      Real khatMeasure = 0;
      Real integralAbsDetA = 0;
      Math::SpatialMatrix<Real> Aq;
      int sigmaFirst = 0;
      for (size_t q = 0; q < nqp; ++q)
      {
        const auto& pt = quadrature.getPoint(q);
        const auto& rc = pt.getReferenceCoordinates();
        const Real wq = qf.getWeight(q);
        khatMeasure += wq;
        cell.getTransformation().jacobian(Aq, rc);
        const Real detAq = Aq.determinant();
        if (detAq == Real(0))
          throw std::runtime_error(
            "precomputeCellGeometry: degenerate A_K(xhat_q) at cell " +
            std::to_string(c.index));
        const int sigmaq = detAq > 0 ? 1 : -1;
        if (q == 0)
          sigmaFirst = sigmaq;
        else if (sigmaq != sigmaFirst)
          throw std::runtime_error("precomputeCellGeometry: sigma_Kq inconsistent across "
                                   "quadrature points at cell " +
            std::to_string(c.index) +
            " (the 2D affine prototype rejects curved/inverted cells).");
        integralAbsDetA += wq * std::abs(detAq);
      }
      c.sigmaK = sigmaFirst;
      c.Jscale = integralAbsDetA / khatMeasure;

      Math::SpatialPoint rcCentroid(2);
      rcCentroid(0) = Real(1) / Real(3);
      rcCentroid(1) = Real(1) / Real(3);
      cell.getTransformation().jacobian(c.A, rcCentroid);
      c.detAK = c.A.determinant();
      c.area = Real(0.5) * std::abs(c.detAK);
      c.Ainv = c.A.inverse();
      c.AinvT = c.Ainv.transpose();
      c.C = c.A * c.A.transpose();

      c.gradN.resize(3, 2);
      const auto col0 = c.AinvT.col(0);
      const auto col1 = c.AinvT.col(1);
      c.gradN(0, 0) = -(col0(0) + col1(0));
      c.gradN(0, 1) = -(col0(1) + col1(1));
      c.gradN(1, 0) = col0(0);
      c.gradN(1, 1) = col0(1);
      c.gradN(2, 0) = col1(0);
      c.gradN(2, 1) = col1(1);

      cellToLocal[c.index] = cache.size();
      cache.push_back(std::move(c));
    }
    return {std::move(cache), std::move(cellToLocal)};
  }
}

#endif
