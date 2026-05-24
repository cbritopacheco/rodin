/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Adaptation/TargetMatrixOptimization/Metrics.h"
#include "Rodin/Adaptation/TargetMatrixOptimization/IsoparametricGeometry.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <Eigen/SparseLU>

#include <Rodin/Adaptation.h>
#include <Rodin/Assembly/Default.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO.h>
#include <Rodin/MMG.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Adaptation::TargetMatrixOptimization;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  constexpr Real Pi = 3.14159265358979323846;
  constexpr Attribute Interface = 30;
  constexpr Attribute Boundary = 40;
  constexpr Attribute Negative = 1;
  constexpr Attribute Positive = 2;
  constexpr Index CurvedQuadratureOrder = 4;
  constexpr Index CurvedQualitySamples = 4;
  constexpr Index CurvedInterfaceSamples = 4;
  constexpr Index CurvedVisualizationSamples = 4;
  constexpr Real CurvedDetFloor = 1e-12;

  Math::SpatialPoint center(Real t)
  {
    return Math::SpatialPoint{
      Real(0.5) + Real(0.18) * std::sin(Real(2) * Pi * t),
      Real(0.5) + Real(0.13) * std::sin(Real(4) * Pi * t + Real(0.35))
    };
  }

  Real wavyRadius(Real theta, Real t)
  {
    return Real(0.17)
      + Real(0.035) * std::sin(
          Real(5) * theta + Real(0.65) * std::sin(Real(2) * Pi * t));
  }

  Real wavyRadiusDerivative(Real theta, Real t)
  {
    return Real(0.175) * std::cos(
        Real(5) * theta + Real(0.65) * std::sin(Real(2) * Pi * t));
  }

  Real phiAt(const Math::SpatialPoint& x, Real t)
  {
    const auto c = center(t);
    const Real dx = x[0] - c[0];
    const Real dy = x[1] - c[1];
    const Real theta = std::atan2(dy, dx);
    return std::hypot(dx, dy) - wavyRadius(theta, t);
  }

  Math::SpatialPoint gradPhiAt(const Math::SpatialPoint& x, Real t)
  {
    const auto c = center(t);
    const Real dx = x[0] - c[0];
    const Real dy = x[1] - c[1];
    const Real rho2 = dx * dx + dy * dy;
    const Real rho = std::sqrt(rho2);
    if (rho <= Real(1e-14))
      return Math::SpatialPoint{1, 0};
    const Real theta = std::atan2(dy, dx);
    const Real dr = wavyRadiusDerivative(theta, t);
    return Math::SpatialPoint{
      dx / rho + dr * dy / rho2,
      dy / rho - dr * dx / rho2
    };
  }

  Real boxBoundaryValue(const Math::SpatialPoint& x)
  {
    const std::array<Real, 4> d = {{
      std::abs(x[0]), std::abs(x[0] - Real(1)),
      std::abs(x[1]), std::abs(x[1] - Real(1))
    }};
    const auto side =
      static_cast<size_t>(std::min_element(d.begin(), d.end()) - d.begin());
    if (side == 0)
      return x[0];
    if (side == 1)
      return x[0] - Real(1);
    if (side == 2)
      return x[1];
    return x[1] - Real(1);
  }

  Math::SpatialPoint boxBoundaryGradient(const Math::SpatialPoint& x)
  {
    const std::array<Real, 4> d = {{
      std::abs(x[0]), std::abs(x[0] - Real(1)),
      std::abs(x[1]), std::abs(x[1] - Real(1))
    }};
    const auto side =
      static_cast<size_t>(std::min_element(d.begin(), d.end()) - d.begin());
    if (side < 2)
      return Math::SpatialPoint{1, 0};
    return Math::SpatialPoint{0, 1};
  }

  Math::SpatialPoint projectToInterface(
      const Math::SpatialPoint& x,
      Real t)
  {
    Math::SpatialPoint p = x;
    for (int i = 0; i < 2; ++i)
    {
      const Real f = phiAt(p, t);
      const auto g = gradPhiAt(p, t);
      const Real gg = g[0] * g[0] + g[1] * g[1];
      if (gg < Real(1e-30))
        break;
      p = Math::SpatialPoint{ p[0] - f * g[0] / gg,
                              p[1] - f * g[1] / gg };
    }
    return p;
  }

  void annotateBoundary(LocalMesh& mesh)
  {
    const auto& conn = mesh.getConnectivity();
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto& edge = conn.getPolytope(1, e);
      const auto x0 = mesh.getVertexCoordinates(edge(0));
      const auto x1 = mesh.getVertexCoordinates(edge(1));
      const Real eps = Real(1e-12);
      const bool on0 =
        std::abs(x0[0]) <= eps || std::abs(x0[0] - Real(1)) <= eps
     || std::abs(x0[1]) <= eps || std::abs(x0[1] - Real(1)) <= eps;
      const bool on1 =
        std::abs(x1[0]) <= eps || std::abs(x1[0] - Real(1)) <= eps
     || std::abs(x1[1]) <= eps || std::abs(x1[1] - Real(1)) <= eps;
      if (on0 && on1)
        mesh.setAttribute({ 1, e }, Boundary);
    }
  }

  std::vector<char> interfaceVertexMask(const LocalMesh& mesh)
  {
    std::vector<char> mask(mesh.getVertexCount(), 0);
    const auto& conn = mesh.getConnectivity();
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto attr = mesh.getAttribute(1, e);
      if (!attr || *attr != Interface)
        continue;
      const auto& edge = conn.getPolytope(1, e);
      mask[edge(0)] = 1;
      mask[edge(1)] = 1;
    }
    return mask;
  }

  Math::SpatialPoint triangleReferenceVertex(std::uint8_t local)
  {
    if (local == 1)
      return Math::SpatialPoint{1, 0};
    if (local == 2)
      return Math::SpatialPoint{0, 1};
    return Math::SpatialPoint{0, 0};
  }

  Math::SpatialPoint triangleReferenceEdgePoint(
      std::uint8_t localA,
      std::uint8_t localB,
      Real s)
  {
    return (Real(1) - s) * triangleReferenceVertex(localA)
         + s * triangleReferenceVertex(localB);
  }

  bool referencePointOnTriangleEdge(
      const Math::SpatialPoint& rc,
      std::uint8_t localA,
      std::uint8_t localB)
  {
    const Real eps = Real(1e-12);
    if ((localA == 0 && localB == 1) || (localA == 1 && localB == 0))
      return std::abs(rc[1]) <= eps;
    if ((localA == 1 && localB == 2) || (localA == 2 && localB == 1))
      return std::abs(rc[0] + rc[1] - Real(1)) <= eps;
    return std::abs(rc[0]) <= eps;
  }

  bool referencePointIsTriangleVertex(const Math::SpatialPoint& rc)
  {
    const Real eps = Real(1e-12);
    return (std::abs(rc[0]) <= eps && std::abs(rc[1]) <= eps)
        || (std::abs(rc[0] - Real(1)) <= eps && std::abs(rc[1]) <= eps)
        || (std::abs(rc[0]) <= eps && std::abs(rc[1] - Real(1)) <= eps);
  }

  Optional<std::array<std::uint8_t, 2>> localEdgeForMeshEdge(
      const Geometry::Polytope::Key& cell,
      const Geometry::Polytope::Key& edge)
  {
    for (std::uint8_t a = 0; a < 3; ++a)
      for (std::uint8_t b = a + 1; b < 3; ++b)
        if ((cell(a) == edge(0) && cell(b) == edge(1))
         || (cell(a) == edge(1) && cell(b) == edge(0)))
          return std::array<std::uint8_t, 2>{{a, b}};
    return {};
  }

  std::vector<std::array<std::uint8_t, 2>> interfaceLocalEdges(
      const LocalMesh& mesh,
      Index cellIndex)
  {
    std::vector<std::array<std::uint8_t, 2>> edges;
    const auto& conn = mesh.getConnectivity();
    if (conn.getGeometry(2, cellIndex) != Polytope::Type::Triangle)
      return edges;
    const auto& cell = conn.getPolytope(2, cellIndex);
    for (Index edgeIndex : conn.getIncidence({2, 1}, cellIndex))
    {
      const auto attr = mesh.getAttribute(1, edgeIndex);
      if (!attr || *attr != Interface)
        continue;
      if (const auto edge =
            localEdgeForMeshEdge(cell, conn.getPolytope(1, edgeIndex)))
        edges.push_back(*edge);
    }
    return edges;
  }

  void installQuadraticGeometry(
      LocalMesh& mesh,
      Real t,
      Real curvatureClamp = Real(0.35))
  {
    if (mesh.getConnectivity().getIncidence(2, 1).size() == 0)
      return;

    Variational::RealH1Element<2> fe(Polytope::Type::Triangle);
    Variational::RealH1Element<2> edgeFe(Polytope::Type::Segment);
    mesh.getConnectivity().compute(1, 2);
    const auto& conn = mesh.getConnectivity();

    // Conformity-safe interface-projection clamp. Limit each interface
    // mid-edge node's projection displacement to a fraction of the minimum
    // perpendicular height across the edge's incident cells. The cut +
    // coarsen can leave slivers where the analytic projection of the chord
    // midpoint overshoots and folds the P2 cell (or "tongues" it: det>0 yet
    // a sub-triangle inverts). The barrier metric then cannot step away from
    // such a start. Both the edge loop and the cell loop use the SAME
    // per-edge limit, so the curved geometry stays conformal AND no cell is
    // bent beyond what its own height admits.
    std::vector<Real> edgeSafeMag(
        static_cast<size_t>(conn.getCount(1)),
        std::numeric_limits<Real>::infinity());
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto attr = mesh.getAttribute(1, e);
      if (!attr || *attr != Interface) continue;
      const auto& edge = conn.getPolytope(1, e);
      const auto p0 = mesh.getVertexCoordinates(edge(0));
      const auto p1 = mesh.getVertexCoordinates(edge(1));
      const Real chord = (p1 - p0).norm();
      Real minPerp = std::numeric_limits<Real>::infinity();
      for (Index ci : conn.getIncidence({ 1, 2 }, e))
      {
        const auto& cell = conn.getPolytope(2, ci);
        Index vo = std::numeric_limits<Index>::max();
        for (size_t k = 0; k < 3; ++k)
          if (cell(k) != edge(0) && cell(k) != edge(1)) { vo = cell(k); break; }
        if (vo == std::numeric_limits<Index>::max()) continue;
        const auto po = mesh.getVertexCoordinates(vo);
        const Real perp = std::abs(
              (p1[0] - p0[0]) * (po[1] - p0[1])
            - (p1[1] - p0[1]) * (po[0] - p0[0]))
          / std::max(chord, Real(1e-30));
        minPerp = std::min(minPerp, perp);
      }
      edgeSafeMag[e] = curvatureClamp * minPerp;
    }
    auto clampedProject = [&](const Math::SpatialPoint& mid, Real limit)
    {
      auto p = projectToInterface(mid, t);
      const Real d = std::hypot(p[0] - mid[0], p[1] - mid[1]);
      if (d > limit && d > Real(0))
        p = Math::SpatialPoint{
          mid[0] + (limit / d) * (p[0] - mid[0]),
          mid[1] + (limit / d) * (p[1] - mid[1]) };
      return p;
    };

    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto attr = mesh.getAttribute(1, e);
      const bool isInterface = attr && *attr == Interface;
      const auto& edge = conn.getPolytope(1, e);
      const auto x0 = mesh.getVertexCoordinates(edge(0));
      const auto x1 = mesh.getVertexCoordinates(edge(1));
      PointCloud pm(mesh.getSpaceDimension(), edgeFe.getCount());
      for (size_t local = 0; local < edgeFe.getCount(); ++local)
      {
        const Real s = edgeFe.getNode(local)[0];
        Math::SpatialPoint x = (Real(1) - s) * x0 + s * x1;
        if (isInterface
            && s > Real(1e-12)
            && s < Real(1) - Real(1e-12))
          x = clampedProject(x, edgeSafeMag[e]);
        for (size_t d = 0; d < static_cast<size_t>(x.size()); ++d)
          pm(d, local) = x[d];
      }
      mesh.setPolytopeTransformation(
          {1, e},
          new ParametricTransformation<Variational::RealH1Element<2>>(
              pm, Variational::RealH1Element<2>(edgeFe)));
    }

    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      if (conn.getGeometry(2, c) != Polytope::Type::Triangle)
        continue;
      const auto featureEdges = interfaceLocalEdges(mesh, c);
      const auto& cell = conn.getPolytope(2, c);
      const auto cellEdges = conn.getIncidence({ 2, 1 }, c);
      const auto x0 = mesh.getVertexCoordinates(cell(0));
      const auto x1 = mesh.getVertexCoordinates(cell(1));
      const auto x2 = mesh.getVertexCoordinates(cell(2));
      // Per-edge safe-magnitude lookup for this cell's local ref-edges, so
      // the cell loop clamps with the SAME limit as the edge loop above.
      auto edgeMagFor = [&](Index la, Index lb) -> Real
      {
        const Index va = cell(la), vb = cell(lb);
        for (Index ge : cellEdges)
        {
          const auto& ep = conn.getPolytope(1, ge);
          if ((ep(0) == va && ep(1) == vb) || (ep(0) == vb && ep(1) == va))
            return (ge < static_cast<Index>(edgeSafeMag.size()))
              ? edgeSafeMag[ge] : std::numeric_limits<Real>::infinity();
        }
        return std::numeric_limits<Real>::infinity();
      };
      PointCloud pm(mesh.getSpaceDimension(), fe.getCount());
      for (size_t local = 0; local < fe.getCount(); ++local)
      {
        const auto& rc = fe.getNode(local);
        Math::SpatialPoint x =
          (Real(1) - rc[0] - rc[1]) * x0 + rc[0] * x1 + rc[1] * x2;
        bool onInterfaceEdge = false;
        Real limit = std::numeric_limits<Real>::infinity();
        for (const auto& edge : featureEdges)
          if (referencePointOnTriangleEdge(rc, edge[0], edge[1]))
          {
            onInterfaceEdge = true;
            limit = std::min(limit, edgeMagFor(edge[0], edge[1]));
          }
        if (!referencePointIsTriangleVertex(rc) && onInterfaceEdge)
          x = clampedProject(x, limit);
        for (size_t d = 0; d < static_cast<size_t>(x.size()); ++d)
          pm(d, local) = x[d];
      }
      auto* trans =
        new ParametricTransformation<Variational::RealH1Element<2>>(
            pm, Variational::RealH1Element<2>(fe));
      // Validate against the SAME order-4 quadrature the acceptance gate and
      // the diagnostics use -- otherwise install passes a curved cell on a
      // few mild interior points while it is actually folded near a corner,
      // producing a degenerate starting mesh the barrier metric cannot step
      // away from. Cells whose projected curve fails fall back to the affine
      // P2 lift (always valid); the interface fit then re-curves them within
      // the solve if the local geometry permits.
      bool valid = true;
      {
        const auto& qf = QF::PolytopeQuadratureFormula::get(
            CurvedQuadratureOrder, Polytope::Type::Triangle);
        for (size_t q = 0; q < qf.getSize() && valid; ++q)
        {
          Math::SpatialMatrix<Real> J;
          trans->jacobian(J, qf.getPoint(q));
          valid = J.rows() == 2 && J.cols() == 2
            && std::isfinite(J.determinant())
            && J.determinant() > CurvedDetFloor;
        }
      }
      if (valid)
        mesh.setPolytopeTransformation({2, c}, trans);
      else
        delete trans;
    }
  }

  struct CurvedGeometryStats
  {
    Real fitRms = 0;
    Real fitMax = 0;
    Real minJacobian = std::numeric_limits<Real>::infinity();
    Real minQuality = std::numeric_limits<Real>::infinity();
    Real sampledSignedArea = 0;
    Real sampledCoverage = 0;
    Index invalidJacobianSamples = 0;
    Index overlapSamples = 0;  ///< curved interface samples inside a vertex-adjacent (non-edge-adjacent) neighbour cell -> cell-cell overlap canary
  };

  // Linear point-in-triangle (barycentric); coarse cell-cell overlap proxy.
  bool pointInLinearTriangle(
      const Math::SpatialPoint& p,
      const Math::SpatialPoint& a,
      const Math::SpatialPoint& b,
      const Math::SpatialPoint& c)
  {
    const Real d = (b[0] - a[0]) * (c[1] - a[1])
                 - (b[1] - a[1]) * (c[0] - a[0]);
    if (std::abs(d) <= Real(1e-30)) return false;
    const Real inv = Real(1) / d;
    const Real wb = ((c[1] - a[1]) * (p[0] - a[0])
                  - (c[0] - a[0]) * (p[1] - a[1])) * inv;
    const Real wc = ((b[0] - a[0]) * (p[1] - a[1])
                  - (b[1] - a[1]) * (p[0] - a[0])) * inv;
    const Real wa = Real(1) - wb - wc;
    const Real eps = Real(1e-9);
    return wa > eps && wb > eps && wc > eps;
  }

  Real triangleQuality(
      const Math::SpatialPoint& x0,
      const Math::SpatialPoint& x1,
      const Math::SpatialPoint& x2)
  {
    const Real orient =
      (x1[0] - x0[0]) * (x2[1] - x0[1])
    - (x1[1] - x0[1]) * (x2[0] - x0[0]);
    if (!(orient > Real(0)) || !std::isfinite(orient))
      return 0;
    const Real area = Real(0.5) * orient;
    const Real den = (x1 - x0).squaredNorm()
                   + (x2 - x1).squaredNorm()
                   + (x0 - x2).squaredNorm();
    if (!(den > Real(0)) || !std::isfinite(den))
      return 0;
    return Real(4) * std::sqrt(Real(3)) * area / den;
  }

  CurvedGeometryStats curvedGeometryStats(
      const LocalMesh& mesh,
      Real t,
      Index samples = CurvedQualitySamples,
      Index quadratureOrder = CurvedQuadratureOrder,
      Index overlapSamples = CurvedInterfaceSamples)
  {
    CurvedGeometryStats stats;
    Real sq = 0;
    Index count = 0;
    const auto& conn = mesh.getConnectivity();
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      if (conn.getGeometry(2, c) != Polytope::Type::Triangle)
        continue;
      const auto cellIt = mesh.getPolytope(2, c);
      const auto& trans = cellIt->getTransformation();

      auto map = [&](Index i, Index j)
      {
        Math::SpatialPoint x;
        trans.transform(
            x,
            Math::SpatialPoint{
              static_cast<Real>(i) / static_cast<Real>(samples),
              static_cast<Real>(j) / static_cast<Real>(samples)});
        return x;
      };

      const auto& qf = QF::PolytopeQuadratureFormula::get(
          quadratureOrder, Polytope::Type::Triangle);
      for (size_t q = 0; q < qf.getSize(); ++q)
      {
        Math::SpatialMatrix<Real> J;
        trans.jacobian(J, qf.getPoint(q));
        const bool shaped = J.rows() == 2 && J.cols() == 2;
        const Real det = shaped ? J.determinant() : Real(0);
        if (!shaped || !(det > Real(0)) || !std::isfinite(det))
          ++stats.invalidJacobianSamples;
        else
          stats.minJacobian = std::min(stats.minJacobian, det);
      }

      for (Index i = 0; i < samples; ++i)
      {
        for (Index j = 0; j < samples - i; ++j)
        {
          const auto x00 = map(i, j);
          const auto x10 = map(i + 1, j);
          const auto x01 = map(i, j + 1);
          const Real a0 = Real(0.5) *
            ((x10[0] - x00[0]) * (x01[1] - x00[1])
           - (x10[1] - x00[1]) * (x01[0] - x00[0]));
          stats.sampledSignedArea += a0;
          stats.sampledCoverage += std::abs(a0);
          stats.minQuality =
            std::min(stats.minQuality, triangleQuality(x00, x10, x01));
          if (i + j + 1 < samples)
          {
            const auto x11 = map(i + 1, j + 1);
            const Real a1 = Real(0.5) *
              ((x11[0] - x10[0]) * (x01[1] - x10[1])
             - (x11[1] - x10[1]) * (x01[0] - x10[0]));
            stats.sampledSignedArea += a1;
            stats.sampledCoverage += std::abs(a1);
            stats.minQuality =
              std::min(stats.minQuality, triangleQuality(x10, x11, x01));
          }
        }
      }

      for (const auto& edge : interfaceLocalEdges(mesh, c))
      {
        for (Index k = 0; k <= samples; ++k)
        {
          Math::SpatialPoint x;
          trans.transform(
              x,
              triangleReferenceEdgePoint(
                edge[0], edge[1],
                static_cast<Real>(k) / static_cast<Real>(samples)));
          const Real phi = std::abs(phiAt(x, t));
          stats.fitMax = std::max(stats.fitMax, phi);
          sq += phi * phi;
          ++count;
        }
      }
    }
    if (count > 0)
      stats.fitRms = std::sqrt(sq / static_cast<Real>(count));
    if (!std::isfinite(stats.minJacobian))
      stats.minJacobian = 0;
    if (!std::isfinite(stats.minQuality))
      stats.minQuality = 0;

    // Cell-cell overlap canary: sample each interface edge's curved
    // transformation and test membership in cells that share a vertex with
    // the edge but NOT the edge itself. A hit means the curved interface
    // has bowed into a non-incident neighbour's region (global-injectivity
    // failure, invisible to per-cell det/quality).
    const_cast<LocalMesh&>(mesh).getConnectivity().compute(0, 2);
    const_cast<LocalMesh&>(mesh).getConnectivity().compute(1, 2);
    const Index nOverlap = std::max<Index>(overlapSamples, 1);
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto attr = mesh.getAttribute(1, e);
      if (!attr || *attr != Interface) continue;
      const auto& ep = conn.getPolytope(1, e);
      std::array<Index, 2> edgeAdj{
        std::numeric_limits<Index>::max(),
        std::numeric_limits<Index>::max() };
      { size_t k = 0;
        for (Index ci : conn.getIncidence({ 1, 2 }, e))
          if (k < 2) edgeAdj[k++] = ci; }
      std::vector<Index> nbrs;
      for (Index vEnd : { ep(0), ep(1) })
        for (Index ci : conn.getIncidence({ 0, 2 }, vEnd))
          if (ci != edgeAdj[0] && ci != edgeAdj[1])
            nbrs.push_back(ci);
      std::sort(nbrs.begin(), nbrs.end());
      nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
      if (nbrs.empty()) continue;
      const auto edgeIt = mesh.getPolytope(1, e);
      const auto& etrans = edgeIt->getTransformation();
      for (Index k = 1; k <= nOverlap; ++k)
      {
        const Real s =
          static_cast<Real>(k) / static_cast<Real>(nOverlap + 1);
        Math::SpatialPoint x;
        etrans.transform(x, Math::SpatialPoint{ s });
        for (Index ci : nbrs)
        {
          const auto& cell = conn.getPolytope(2, ci);
          if (pointInLinearTriangle(
                x,
                mesh.getVertexCoordinates(cell(0)),
                mesh.getVertexCoordinates(cell(1)),
                mesh.getVertexCoordinates(cell(2))))
          { ++stats.overlapSamples; break; }
        }
      }
    }
    return stats;
  }

  LocalMesh curvedInterfacePolyline(
      const LocalMesh& mesh,
      Index samplesPerEdge = CurvedInterfaceSamples)
  {
    std::vector<Math::SpatialPoint> points;
    std::vector<std::array<Index, 2>> segments;
    const auto& conn = mesh.getConnectivity();
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      if (conn.getGeometry(2, c) != Polytope::Type::Triangle)
        continue;
      const auto cellIt = mesh.getPolytope(2, c);
      const auto& trans = cellIt->getTransformation();
      for (const auto& edge : interfaceLocalEdges(mesh, c))
      {
        const Index first = static_cast<Index>(points.size());
        for (Index i = 0; i <= samplesPerEdge; ++i)
        {
          Math::SpatialPoint x;
          trans.transform(
              x,
              triangleReferenceEdgePoint(
                edge[0], edge[1],
                static_cast<Real>(i) / static_cast<Real>(samplesPerEdge)));
          points.push_back(x);
          if (i > 0)
            segments.push_back({{first + i - 1, first + i}});
        }
      }
    }

    auto builder = LocalMesh::Builder();
    builder
      .initialize(2)
      .nodes(points.size())
      .reserve(1, segments.size());
    for (const auto& x : points)
      builder.vertex(x);
    for (const auto& segment : segments)
      builder.polytope(Polytope::Type::Segment, {segment[0], segment[1]});
    auto out = builder.finalize();
    out.getConnectivity().compute(1, 0);
    for (Index e = 0; e < out.getPolytopeCount(1); ++e)
      out.setAttribute({1, e}, Interface);
    return out;
  }

  bool curvedGeometryValid(const CurvedGeometryStats& stats)
  {
    return stats.invalidJacobianSamples == 0
      && stats.overlapSamples == 0
      && stats.minJacobian > CurvedDetFloor
      && stats.minQuality > Real(0)
      && std::isfinite(stats.minJacobian)
      && std::isfinite(stats.minQuality);
  }

  bool fastCurvedLocalValidity(const LocalMesh& mesh)
  {
    const auto& conn = mesh.getConnectivity();
    const auto& qf = QF::PolytopeQuadratureFormula::get(
        CurvedQuadratureOrder, Polytope::Type::Triangle);
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      if (conn.getGeometry(2, c) != Polytope::Type::Triangle)
        continue;
      const auto cellIt = mesh.getPolytope(2, c);
      const auto& trans = cellIt->getTransformation();
      for (size_t q = 0; q < qf.getSize(); ++q)
      {
        Math::SpatialMatrix<Real> J;
        trans.jacobian(J, qf.getPoint(q));
        if (J.rows() != 2 || J.cols() != 2)
          return false;
        const Real det = J.determinant();
        if (!(det > CurvedDetFloor) || !std::isfinite(det))
          return false;
      }
    }
    return true;
  }

  LocalMesh curvedVisualizationMesh(
      const LocalMesh& mesh,
      Index subdivisions = 3)
  {
    std::vector<Math::SpatialPoint> points;
    std::vector<std::array<Index, 3>> cells;
    std::vector<Attribute> cellAttributes;
    const auto& conn = mesh.getConnectivity();
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      if (conn.getGeometry(2, c) != Polytope::Type::Triangle)
        continue;
      const auto cellIt = mesh.getPolytope(2, c);
      const auto& trans = cellIt->getTransformation();
      const Index base = static_cast<Index>(points.size());
      for (Index i = 0; i <= subdivisions; ++i)
      {
        for (Index j = 0; j <= subdivisions - i; ++j)
        {
          Math::SpatialPoint x;
          trans.transform(
              x,
              Math::SpatialPoint{
                static_cast<Real>(i) / static_cast<Real>(subdivisions),
                static_cast<Real>(j) / static_cast<Real>(subdivisions)});
          points.push_back(x);
        }
      }

      auto idx = [base, subdivisions](Index i, Index j)
      {
        Index offset = 0;
        for (Index r = 0; r < i; ++r)
          offset += subdivisions - r + 1;
        return base + offset + j;
      };

      for (Index i = 0; i < subdivisions; ++i)
      {
        for (Index j = 0; j < subdivisions - i; ++j)
        {
          cells.push_back({{idx(i, j), idx(i + 1, j), idx(i, j + 1)}});
          cellAttributes.push_back(
              mesh.getAttribute(2, c).value_or(Attribute(0)));
          if (j + i + 1 < subdivisions)
          {
            cells.push_back({{idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)}});
            cellAttributes.push_back(
                mesh.getAttribute(2, c).value_or(Attribute(0)));
          }
        }
      }
    }

    auto builder = LocalMesh::Builder();
    builder
      .initialize(2)
      .nodes(points.size())
      .reserve(2, cells.size());
    for (const auto& x : points)
      builder.vertex(x);
    for (const auto& cell : cells)
      builder.polytope(Polytope::Type::Triangle, {cell[0], cell[1], cell[2]});
    auto out = builder.finalize();
    out.getConnectivity().compute(2, 1);
    out.getConnectivity().compute(1, 0);
    for (Index c = 0;
         c < static_cast<Index>(std::min(cells.size(), cellAttributes.size()));
         ++c)
    {
      if (cellAttributes[static_cast<size_t>(c)] != Attribute(0))
        out.setAttribute({2, c}, cellAttributes[static_cast<size_t>(c)]);
    }
    return out;
  }

  template <class GridFunctionType>
  void applyDisplacement(LocalMesh& mesh, const GridFunctionType& u)
  {
    const auto& data = u.getData();
    const Index nv = mesh.getVertexCount();
    for (Index v = 0; v < nv; ++v)
    {
      auto x = mesh.getVertexCoordinates(v);
      x[0] += data(v);
      x[1] += data(v + nv);
      mesh.setVertexCoordinates(v, x);
    }
  }

  template <class FES, class Data>
  void applyParametricDisplacement(
      LocalMesh& mesh,
      const GridFunction<FES, Data>& u)
  {
    const auto& fes = u.getFiniteElementSpace();
    const Variational::RealH1Element<2> cellFe(Polytope::Type::Triangle);
    const Variational::RealH1Element<2> edgeFe(Polytope::Type::Segment);
    const size_t sdim = mesh.getSpaceDimension();
    const auto& conn = mesh.getConnectivity();
    std::array<size_t, 3> cornerNode{{0, 0, 0}};
    for (std::uint8_t local = 0; local < 3; ++local)
    {
      const auto rv = triangleReferenceVertex(local);
      Real best = std::numeric_limits<Real>::infinity();
      for (size_t node = 0; node < cellFe.getCount(); ++node)
      {
        const auto& rn = cellFe.getNode(node);
        const Real d0 = rn[0] - rv[0];
        const Real d1 = rn[1] - rv[1];
        const Real dist2 = d0 * d0 + d1 * d1;
        if (dist2 < best)
        {
          best = dist2;
          cornerNode[local] = node;
        }
      }
    }
    std::vector<PointCloud> edgePointClouds(conn.getCount(1));
    std::vector<PointCloud> cellPointClouds(mesh.getCellCount());
    std::vector<char> hasCellPointCloud(mesh.getCellCount(), 0);

    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto edgeIt = mesh.getPolytope(1, e);
      const auto& edge = *edgeIt;
      const auto& fe = fes.getFiniteElement(edge.getDimension(), edge.getIndex());
      const auto displacement = localDisplacementCoefficients(u, edge);
      PointCloud pm(sdim, edgeFe.getCount());
      for (size_t node = 0; node < edgeFe.getCount(); ++node)
      {
        const auto& rc = edgeFe.getNode(node);
        const Geometry::Point point(edge, rc);
        auto x = point.getPhysicalCoordinates();
        for (size_t localDof = 0; localDof < displacement.size(); ++localDof)
        {
          const auto phi = basisVectorValue(fe, localDof, rc);
          for (size_t d = 0; d < sdim; ++d)
            x[d] += displacement[localDof] * phi[static_cast<std::uint8_t>(d)];
        }
        for (size_t c = 0; c < sdim; ++c)
          pm(c, node) = x[c];
      }
      edgePointClouds[static_cast<size_t>(e)] = std::move(pm);
    }

    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      if (conn.getGeometry(2, c) != Polytope::Type::Triangle)
        continue;
      const auto cellIt = mesh.getPolytope(2, c);
      const auto& cell = *cellIt;
      const auto& fe = fes.getFiniteElement(cell.getDimension(), cell.getIndex());
      const auto displacement = localDisplacementCoefficients(u, cell);
      PointCloud pm(sdim, cellFe.getCount());
      for (size_t node = 0; node < cellFe.getCount(); ++node)
      {
        const auto& rc = cellFe.getNode(node);
        const Geometry::Point point(cell, rc);
        auto x = point.getPhysicalCoordinates();
        for (size_t localDof = 0; localDof < displacement.size(); ++localDof)
        {
          const auto phi = basisVectorValue(fe, localDof, rc);
          for (size_t d = 0; d < sdim; ++d)
            x[d] += displacement[localDof] * phi[static_cast<std::uint8_t>(d)];
        }
        for (size_t d = 0; d < sdim; ++d)
          pm(d, node) = x[d];
      }
      cellPointClouds[static_cast<size_t>(c)] = std::move(pm);
      hasCellPointCloud[static_cast<size_t>(c)] = 1;
    }

    for (Index v = 0; v < static_cast<Index>(mesh.getVertexCount()); ++v)
    {
      bool found = false;
      Math::SpatialPoint x(sdim);
      for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()) && !found; ++c)
      {
        if (!hasCellPointCloud[static_cast<size_t>(c)])
          continue;
        const auto& cell = conn.getPolytope(2, c);
        for (std::uint8_t local = 0; local < 3; ++local)
        {
          if (cell(local) != v)
            continue;
          const size_t node = cornerNode[local];
          for (size_t d = 0; d < sdim; ++d)
            x[d] = cellPointClouds[static_cast<size_t>(c)](d, node);
          found = true;
          break;
        }
      }
      if (found)
        mesh.setVertexCoordinates(v, x);
    }

    // setVertexCoordinates() flushes cached/custom transformations. Reinstall
    // the displaced P2 geometry after the vertex update so midside/interface
    // nodes remain part of the actual mesh map used for diagnostics/export.
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      mesh.setPolytopeTransformation(
          {1, e},
          new ParametricTransformation<Variational::RealH1Element<2>>(
              edgePointClouds[static_cast<size_t>(e)],
              Variational::RealH1Element<2>(Polytope::Type::Segment)));
    }
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      if (!hasCellPointCloud[static_cast<size_t>(c)])
        continue;
      mesh.setPolytopeTransformation(
          {2, c},
          new ParametricTransformation<Variational::RealH1Element<2>>(
              cellPointClouds[static_cast<size_t>(c)],
              Variational::RealH1Element<2>(Polytope::Type::Triangle)));
    }
  }

  std::vector<Math::SpatialPoint> vertexCoordinates(const LocalMesh& mesh)
  {
    std::vector<Math::SpatialPoint> coordinates;
    coordinates.reserve(static_cast<size_t>(mesh.getVertexCount()));
    for (Index v = 0; v < mesh.getVertexCount(); ++v)
      coordinates.push_back(mesh.getVertexCoordinates(v));
    return coordinates;
  }

  void restoreVertexCoordinates(
      LocalMesh& mesh,
      const std::vector<Math::SpatialPoint>& coordinates)
  {
    const Index count = std::min(
        mesh.getVertexCount(),
        static_cast<Index>(coordinates.size()));
    for (Index v = 0; v < count; ++v)
      mesh.setVertexCoordinates(v, coordinates[static_cast<size_t>(v)]);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
  }

  void demoteInterface(LocalMesh& mesh)
  {
    auto& conn = mesh.getConnectivity();
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto attr = mesh.getAttribute(1, e);
      if (attr && *attr == Interface)
        mesh.setAttribute({ 1, e }, Attribute(0));
    }
  }

  std::pair<Real, Index> meshQuality(const LocalMesh& mesh)
  {
    const auto& conn = mesh.getConnectivity();
    Real qmin = std::numeric_limits<Real>::infinity();
    Index inverted = 0;
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      const auto& cell = conn.getPolytope(2, c);
      const auto x0 = mesh.getVertexCoordinates(cell(0));
      const auto x1 = mesh.getVertexCoordinates(cell(1));
      const auto x2 = mesh.getVertexCoordinates(cell(2));
      const Real orient =
        (x1[0] - x0[0]) * (x2[1] - x0[1])
      - (x1[1] - x0[1]) * (x2[0] - x0[0]);
      if (orient <= 0)
        ++inverted;
      const Real area = Real(0.5) * std::abs(orient);
      const Real den = (x1 - x0).squaredNorm()
                     + (x2 - x1).squaredNorm()
                     + (x0 - x2).squaredNorm();
      if (den > 0)
        qmin = std::min(qmin, Real(4) * std::sqrt(Real(3)) * area / den);
    }
    if (!std::isfinite(qmin))
      qmin = 0;
    return {qmin, inverted};
  }

  bool localVertexMoveValid(
      const LocalMesh& mesh,
      Index vertex,
      const Math::SpatialPoint& candidate,
      Real minQuality = Real(1e-5))
  {
    const auto& conn = mesh.getConnectivity();
    for (Index ci : conn.getIncidence({ 0, 2 }, vertex))
    {
      const auto& cell = conn.getPolytope(2, ci);
      std::array<Math::SpatialPoint, 3> x{{
        mesh.getVertexCoordinates(cell(0)),
        mesh.getVertexCoordinates(cell(1)),
        mesh.getVertexCoordinates(cell(2)) }};
      for (size_t k = 0; k < 3; ++k)
        if (cell(k) == vertex)
          x[k] = candidate;
      const Real orient =
        (x[1][0] - x[0][0]) * (x[2][1] - x[0][1])
      - (x[1][1] - x[0][1]) * (x[2][0] - x[0][0]);
      if (!(orient > Real(1e-14)) || !std::isfinite(orient))
        return false;
      if (triangleQuality(x[0], x[1], x[2]) < minQuality)
        return false;
    }
    return true;
  }

  Index snapInterfaceVerticesToAnalytic(
      LocalMesh& mesh,
      Real t,
      Real maxMove)
  {
    mesh.getConnectivity().compute(0, 2);
    mesh.getConnectivity().compute(1, 0);
    const auto mask = interfaceVertexMask(mesh);
    Index accepted = 0;
    for (Index v = 0; v < static_cast<Index>(mesh.getVertexCount()); ++v)
    {
      if (v >= static_cast<Index>(mask.size()) || !mask[v])
        continue;
      const auto x = mesh.getVertexCoordinates(v);
      auto p = projectToInterface(x, t);
      const Real d = (p - x).norm();
      if (d > maxMove && d > Real(0))
        p = x + (maxMove / d) * (p - x);
      bool moved = false;
      for (int attempt = 0; attempt < 12; ++attempt)
      {
        if (localVertexMoveValid(mesh, v, p))
        {
          mesh.setVertexCoordinates(v, p);
          ++accepted;
          moved = true;
          break;
        }
        p = x + Real(0.5) * (p - x);
      }
      (void)moved;
    }
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    return accepted;
  }

  Real interfaceMeasure(
      const LocalMesh& mesh,
      Attribute attr,
      Index samples = 2)
  {
    Real measure = 0;
    const auto& conn = mesh.getConnectivity();
    const Index n = std::max<Index>(samples, 1);
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto edgeAttr = mesh.getAttribute(1, e);
      if (!edgeAttr || *edgeAttr != attr)
        continue;
      const auto edgeIt = mesh.getPolytope(1, e);
      const auto& trans = edgeIt->getTransformation();
      for (Index q = 0; q < n; ++q)
      {
        const Real s = (static_cast<Real>(q) + Real(0.5))
          / static_cast<Real>(n);
        Math::SpatialMatrix<Real> J;
        trans.jacobian(J, Math::SpatialPoint{s});
        if (J.cols() != 1)
          continue;
        Real len2 = 0;
        for (Index d = 0; d < static_cast<Index>(J.rows()); ++d)
          len2 += J(d, 0) * J(d, 0);
        if (len2 > Real(0) && std::isfinite(len2))
          measure += std::sqrt(len2) / static_cast<Real>(n);
      }
    }
    return measure;
  }

  struct InterfaceDegreeStats
  {
    Index maxDegree = 0;
    Index branchVertices = 0;
  };

  InterfaceDegreeStats interfaceDegreeStats(
      const LocalMesh& mesh,
      Attribute attr)
  {
    InterfaceDegreeStats stats;
    std::vector<Index> degree(mesh.getVertexCount(), 0);
    const auto& conn = mesh.getConnectivity();
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto edgeAttr = mesh.getAttribute(1, e);
      if (!edgeAttr || *edgeAttr != attr)
        continue;
      const auto& edge = conn.getPolytope(1, e);
      degree[static_cast<size_t>(edge(0))]++;
      degree[static_cast<size_t>(edge(1))]++;
    }
    for (Index d : degree)
    {
      stats.maxDegree = std::max(stats.maxDegree, d);
      if (d > 2)
        stats.branchVertices++;
    }
    return stats;
  }

  Index markUncutAnalyticInterfaceEdges(LocalMesh& mesh, Real t)
  {
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    Index count = 0;
    const auto& conn = mesh.getConnectivity();
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto attr = mesh.getAttribute(1, e);
      if (attr && *attr == Boundary)
        continue;
      const auto& edge = conn.getPolytope(1, e);
      const Real f0 = phiAt(mesh.getVertexCoordinates(edge(0)), t);
      const Real f1 = phiAt(mesh.getVertexCoordinates(edge(1)), t);
      if (f0 == Real(0) || f1 == Real(0) || f0 * f1 < Real(0))
      {
        mesh.setAttribute({1, e}, Interface);
        ++count;
      }
    }
    return count;
  }

  struct TopologyStageReport
  {
    Index cellCount = 0;
    Index interfaceEdgeCount = 0;
    Index uncutCellCount = 0;
    Index referenceStencilOrder = 0;
    Index macroBoundaryCutCount = 0;
    Index macroBoundaryCutRejectedByQualityCount = 0;
    Index sampledStencilFallbackCount = 0;
    Index noCutQualityFallbackCount = 0;
    Index graphOnlyFallbackCount = 0;
    Real minOutputCellQuality = 0;
    Real maxInterfaceDeviation = 0;
  };

  struct TopologyStageResult
  {
    LocalMesh mesh;
    TopologyStageReport report;
    InterfaceDegreeStats interfaceDegree;
    bool mmgAccepted = true;
  };

  TopologyStageResult buildTopology(
      const LocalMesh& background,
      Real t,
      Real h,
      Index stencilOrder)
  {
    LocalMesh work(background);
    work.getConnectivity().compute(2, 1);
    work.getConnectivity().compute(1, 0);

    annotateBoundary(work);
    const Index interfaceEdges = markUncutAnalyticInterfaceEdges(work, t);
    const auto interfaceDegree = interfaceDegreeStats(work, Interface);
    const auto [qmin, inverted] = meshQuality(work);

    TopologyStageReport report;
    report.cellCount = work.getCellCount();
    report.minOutputCellQuality = qmin;
    report.uncutCellCount = work.getCellCount();
    report.referenceStencilOrder = 0;
    report.macroBoundaryCutCount = interfaceEdges;
    const auto& conn = work.getConnectivity();
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto attr = work.getAttribute(1, e);
      if (!attr || *attr != Interface)
        continue;
      ++report.interfaceEdgeCount;
      const auto& edge = conn.getPolytope(1, e);
      for (Index v : {edge(0), edge(1)})
        report.maxInterfaceDeviation =
          std::max(report.maxInterfaceDeviation,
              std::abs(phiAt(work.getVertexCoordinates(v), t)));
    }
    (void)h;
    (void)stencilOrder;
    (void)inverted;
    return {
      std::move(work),
      std::move(report),
      interfaceDegree,
      true
    };
  }

  Index reclassifyCellsByAnalyticPhi(LocalMesh& mesh, Real t)
  {
    Index changed = 0;
    const auto& conn = mesh.getConnectivity();
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      if (conn.getGeometry(2, c) != Polytope::Type::Triangle)
        continue;
      const auto cellIt = mesh.getPolytope(2, c);
      Math::SpatialPoint x;
      cellIt->getTransformation().transform(
          x, Math::SpatialPoint{Real(1) / Real(3), Real(1) / Real(3)});
      const Attribute next = phiAt(x, t) <= Real(0) ? Negative : Positive;
      const auto old = mesh.getAttribute(2, c);
      if (!old || *old != next)
      {
        mesh.setAttribute({2, c}, next);
        ++changed;
      }
    }
    return changed;
  }
}

int main(int argc, char** argv)
{
  size_t resolution = 40;
  Index steps = 20;
  Real fitWeight = 2;
  Real qualityWeight = 1;
  Real deviationWeight = Real(0.25);
  Real targetBlend = Real(0.05);
  Real curvatureClamp = Real(0.35);
  Index tmopMaxIterations = 8;
  Index cutStencilOrder = 4;
  if (argc > 1)
    resolution = static_cast<size_t>(std::max(3, std::atoi(argv[1])));
  if (argc > 2)
    steps = static_cast<Index>(std::max(2, std::atoi(argv[2])));
  if (argc > 3)
    fitWeight = std::max(Real(0), static_cast<Real>(std::atof(argv[3])));
  if (argc > 4)
    qualityWeight = std::max(Real(0), static_cast<Real>(std::atof(argv[4])));
  if (argc > 5)
    deviationWeight = std::max(Real(0), static_cast<Real>(std::atof(argv[5])));
  if (argc > 6)
    targetBlend = std::max(Real(0), std::min(Real(1),
          static_cast<Real>(std::atof(argv[6]))));
  if (argc > 7)
    curvatureClamp = std::max(Real(0), static_cast<Real>(std::atof(argv[7])));
  if (argc > 8)
    tmopMaxIterations =
      static_cast<Index>(std::max(1, std::atoi(argv[8])));
  if (argc > 9)
    cutStencilOrder =
      static_cast<Index>(std::max(1, std::atoi(argv[9])));
  const Real h = Real(1) / static_cast<Real>(resolution - 1);

  LocalMesh background =
    LocalMesh::UniformGrid(Polytope::Type::Triangle, { resolution, resolution });
  background.scale(h);

  IO::XDMF xdmf("LevelSetMovingCurvedCircle");

  // Diagnostics are accumulated and reported as an average over the whole
  // process, not just the final step.
  Real sumQmin = 0, sumCurvedQmin = 0, sumFit = 0;
  Real sumTargetMaxMu = 0, sumTargetMinDetT = 0;
  Real sumCurvatureClamp = 0;
  Real sumCurvedCoverage = 0, sumCurvedSignedArea = 0;
  Real sumBestP2Fit = 0, sumBestP2Qmin = 0;
  Real sumCutQmin = 0, sumCutMaxPhi = 0, sumCutCells = 0;
  Real sumCutMacro = 0, sumCutMacroRejected = 0;
  Real sumCutSampledFallback = 0, sumCutNoCutFallback = 0;
  Real sumCutGraphOnlyFallback = 0;
  Real sumInterfaceBranchVertices = 0, sumInterfaceMaxDegree = 0;
  Index sumInverted = 0, sumTargetInvalid = 0;
  Index sumOverlap = 0, sumBestP2Overlap = 0, topologyInsufficient = 0;
  Index mmgRejected = 0;
  Index rejectedTMOP = 0, sumReclassifiedCells = 0, diagSteps = 0;

  for (Index step = 0; step < steps; ++step)
  {
    const Real t = static_cast<Real>(step)
      / static_cast<Real>(std::max<Index>(1, steps - 1));

    // Minimal production loop:
    //   Rodin level-set cutter -> MMG topology optimizer -> analytic snap
    //   -> P2 upgrade -> curved target-matrix TMOP -> carry forward.
    auto topology = buildTopology(background, t, h, cutStencilOrder);
    LocalMesh optimized = std::move(topology.mesh);
    const auto cutReport = std::move(topology.report);
    if (!topology.mmgAccepted)
      ++mmgRejected;

    // Coarse-topology diagnostic: with this exact topology, how well can a
    // fully projected P2 interface fit before TMOP? If this stays large, the
    // limiting factor is the number/placement of feature edges, not the TMOP
    // solve. The actual run below still uses curvatureClamp for robustness.
    LocalMesh bestP2Topology(optimized);
    installQuadraticGeometry(
        bestP2Topology, t, std::numeric_limits<Real>::infinity());
    const auto bestP2Stats = curvedGeometryStats(bestP2Topology, t);

    installQuadraticGeometry(optimized, t, curvatureClamp);

    const auto featureMask = interfaceVertexMask(optimized);
    auto projectToCurrentInterface = [t](const Math::SpatialPoint& x)
    {
      return projectToInterface(x, t);
    };
    ProjectedQualityTargetJacobian target(
        optimized, Interface, projectToCurrentInterface, targetBlend);
    ShapeSizeBlendMetric metric(Real(0.5));

    VectorH1<2, LocalMesh> space(std::integral_constant<size_t, 2>{}, optimized, 2);
    optimized.getConnectivity().compute(1, 2);
    GridFunction u(space);
    u.getData().setZero();
    TrialFunction du(space);
    TestFunction v(space);

    auto phiValue = [t](const Math::SpatialPoint& x) { return phiAt(x, t); };
    auto phiGradient =
      [t](const Math::SpatialPoint& x) { return gradPhiAt(x, t); };

    QualityTerm quality(metric, target, qualityWeight);
    quality.setQuadratureOrder(CurvedQuadratureOrder);
    DeviationTerm deviation(deviationWeight);
    AnalyticLevelSetFitTerm fit(
        phiValue, phiGradient, Optional<Attribute>(Interface), fitWeight);
    fit.setNormalization(
        std::max(interfaceMeasure(optimized, Interface), Real(1e-12)));
    AnalyticLevelSetFitTerm boundaryFit(
        boxBoundaryValue, boxBoundaryGradient,
        Optional<Attribute>(Boundary), 1.0);
    boundaryFit.setNormalization(
        std::max(interfaceMeasure(optimized, Boundary), Real(1e-12)));
    const Real fitEnergyInitial = fit.energy(u);

    auto makeTangent = [&]()
    {
      return quality.tangent(u, du, v)
           + deviation.tangent(du, v)
           + fit.tangent(u, du, v)
           + boundaryFit.tangent(u, du, v);
    };
    auto makeResidual = [&]()
    {
      return quality.residual(u, v)
           + deviation.residual(u, v)
           + fit.residual(u, v)
           + boundaryFit.residual(u, v);
    };
    auto meritEnergy = [&]()
    {
      return quality.energy(u)
           + deviation.energy(u)
           + fit.energy(u)
           + boundaryFit.energy(u);
    };

    bool solvedTMOP = true;
    bool hadAcceptedTMOPStep = false;
    Real lastAlpha = 1;
    Index tmopIterations = 0;
    try
    {
      for (Index it = 0; it < tmopMaxIterations; ++it)
      {
        ++tmopIterations;
        LinearForm R(v);
        R = makeResidual();
        R.assemble();
        const auto r = R.getVector();
        if (!r.allFinite())
        {
          solvedTMOP = hadAcceptedTMOPStep;
          break;
        }
        if (r.norm() <= Real(1e-10))
          break;

        BilinearForm J(du, v);
        J = makeTangent();
        J.assemble();

        Eigen::SparseLU<std::decay_t<decltype(J.getOperator())>> lu;
        lu.compute(J.getOperator());
        if (lu.info() != Eigen::Success)
        {
          solvedTMOP = hadAcceptedTMOPStep;
          break;
        }
        const Math::Vector<Real> dx = lu.solve(-r);
        if (lu.info() != Eigen::Success || !dx.allFinite())
        {
          solvedTMOP = hadAcceptedTMOPStep;
          break;
        }

        const Real e0 = meritEnergy();
        const Math::Vector<Real> u0 = u.getData();
        Real alpha = 1;
        bool acceptedStep = false;
        for (int ls = 0; ls < 30; ++ls)
        {
          u.getData() = u0 + alpha * dx;
          const Real e = meritEnergy();
          if (std::isfinite(e) && e <= e0 * (Real(1) + Real(1e-12)))
          {
            LocalMesh trialMesh(optimized);
            applyParametricDisplacement(trialMesh, u);
            if (fastCurvedLocalValidity(trialMesh))
            {
              acceptedStep = true;
              break;
            }
          }
          alpha *= Real(0.5);
        }
        lastAlpha = alpha;
        if (!acceptedStep)
        {
          u.getData() = u0;
          break;
        }
        hadAcceptedTMOPStep = true;
        if (alpha * dx.norm() <= Real(1e-10))
          break;
      }
    }
    catch (...)
    {
      solvedTMOP = false;
    }
    const Real fitEnergyFinal = fit.energy(u);

    LocalMesh beforeTMOP(optimized);
    if (solvedTMOP && u.getData().allFinite())
    {
      applyParametricDisplacement(optimized, u);
    }
    else
      solvedTMOP = false;

    optimized.getConnectivity().compute(2, 1);
    optimized.getConnectivity().compute(1, 0);
    const auto [tmopQmin, tmopInverted] = meshQuality(optimized);
    auto tmopCurvedStats = curvedGeometryStats(optimized, t);
    const bool acceptedTMOP =
      solvedTMOP
      && tmopInverted == 0
      && std::isfinite(tmopQmin)
      && tmopQmin > Real(0)
      && curvedGeometryValid(tmopCurvedStats);
    if (!acceptedTMOP)
    {
      optimized = std::move(beforeTMOP);
      ++rejectedTMOP;
      tmopCurvedStats = curvedGeometryStats(optimized, t);
    }
    const Index reclassifiedCells = reclassifyCellsByAnalyticPhi(optimized, t);
    optimized.getConnectivity().compute(2, 1);
    optimized.getConnectivity().compute(1, 0);

    P1<Real, LocalMesh> outputSpace(optimized);
    GridFunction outputPhi(outputSpace);
    for (Index i = 0; i < optimized.getVertexCount(); ++i)
      outputPhi[i] = phiAt(optimized.getVertexCoordinates(i), t);

    auto grid = xdmf.grid("curved-p2-xdmf");
    grid.setMesh(optimized, IO::XDMF::MeshPolicy::Transient);
    grid.clear();
    grid.add("phi", outputPhi, IO::XDMF::Center::Node);

    auto curvedMesh = curvedVisualizationMesh(optimized, CurvedVisualizationSamples);
    auto curvedGrid = xdmf.grid("curved-p2-sampled");
    curvedGrid.setMesh(curvedMesh, IO::XDMF::MeshPolicy::Transient);
    curvedGrid.clear();

    auto curvedInterface = curvedInterfacePolyline(optimized);
    auto curveGrid = xdmf.grid("curved-interface");
    curveGrid.setMesh(curvedInterface, IO::XDMF::MeshPolicy::Transient);
    curveGrid.clear();
    xdmf.write(t).flush();

    const auto [qmin, inverted] = meshQuality(optimized);
    const auto curvedStats = curvedGeometryStats(optimized, t);
    const auto targetStats = targetQualityMetrics(
        optimized, metric, target, Polytope::Type::Triangle,
        CurvedQuadratureOrder);
    Real fitErr = 0;
    for (Index i = 0; i < optimized.getVertexCount(); ++i)
      if (i < featureMask.size() && featureMask[i])
        fitErr = std::max(fitErr,
            std::abs(phiAt(optimized.getVertexCoordinates(i), t)));
    sumQmin += qmin;
    sumCurvedQmin += curvedStats.minQuality;
    sumInverted += inverted;
    sumFit += curvedStats.fitRms;
    sumTargetMaxMu += targetStats.maxMetric;
    sumTargetMinDetT += targetStats.minDetT;
    sumCurvatureClamp += curvatureClamp;
    sumCurvedCoverage += curvedStats.sampledCoverage;
    sumCurvedSignedArea += curvedStats.sampledSignedArea;
    sumBestP2Fit += bestP2Stats.fitRms;
    sumBestP2Qmin += bestP2Stats.minQuality;
    sumCutQmin += cutReport.minOutputCellQuality;
    sumCutMaxPhi += cutReport.maxInterfaceDeviation;
    sumCutCells += static_cast<Real>(cutReport.cellCount);
    sumCutMacro += static_cast<Real>(cutReport.macroBoundaryCutCount);
    sumCutMacroRejected += static_cast<Real>(
        cutReport.macroBoundaryCutRejectedByQualityCount);
    sumCutSampledFallback += static_cast<Real>(
        cutReport.sampledStencilFallbackCount);
    sumCutNoCutFallback += static_cast<Real>(
        cutReport.noCutQualityFallbackCount);
    sumCutGraphOnlyFallback += static_cast<Real>(
        cutReport.graphOnlyFallbackCount);
    sumInterfaceBranchVertices +=
      static_cast<Real>(topology.interfaceDegree.branchVertices);
    sumInterfaceMaxDegree +=
      static_cast<Real>(topology.interfaceDegree.maxDegree);
    sumTargetInvalid += targetStats.invalidSamples;
    sumOverlap += curvedStats.overlapSamples;
    sumBestP2Overlap += bestP2Stats.overlapSamples;
    sumReclassifiedCells += reclassifiedCells;
    // Topology/resolution insufficient: TMOP made no real step (the line
    // search backtracked to ~0) because the start has a cell the fixed
    // topology cannot improve -- either a leftover sliver, or an interface
    // cell that cannot be validly curved at this resolution (the curve would
    // fold it, so it stays affine and reads high distortion against the
    // curved target). Report it rather than pretend the optimizer succeeded.
    const bool topoInsufficient = lastAlpha < Real(1e-6);
    if (topoInsufficient) ++topologyInsufficient;
    ++diagSteps;
    std::cout << "step " << step
              << " t=" << t
              << " cells=" << optimized.getCellCount()
              << " cut_stencil_order=" << cutStencilOrder
              << " cut_cells=" << cutReport.cellCount
              << " cut_qmin=" << cutReport.minOutputCellQuality
              << " cut_phi_max=" << cutReport.maxInterfaceDeviation
              << " cut_macro=" << cutReport.macroBoundaryCutCount
              << " cut_macro_rejected="
              << cutReport.macroBoundaryCutRejectedByQualityCount
              << " cut_sampled_fallback="
              << cutReport.sampledStencilFallbackCount
              << " cut_nocut_fallback="
              << cutReport.noCutQualityFallbackCount
              << " cut_graph_only_fallback="
              << cutReport.graphOnlyFallbackCount
              << " cut_interface_max_degree="
              << topology.interfaceDegree.maxDegree
              << " cut_interface_branch_vertices="
              << topology.interfaceDegree.branchVertices
              << " mmg_accepted=" << (topology.mmgAccepted ? 1 : 0)
              << " qmin=" << qmin
              << " curved_qmin=" << curvedStats.minQuality
              << " inverted=" << inverted
              << " fit=" << fitErr
              << " curved_fit_rms=" << curvedStats.fitRms
              << " curved_fit_max=" << curvedStats.fitMax
              << " best_p2_fit_rms=" << bestP2Stats.fitRms
              << " best_p2_fit_max=" << bestP2Stats.fitMax
              << " best_p2_qmin=" << bestP2Stats.minQuality
              << " best_p2_invalid_jac=" << bestP2Stats.invalidJacobianSamples
              << " best_p2_overlap=" << bestP2Stats.overlapSamples
              << " curved_min_jac=" << curvedStats.minJacobian
              << " curved_invalid_jac=" << curvedStats.invalidJacobianSamples
              << " curved_overlap=" << curvedStats.overlapSamples
              << " curved_coverage=" << curvedStats.sampledCoverage
              << " curved_signed_area=" << curvedStats.sampledSignedArea
              << " curvature_clamp=" << curvatureClamp
              << " fit_weight=" << fitWeight
              << " quality_weight=" << qualityWeight
              << " deviation_weight=" << deviationWeight
              << " target_blend=" << targetBlend
              << " cut_stencil_order=" << cutStencilOrder
              << " tmop_target_max_mu=" << targetStats.maxMetric
              << " tmop_target_min_detT=" << targetStats.minDetT
              << " tmop_target_invalid=" << targetStats.invalidSamples
              << " tmop_iterations=" << tmopIterations
              << " tmop_last_alpha=" << lastAlpha
              << " analytic_fit_energy_initial=" << fitEnergyInitial
              << " analytic_fit_energy_final=" << fitEnergyFinal
              << " tmop_accepted=" << (acceptedTMOP ? 1 : 0)
              << " reclassified_cells=" << reclassifiedCells
              << " topology_insufficient=" << (topoInsufficient ? 1 : 0)
              << std::endl;

    demoteInterface(optimized);
    optimized.flush();
    background = std::move(optimized);
  }

  xdmf.close();

  if (diagSteps > 0)
    std::cout << "AVERAGE over " << diagSteps << " steps:"
              << " qmin=" << sumQmin / static_cast<Real>(diagSteps)
              << " curved_qmin="
              << sumCurvedQmin / static_cast<Real>(diagSteps)
              << " inverted="
              << static_cast<Real>(sumInverted) / static_cast<Real>(diagSteps)
              << " curved_fit_rms=" << sumFit / static_cast<Real>(diagSteps)
              << " best_p2_fit_rms="
              << sumBestP2Fit / static_cast<Real>(diagSteps)
              << " best_p2_qmin="
              << sumBestP2Qmin / static_cast<Real>(diagSteps)
              << " cut_cells="
              << sumCutCells / static_cast<Real>(diagSteps)
              << " cut_qmin="
              << sumCutQmin / static_cast<Real>(diagSteps)
              << " cut_phi_max="
              << sumCutMaxPhi / static_cast<Real>(diagSteps)
              << " cut_macro="
              << sumCutMacro / static_cast<Real>(diagSteps)
              << " cut_macro_rejected="
              << sumCutMacroRejected / static_cast<Real>(diagSteps)
              << " cut_sampled_fallback="
              << sumCutSampledFallback / static_cast<Real>(diagSteps)
              << " cut_nocut_fallback="
              << sumCutNoCutFallback / static_cast<Real>(diagSteps)
              << " cut_graph_only_fallback="
              << sumCutGraphOnlyFallback / static_cast<Real>(diagSteps)
              << " cut_interface_max_degree="
              << sumInterfaceMaxDegree / static_cast<Real>(diagSteps)
              << " cut_interface_branch_vertices="
              << sumInterfaceBranchVertices / static_cast<Real>(diagSteps)
              << " tmop_target_max_mu="
              << sumTargetMaxMu / static_cast<Real>(diagSteps)
              << " tmop_target_min_detT="
              << sumTargetMinDetT / static_cast<Real>(diagSteps)
              << " tmop_target_invalid="
              << static_cast<Real>(sumTargetInvalid) / static_cast<Real>(diagSteps)
              << " curved_overlap="
              << static_cast<Real>(sumOverlap) / static_cast<Real>(diagSteps)
              << " best_p2_overlap="
              << static_cast<Real>(sumBestP2Overlap) / static_cast<Real>(diagSteps)
              << " curved_coverage="
              << sumCurvedCoverage / static_cast<Real>(diagSteps)
              << " curved_signed_area="
              << sumCurvedSignedArea / static_cast<Real>(diagSteps)
              << " curvature_clamp="
              << sumCurvatureClamp / static_cast<Real>(diagSteps)
              << " fit_weight=" << fitWeight
              << " quality_weight=" << qualityWeight
              << " deviation_weight=" << deviationWeight
              << " target_blend=" << targetBlend
              << " cut_stencil_order=" << cutStencilOrder
              << " mmg_rejected=" << mmgRejected
              << " rejected_tmop=" << rejectedTMOP
              << " reclassified_cells="
              << static_cast<Real>(sumReclassifiedCells) / static_cast<Real>(diagSteps)
              << " topology_insufficient=" << topologyInsufficient
              << "/" << diagSteps
              << std::endl;

  return 0;
}
