/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_ISOPARAMETRICGEOMETRY_H
#define RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_ISOPARAMETRICGEOMETRY_H

/**
 * @file
 * @brief Free functions for isoparametric mesh adaptation.
 *
 * The whole curved-mesh pipeline reduces to a handful of stateless
 * operations on a LocalMesh and an arbitrary-FES displacement
 * GridFunction:
 *
 *   upgradeTransformations<Element>(mesh)  // linear -> affine Pk parametric
 *   moveMesh         (mesh, u)             // FES-INDEPENDENT displace
 *   syncLinearBackbone<Element>(mesh)      // P2 corners -> linear vertex coords
 *   demoteTransformations(mesh)            // strip Pk
 *   curvedMetrics    (mesh, phi, attr)     // diagnostics
 *   enforceInterfaceTangentOnU(mesh, u, projector, attr)
 *   isCurvedMoveValid(mesh, u, floor)
 *
 * moveMesh is FES-independent: it samples u(p) at each parametric node's
 * current physical position via GridFunction::getValue and adds it to the
 * node. P1 displacement moves P2 elements, P2 moves P1, etc. The only
 * requirement is that the mesh element be Lagrange-nodal at its reference
 * nodes (which RealH1Element<K> is, by construction on Fekete points).
 */

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/ParametricTransformation.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/Geometry/PointCloud.h"
#include "Rodin/Math.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/Types.h"
#include "Rodin/Variational/GridFunction.h"

namespace Rodin::Adaptation::TargetMatrixOptimization
{
  /// Curved-truth diagnostics sampled directly through the cell + edge
  /// parametric transformations.
  struct CurvedMetrics
  {
    Real fitRms = 0;
    Real fitMax = 0;
    Real qmin = std::numeric_limits<Real>::infinity();
    Real minDet = std::numeric_limits<Real>::infinity();
    Index overlapSamples = 0;
  };

  namespace Detail
  {
    inline Real triQuality(
        const Math::SpatialPoint& a,
        const Math::SpatialPoint& b,
        const Math::SpatialPoint& c)
    {
      const Real orient = (b[0] - a[0]) * (c[1] - a[1])
                        - (b[1] - a[1]) * (c[0] - a[0]);
      const Real area = Real(0.5) * std::abs(orient);
      const Real den = (b - a).squaredNorm()
                     + (c - b).squaredNorm()
                     + (a - c).squaredNorm();
      if (den <= Real(0)) return Real(0);
      return Real(4) * std::sqrt(Real(3)) * area / den;
    }

    inline bool pointInLinearTriangle(
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

    /// Build the affine PointCloud for one cell: each scalar reference node
    /// j is mapped through the cell's three linear vertices.
    template <class Element>
    Geometry::PointCloud affineCellPointCloud(
        const Geometry::LocalMesh& mesh,
        Index ci,
        const Element& /*fe*/)
    {
      const auto& conn = mesh.getConnectivity();
      const auto& cell = conn.getPolytope(2, ci);
      const auto Va = mesh.getVertexCoordinates(cell(0));
      const auto Vb = mesh.getVertexCoordinates(cell(1));
      const auto Vc = mesh.getVertexCoordinates(cell(2));
      const auto& ref =
        Element::getNodes(Geometry::Polytope::Type::Triangle);
      Geometry::PointCloud pc(2, ref.size());
      for (size_t j = 0; j < ref.size(); ++j)
      {
        const Real r = ref[j][0];
        const Real s = ref[j][1];
        pc(0, j) = (1 - r - s) * Va[0] + r * Vb[0] + s * Vc[0];
        pc(1, j) = (1 - r - s) * Va[1] + r * Vb[1] + s * Vc[1];
      }
      return pc;
    }

    /// Build the affine PointCloud for one edge using the segment element's
    /// reference nodes mapped through the edge's two linear endpoints.
    template <class Element>
    Geometry::PointCloud affineEdgePointCloud(
        const Geometry::LocalMesh& mesh,
        Index e,
        const Element& /*feSeg*/)
    {
      const auto& conn = mesh.getConnectivity();
      const auto& edge = conn.getPolytope(1, e);
      const auto a = mesh.getVertexCoordinates(edge(0));
      const auto b = mesh.getVertexCoordinates(edge(1));
      const auto& ref =
        Element::getNodes(Geometry::Polytope::Type::Segment);
      Geometry::PointCloud pc(2, ref.size());
      for (size_t j = 0; j < ref.size(); ++j)
      {
        const Real r = ref[j][0];
        pc(0, j) = (Real(1) - r) * a[0] + r * b[0];
        pc(1, j) = (Real(1) - r) * a[1] + r * b[1];
      }
      return pc;
    }

    /// Identify which scalar local of `Element` corresponds to each of the
    /// three triangle corners (refs (0,0), (1,0), (0,1)) by closest-match.
    template <class Element>
    std::array<size_t, 3> cornerScalarLocals()
    {
      static const std::array<Math::SpatialPoint, 3> kCorner = {{
        Math::SpatialPoint{Real(0), Real(0)},
        Math::SpatialPoint{Real(1), Real(0)},
        Math::SpatialPoint{Real(0), Real(1)} }};
      const auto& nodes =
        Element::getNodes(Geometry::Polytope::Type::Triangle);
      std::array<size_t, 3> out{nodes.size(), nodes.size(), nodes.size()};
      for (size_t k = 0; k < 3; ++k)
      {
        Real best = std::numeric_limits<Real>::infinity();
        for (size_t j = 0; j < nodes.size(); ++j)
        {
          const Real dr = nodes[j][0] - kCorner[k][0];
          const Real ds = nodes[j][1] - kCorner[k][1];
          const Real d = dr * dr + ds * ds;
          if (d < best) { best = d; out[k] = j; }
        }
      }
      return out;
    }
  }

  // =========================================================================
  // Stage operations -- free functions, no class state.
  // =========================================================================

  /// Linear-to-affine Pk lift. For every triangle cell, attach a
  /// ParametricTransformation<Element> whose PointCloud is the affine map
  /// of the element's reference nodes through the cell's three linear
  /// vertex coords. The geometry IS the linear mesh, just represented in
  /// the Pk basis. Conformity is trivial -- every parametric node's
  /// position is a pure function of the shared linear vertex coords.
  ///
  /// Also attaches an Edge_K ParametricTransformation on every edge tagged
  /// with @p interfaceAttribute, so downstream curved-edge integrators have
  /// a P2 edge to integrate along.
  template <class Element>
  void upgradeTransformations(
      Geometry::LocalMesh& mesh,
      const Element& fe,
      Geometry::Attribute interfaceAttribute = 0)
  {
    const auto& conn = mesh.getConnectivity();
    const Index nc = static_cast<Index>(mesh.getCellCount());
    for (Index ci = 0; ci < nc; ++ci)
    {
      if (conn.getGeometry(2, ci) != Geometry::Polytope::Type::Triangle)
        continue;
      mesh.setPolytopeTransformation({ size_t(2), ci },
          new Geometry::ParametricTransformation<Element>(
            Detail::affineCellPointCloud(mesh, ci, fe), fe));
    }

    if (interfaceAttribute == 0)
      return;
    Element feSeg(Geometry::Polytope::Type::Segment);
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto attr = mesh.getAttribute(1, e);
      if (!attr || *attr != interfaceAttribute) continue;
      mesh.setPolytopeTransformation({ size_t(1), e },
          new Geometry::ParametricTransformation<Element>(
            Detail::affineEdgePointCloud(mesh, e, feSeg), feSeg));
    }
  }

  /// Strip every cached parametric transformation. Subsequent
  /// getPolytopeTransformation() calls re-derive the default affine
  /// transform from the current linear vertex coords.
  inline void demoteTransformations(Geometry::LocalMesh& mesh)
  {
    mesh.flush();
  }

  /// FES-INDEPENDENT displacement application. For every triangle cell,
  /// rebuild its parametric transformation by sampling @p u at each
  /// reference node's CURRENT physical position and adding the sampled
  /// displacement to that node.
  ///
  /// The only requirement on @p Element is Lagrange-nodality at its
  /// reference nodes (so that the PointCloud column at node j IS the
  /// physical position at j). @p u may live in any FES (any order, any
  /// vdim>=sdim): GridFunction::getValue does the evaluation.
  ///
  /// Curved interface-edge transformations (if attached) are refreshed
  /// the same way.
  template <class Element, class FES, class Data>
  void moveMesh(
      Geometry::LocalMesh& mesh,
      const Variational::GridFunction<FES, Data>& u,
      const Element& fe,
      Geometry::Attribute interfaceAttribute = 0)
  {
    const auto& conn = mesh.getConnectivity();
    const auto& ref =
      Element::getNodes(Geometry::Polytope::Type::Triangle);
    const Index nc = static_cast<Index>(mesh.getCellCount());
    for (Index ci = 0; ci < nc; ++ci)
    {
      if (conn.getGeometry(2, ci) != Geometry::Polytope::Type::Triangle)
        continue;
      auto cellIt = mesh.getPolytope(2, ci);
      Geometry::PointCloud pc(2, ref.size());
      for (size_t j = 0; j < ref.size(); ++j)
      {
        const Geometry::Point p(*cellIt, ref[j]);
        const auto phys = p.getPhysicalCoordinates();
        const auto uVal = u.getValue(p);
        pc(0, j) = phys[0] + uVal[0];
        pc(1, j) = phys[1] + uVal[1];
      }
      mesh.setPolytopeTransformation({ size_t(2), ci },
          new Geometry::ParametricTransformation<Element>(pc, fe));
    }

    if (interfaceAttribute == 0)
      return;
    Element feSeg(Geometry::Polytope::Type::Segment);
    const auto& refSeg =
      Element::getNodes(Geometry::Polytope::Type::Segment);
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto attr = mesh.getAttribute(1, e);
      if (!attr || *attr != interfaceAttribute) continue;
      auto edgeIt = mesh.getPolytope(1, e);
      Geometry::PointCloud pc(2, refSeg.size());
      for (size_t j = 0; j < refSeg.size(); ++j)
      {
        const Geometry::Point p(*edgeIt, refSeg[j]);
        const auto phys = p.getPhysicalCoordinates();
        const auto uVal = u.getValue(p);
        pc(0, j) = phys[0] + uVal[0];
        pc(1, j) = phys[1] + uVal[1];
      }
      mesh.setPolytopeTransformation({ size_t(1), e },
          new Geometry::ParametricTransformation<Element>(pc, feSeg));
    }
  }

  /// Sync the linear mesh vertex coordinates to the current parametric
  /// corner positions. Use AFTER moveMesh, before any classify / coarsen /
  /// carry-forward / next-step-cut that needs a consistent linear backbone.
  /// setVertexCoordinates flushes the cached transformations, so they are
  /// reattached afterwards from the current PointClouds.
  template <class Element>
  void syncLinearBackbone(Geometry::LocalMesh& mesh, const Element& /*fe*/)
  {
    const auto cornerLocal = Detail::cornerScalarLocals<Element>();
    const auto& conn = mesh.getConnectivity();
    const Index nc = static_cast<Index>(mesh.getCellCount());
    // First read every cell's corner positions (before any mutation).
    std::vector<std::array<Math::SpatialPoint, 3>> cornerPhys(nc);
    for (Index ci = 0; ci < nc; ++ci)
    {
      if (conn.getGeometry(2, ci) != Geometry::Polytope::Type::Triangle)
        continue;
      auto cellIt = mesh.getPolytope(2, ci);
      const auto& ref =
        Element::getNodes(Geometry::Polytope::Type::Triangle);
      for (std::uint8_t k = 0; k < 3; ++k)
        cornerPhys[ci][k] =
          Geometry::Point(*cellIt, ref[cornerLocal[k]])
            .getPhysicalCoordinates();
    }
    // Then overwrite the linear backbone vertices.
    for (Index ci = 0; ci < nc; ++ci)
    {
      if (conn.getGeometry(2, ci) != Geometry::Polytope::Type::Triangle)
        continue;
      const auto& cell = conn.getPolytope(2, ci);
      for (std::uint8_t k = 0; k < 3; ++k)
        mesh.setVertexCoordinates(cell(k), cornerPhys[ci][k]);
    }
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
  }

  /// Curved validity gate: would moving @p mesh by @p u produce a curved
  /// Pk map whose det(grad x_h) on an order-4 quadrature is strictly above
  /// @p qmin_floor everywhere?
  template <class Element, class FES, class Data>
  bool isCurvedMoveValid(
      const Geometry::LocalMesh& mesh,
      const Variational::GridFunction<FES, Data>& u,
      const Element& fe,
      Real qmin_floor = Real(0))
  {
    const auto& qf = QF::PolytopeQuadratureFormula::get(
        4, Geometry::Polytope::Type::Triangle);
    const auto& ref =
      Element::getNodes(Geometry::Polytope::Type::Triangle);
    const auto& conn = mesh.getConnectivity();
    const Index nc = static_cast<Index>(mesh.getCellCount());
    for (Index ci = 0; ci < nc; ++ci)
    {
      if (conn.getGeometry(2, ci) != Geometry::Polytope::Type::Triangle)
        continue;
      auto cellIt = mesh.getPolytope(2, ci);
      Geometry::PointCloud pc(2, ref.size());
      for (size_t j = 0; j < ref.size(); ++j)
      {
        const Geometry::Point p(*cellIt, ref[j]);
        const auto phys = p.getPhysicalCoordinates();
        const auto uVal = u.getValue(p);
        pc(0, j) = phys[0] + uVal[0];
        pc(1, j) = phys[1] + uVal[1];
      }
      Geometry::ParametricTransformation<Element> trial(pc, fe);
      for (size_t q = 0; q < qf.getSize(); ++q)
      {
        Math::SpatialMatrix<Real> J;
        trial.jacobian(J, qf.getPoint(q));
        const Real det = J.determinant();
        if (!(det > qmin_floor) || !std::isfinite(det))
          return false;
      }
    }
    return true;
  }

  /// Curved-truth metrics. phi is supplied by the caller; this function
  /// has no other state.
  template <class PhiFn>
  CurvedMetrics curvedMetrics(
      Geometry::LocalMesh& mesh,
      PhiFn phi,
      Geometry::Attribute interfaceAttribute)
  {
    CurvedMetrics m;
    static const std::array<Real, 5> sSamples = {{
      Real(0.1), Real(0.25), Real(0.5), Real(0.75), Real(0.9) }};
    const auto& conn = mesh.getConnectivity();

    // fit along curved interface edges
    Real sq = 0;
    Index nFit = 0;
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto attr = mesh.getAttribute(1, e);
      if (!attr || *attr != interfaceAttribute) continue;
      auto it = mesh.getPolytope(1, e);
      for (Real s : sSamples)
      {
        const auto x =
          Geometry::Point(*it, Math::SpatialPoint{ s })
            .getPhysicalCoordinates();
        const Real p = std::abs(phi(x));
        m.fitMax = std::max(m.fitMax, p);
        sq += p * p;
        ++nFit;
      }
    }
    if (nFit > 0)
      m.fitRms = std::sqrt(sq / static_cast<Real>(nFit));

    // min det + sub-corner shape quality
    const auto& qf = QF::PolytopeQuadratureFormula::get(
        4, Geometry::Polytope::Type::Triangle);
    constexpr int S = 2;
    const Index nc = static_cast<Index>(mesh.getCellCount());
    for (Index ci = 0; ci < nc; ++ci)
    {
      if (conn.getGeometry(2, ci) != Geometry::Polytope::Type::Triangle)
        continue;
      auto it = mesh.getPolytope(2, ci);
      for (size_t q = 0; q < qf.getSize(); ++q)
        m.minDet = std::min(
            m.minDet,
            Geometry::Point(*it, qf.getPoint(q))
              .getJacobian().determinant());
      auto sample = [&](int i, int j)
      {
        return Geometry::Point(*it, Math::SpatialPoint{
          static_cast<Real>(i) / static_cast<Real>(S),
          static_cast<Real>(j) / static_cast<Real>(S) })
            .getPhysicalCoordinates();
      };
      for (int i = 0; i < S; ++i)
        for (int j = 0; j + i < S; ++j)
        {
          const auto a = sample(i, j);
          const auto b = sample(i + 1, j);
          const auto c = sample(i, j + 1);
          m.qmin = std::min(m.qmin, Detail::triQuality(a, b, c));
          if (i + j + 1 < S)
          {
            const auto d = sample(i + 1, j + 1);
            m.qmin = std::min(m.qmin, Detail::triQuality(b, d, c));
          }
        }
    }
    if (!std::isfinite(m.qmin)) m.qmin = 0;
    if (!std::isfinite(m.minDet)) m.minDet = 0;

    // overlap canary
    mesh.getConnectivity().compute(0, 2);
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto attr = mesh.getAttribute(1, e);
      if (!attr || *attr != interfaceAttribute) continue;
      const auto& ep = conn.getPolytope(1, e);
      std::array<Index, 2> edgeAdj{
        std::numeric_limits<Index>::max(),
        std::numeric_limits<Index>::max() };
      {
        size_t k = 0;
        for (Index ci : conn.getIncidence({ 1, 2 }, e))
          if (k < 2) edgeAdj[k++] = ci;
      }
      std::vector<Index> nbrs;
      for (Index vEnd : { ep(0), ep(1) })
        for (Index ci : conn.getIncidence({ 0, 2 }, vEnd))
          if (ci != edgeAdj[0] && ci != edgeAdj[1])
            nbrs.push_back(ci);
      std::sort(nbrs.begin(), nbrs.end());
      nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
      if (nbrs.empty()) continue;
      auto edgeIt = mesh.getPolytope(1, e);
      for (Real s : sSamples)
      {
        const auto sample =
          Geometry::Point(*edgeIt, Math::SpatialPoint{ s })
            .getPhysicalCoordinates();
        for (Index ci : nbrs)
        {
          const auto& cell = conn.getPolytope(2, ci);
          const auto Va = mesh.getVertexCoordinates(cell(0));
          const auto Vb = mesh.getVertexCoordinates(cell(1));
          const auto Vc = mesh.getVertexCoordinates(cell(2));
          if (Detail::pointInLinearTriangle(sample, Va, Vb, Vc))
          { ++m.overlapSamples; break; }
        }
      }
    }
    return m;
  }
}

#endif
