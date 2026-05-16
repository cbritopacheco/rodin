/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_GEOMETRY_H
#define RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_GEOMETRY_H

#include <array>
#include <limits>
#include <unordered_map>
#include <vector>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Geometry/Connectivity.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Math.h"
#include "Rodin/Types.h"

#include "Metrics.h"

namespace Rodin::Adaptation::TargetMatrixOptimization
{
  struct GeometryNode
  {
    /// Current geometry-node position.
    Math::SpatialPoint x;
    /// Initial position used by deviation penalties.
    Math::SpatialPoint x0;
    /// Fixed nodes are not moved by the optimizer.
    bool fixed = false;
  };

  /**
   * @brief Minimal order-2 triangular geometry container.
   *
   * The topology is fixed. Each cell stores six node indices ordered as
   * `{v0, v1, v2, e01, e12, e20}` using quadratic Lagrange triangle shape
   * functions. This representation deliberately separates geometry nodes from
   * mesh vertices so midside and interface geometry can move without changing
   * the fitted topology.
   */
  struct HighOrderTriangleGeometry
  {
    Index order = 2;
    std::vector<GeometryNode> nodes;
    std::vector<std::array<Index, 6>> cells;
    std::unordered_map<Index, Index> originalVertexToNode;
    std::unordered_map<Index, Index> meshEdgeToNode;
  };

  struct ReferencePoint
  {
    Real xi = Real(1) / Real(3);
    Real eta = Real(1) / Real(3);
  };

  class HighOrderGeometryUpgrade
  {
    public:
      /**
       * @brief Controls whether original mesh vertices are fixed after upgrade.
       */
      HighOrderGeometryUpgrade& setFixOriginalVertices(bool fixed)
      {
        m_fixOriginalVertices = fixed;
        return *this;
      }

      /**
       * @brief Upgrades a linear triangular mesh to order-2 geometry.
       *
       * Mid-edge geometry nodes are initialized exactly at linear edge
       * midpoints. The method requires incidences (2,1) and (1,0) to be
       * precomputed, matching the explicit-connectivity policy used by the
       * level-set topology extractor.
       */
      HighOrderTriangleGeometry upgrade(
          const Geometry::LocalMesh& mesh,
          Index order = 2) const
      {
        if (order != 2)
          Alert::MemberFunctionException(*this, __func__)
            << "Only order-2 triangular geometry is supported for now."
            << Alert::Raise;
        if (mesh.getDimension() != 2)
          Alert::MemberFunctionException(*this, __func__)
            << "Expected a two-dimensional triangular mesh."
            << Alert::Raise;

        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 2, 1);
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 1, 0);

        const auto& conn = mesh.getConnectivity();

        HighOrderTriangleGeometry geometry;
        geometry.order = order;
        geometry.nodes.reserve(mesh.getVertexCount() + conn.getCount(1));
        geometry.cells.reserve(mesh.getCellCount());

        for (Index v = 0; v < mesh.getVertexCount(); ++v)
        {
          GeometryNode node;
          node.x = mesh.getVertexCoordinates(v);
          node.x0 = node.x;
          node.fixed = m_fixOriginalVertices;
          geometry.originalVertexToNode.emplace(v, geometry.nodes.size());
          geometry.nodes.push_back(std::move(node));
        }

        auto findCellEdge =
          [&](const Geometry::Polytope::Key& cell,
              const IndexVector& edges,
              std::uint8_t i,
              std::uint8_t j)
          {
            const Index vi = cell(i);
            const Index vj = cell(j);
            for (Index edge : edges)
            {
              const auto& e = conn.getPolytope(1, edge);
              if ((e(0) == vi && e(1) == vj) ||
                  (e(0) == vj && e(1) == vi))
                return edge;
            }
            return std::numeric_limits<Index>::max();
          };

        auto midpointNode = [&](Index edge)
        {
          const auto it = geometry.meshEdgeToNode.find(edge);
          if (it != geometry.meshEdgeToNode.end())
            return it->second;

          const auto& e = conn.getPolytope(1, edge);
          const auto& a = mesh.getVertexCoordinates(e(0));
          const auto& b = mesh.getVertexCoordinates(e(1));

          GeometryNode node;
          node.x = Math::SpatialPoint({
            (a[0] + b[0]) / Real(2),
            (a[1] + b[1]) / Real(2) });
          node.x0 = node.x;
          node.fixed = false;

          const Index idx = static_cast<Index>(geometry.nodes.size());
          geometry.meshEdgeToNode.emplace(edge, idx);
          geometry.nodes.push_back(std::move(node));
          return idx;
        };

        for (Index c = 0; c < mesh.getCellCount(); ++c)
        {
          if (conn.getGeometry(2, c) != Geometry::Polytope::Type::Triangle)
            Alert::MemberFunctionException(*this, __func__)
              << "Expected a triangular mesh."
              << Alert::Raise;

          const auto& cell = conn.getPolytope(2, c);
          const auto& cellEdges = conn.getIncidence({ 2, 1 }, c);

          const Index e01 = findCellEdge(cell, cellEdges, 0, 1);
          const Index e12 = findCellEdge(cell, cellEdges, 1, 2);
          const Index e20 = findCellEdge(cell, cellEdges, 2, 0);
          if (e01 == std::numeric_limits<Index>::max() ||
              e12 == std::numeric_limits<Index>::max() ||
              e20 == std::numeric_limits<Index>::max())
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Could not recover all triangle edge incidences."
              << Alert::Raise;
          }

          geometry.cells.push_back({{
            geometry.originalVertexToNode.at(cell(0)),
            geometry.originalVertexToNode.at(cell(1)),
            geometry.originalVertexToNode.at(cell(2)),
            midpointNode(e01),
            midpointNode(e12),
            midpointNode(e20) }});
        }

        return geometry;
      }

    private:
      bool m_fixOriginalVertices = true;
  };

  class CurvedTriangleJacobianEvaluator
  {
    public:
      /**
       * @brief Evaluates the 2x2 Jacobian of a P2 triangular geometry map.
       */
      Matrix2 jacobian(
          const HighOrderTriangleGeometry& geometry,
          Index cell,
          const ReferencePoint& p = {}) const
      {
        const auto& c = geometry.cells[cell];
        const Real L0 = Real(1) - p.xi - p.eta;
        const Real L1 = p.xi;
        const Real L2 = p.eta;

        const std::array<Real, 6> dNdxi = {{
          Real(1) - Real(4) * L0,
          Real(4) * L1 - Real(1),
          Real(0),
          Real(4) * (L0 - L1),
          Real(4) * L2,
          -Real(4) * L2 }};

        const std::array<Real, 6> dNdeta = {{
          Real(1) - Real(4) * L0,
          Real(0),
          Real(4) * L2 - Real(1),
          -Real(4) * L1,
          Real(4) * L1,
          Real(4) * (L0 - L2) }};

        Matrix2 J = Matrix2::Zero();
        for (std::uint8_t i = 0; i < 6; ++i)
        {
          const auto& x = geometry.nodes[c[i]].x;
          J(0, 0) += x[0] * dNdxi[i];
          J(1, 0) += x[1] * dNdxi[i];
          J(0, 1) += x[0] * dNdeta[i];
          J(1, 1) += x[1] * dNdeta[i];
        }
        return J;
      }

      /**
       * @brief Evaluates the determinant of the P2 geometry Jacobian.
       */
      Real determinant(
          const HighOrderTriangleGeometry& geometry,
          Index cell,
          const ReferencePoint& p = {}) const
      {
        return jacobian(geometry, cell, p).determinant();
      }

      /**
       * @brief Returns the minimum determinant over all cells and sample points.
       */
      Real minJacobian(
          const HighOrderTriangleGeometry& geometry,
          const std::vector<ReferencePoint>& points) const
      {
        Real min = std::numeric_limits<Real>::infinity();
        for (Index c = 0; c < geometry.cells.size(); ++c)
          for (const auto& p : points)
            min = std::min(min, determinant(geometry, c, p));
        return min;
      }
  };
}

#endif
