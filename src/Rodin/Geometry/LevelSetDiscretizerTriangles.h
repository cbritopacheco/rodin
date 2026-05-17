/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_LEVELSETDISCRETIZERTRIANGLES_H
#define RODIN_GEOMETRY_LEVELSETDISCRETIZERTRIANGLES_H

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Variational/P1/P1.h"

#include "LevelSetInterfaceGraph.h"
#include "Mesh.h"
#include "Polytope.h"
#include "Types.h"

namespace Rodin::Geometry
{
  enum class LevelSetSide
  {
    /// Output cell lies on the negative side of the P1 level set.
    Negative,
    /// Output cell lies on the positive side of the P1 level set.
    Positive,
    /// Reserved for lower-dimensional interface entities.
    Interface,
    /// Output cell came from an invalid or all-zero parent cell.
    Degenerate,
    /// Side could not be assigned.
    Unknown
  };

  /**
   * @brief Describes whether an output vertex is inherited or cut-created.
   *
   * The triangle cutter never creates vertices from quality heuristics. Output
   * vertices are either original background mesh vertices or vertices already
   * present in the InterfaceGraph extracted from the P1 level set.
   */
  enum class OutputVertexOriginKind
  {
    OriginalMeshVertex,
    InterfaceGraphVertex
  };

  /**
   * @brief Provenance for one output mesh vertex.
   */
  struct OutputVertexOrigin
  {
    OutputVertexOriginKind kind = OutputVertexOriginKind::OriginalMeshVertex;
    Optional<Index> originalVertex;
    Optional<Index> interfaceGraphVertex;
  };

  /**
   * @brief Provenance and side classification for one output triangle.
   */
  struct OutputCellProvenance
  {
    Index parentCell = std::numeric_limits<Index>::max();
    LevelSetSide side = LevelSetSide::Unknown;
  };

  /**
   * @brief Parent-cell reference coordinates for one output triangle.
   *
   * vertexBarycentric[i] stores the barycentric coordinates, in the parent
   * background triangle, of local output vertex i. This is the cell-local map
   * needed to transfer a background FunctionBase/GridFunction onto the fitted
   * mesh without guessing from global vertex provenance alone.
   */
  struct OutputCellReference
  {
    Index parentCell = std::numeric_limits<Index>::max();
    std::array<std::array<Real, 3>, 3> vertexBarycentric = {{
      {{Real(0), Real(0), Real(0)}},
      {{Real(0), Real(0), Real(0)}},
      {{Real(0), Real(0), Real(0)}} }};
  };

  /**
   * @brief Provenance for one output mesh edge constrained to the interface.
   *
   * The map in LevelSetDiscretizerTrianglesReport is keyed by the output mesh
   * edge index. The value points back to the source InterfaceGraph edge and, if
   * available, the first parent background cell that generated that segment.
   */
  struct OutputInterfaceEdgeProvenance
  {
    Index sourceInterfaceGraphEdge = std::numeric_limits<Index>::max();
    Optional<Index> parentCell;
  };

  /**
   * @brief Audit trail and diagnostics produced by LevelSetDiscretizerTriangles.
   *
   * The report is meant to feed later geometry optimization and fallback logic:
   * it records how vertices and cells relate to the background mesh and
   * InterfaceGraph, which output edges are protected interface constraints, and
   * simple quality/pathology indicators. Diagnostics are non-fatal; this class
   * is a topology generator and does not try to repair poor cuts.
   */
  struct LevelSetDiscretizerTrianglesReport
  {
    std::vector<OutputVertexOrigin> vertexOrigins;
    std::vector<OutputCellProvenance> cellProvenance;
    std::vector<OutputCellReference> cellReferences;

    std::unordered_map<Index, Index> originalVertexToOutputVertex;
    std::unordered_map<Index, Index> graphVertexToOutputVertex;
    std::unordered_map<Index, OutputInterfaceEdgeProvenance> interfaceEdgeProvenance;

    Index degenerateCellCount = 0;
    Index pathologicalCutCount = 0;
    Index uncutCellCount = 0;
    Index nearVertexCrossingCount = 0;
    Index snappedCrossingCount = 0;
    Index improvedPolygonTriangulationCount = 0;
    Real minOutputCellArea = 0;
    Real minOutputCellQuality = 0;
    Real maxSnapDistance = 0;
    Real maxInterfaceDeviation = 0;
  };

  struct LevelSetDiscretizerTrianglesResult
  {
    /// Body-fitted linear mesh with explicit interface edges.
    LocalMesh mesh;
    /// P1 interface graph used as the source of fitted topology.
    InterfaceGraph interfaceGraph;
    /// Provenance and cut diagnostics for the output mesh.
    LevelSetDiscretizerTrianglesReport report;
  };

  template <class... Params>
  class LevelSetDiscretizerTriangles;

  /**
   * Builds a linear body-fitted triangular topology from a P1 level set.
   *
   * This class is intentionally a topology generator, not a quality optimizer.
   * Poor-quality triangles can be produced near small cuts; they are reported so
   * a later high-order geometry upgrade and TMOP-style optimizer can recover
   * quality while preserving the P1 interface topology.
   *
   * All-zero and invalid parent triangles are preserved as degenerate output
   * cells and counted in the report. This keeps the output mesh auditable while
   * avoiding an unsupported topological decision in this local cutter.
   */
  template <class Mesh, class Data>
  class LevelSetDiscretizerTriangles<
    Variational::GridFunction<Variational::P1<Real, Mesh>, Data>>
  {
    public:
      using MeshType = Mesh;
      using FESType = Variational::P1<Real, Mesh>;
      using GridFunctionType = Variational::GridFunction<FESType, Data>;
      using ResultType = LevelSetDiscretizerTrianglesResult;

      static constexpr Index InvalidIndex = std::numeric_limits<Index>::max();

      explicit LevelSetDiscretizerTriangles(const GridFunctionType& phi)
        : m_phi(phi)
      {}

      LevelSetDiscretizerTriangles& setSignTolerance(Real tol)
      {
        m_signTolerance = std::max(Real(0), tol);
        return *this;
      }

      /**
       * @brief Sets the optional attribute assigned to generated interface edges.
       */
      LevelSetDiscretizerTriangles& setInterfaceAttribute(
          const Optional<Attribute>& attr)
      {
        m_interfaceAttribute = attr;
        return *this;
      }

      LevelSetDiscretizerTriangles& setNegativeCellAttribute(
          const Optional<Attribute>& attr)
      {
        m_negativeCellAttribute = attr;
        return *this;
      }

      LevelSetDiscretizerTriangles& setPositiveCellAttribute(
          const Optional<Attribute>& attr)
      {
        m_positiveCellAttribute = attr;
        return *this;
      }

      LevelSetDiscretizerTriangles& setDiagnosticTolerance(Real tol)
      {
        m_diagnosticTolerance = std::max(Real(0), tol);
        return *this;
      }

      /**
       * @brief Snaps strict edge crossings close to a background vertex.
       *
       * The tolerance is measured as an edge interpolation fraction. The
       * default value zero preserves the exact P1 cut. Positive values collapse
       * near-vertex crossings onto the existing endpoint and report the bounded
       * interface deviation, eliminating the most common sliver topology before
       * TMOP sees the mesh.
       */
      LevelSetDiscretizerTriangles& setCrossingSnapTolerance(Real tol)
      {
        m_crossingSnapTolerance = std::max(Real(0), tol);
        return *this;
      }

      LevelSetDiscretizerTriangles& setAreaTolerance(Real tol)
      {
        m_areaTolerance = std::max(Real(0), tol);
        return *this;
      }

      LevelSetDiscretizerTriangles& setQualityTolerance(Real tol)
      {
        m_qualityTolerance = std::max(Real(0), tol);
        return *this;
      }

      /// Triangle quality measure: (p0, p1, p2) -> quality in [0, 1].
      using QualityMeasure = std::function<Real(
          const Math::SpatialPoint&,
          const Math::SpatialPoint&,
          const Math::SpatialPoint&)>;

      /**
       * @brief Overrides the triangle quality measure used by the robust cut.
       *
       * Defaults to the normalized shape measure
       * @f$4\sqrt3\,A/\sum \ell^2 @f$. The measure drives the quad-diagonal
       * choice, the reported min quality, and the no-cut quality fallback.
       */
      LevelSetDiscretizerTriangles& setQualityMeasure(QualityMeasure measure)
      {
        m_qualityMeasure = std::move(measure);
        return *this;
      }

      /**
       * @brief Minimum child-triangle quality required to perform a cut.
       *
       * Default zero preserves the exact cut (every crossed cell is split).
       * When positive, a cell whose split would produce a child triangle below
       * this quality is left whole on its dominant sign side instead of being
       * cut. The interface is then not body-fitted through that cell (a later
       * TargetMatrixOptimization fit recovers it); the cell is counted in
       * report.uncutCellCount.
       */
      LevelSetDiscretizerTriangles& setMinCutQuality(Real q)
      {
        m_minCutQuality = std::max(Real(0), q);
        return *this;
      }

      /**
       * @brief Builds the fitted linear mesh.
       *
       * Preconditions:
       * - the input mesh must be two-dimensional and triangular;
       * - incidences (2,1) and (1,0) must already be computed by the caller.
       *
       * Guarantees:
       * - original vertices are reused;
       * - InterfaceGraph cut vertices are reused globally;
       * - every InterfaceGraph edge is represented as an output mesh edge;
       * - output triangles are classified by the P1 sign side of their parent
       *   cell;
       * - the report contains enough provenance for later fixed-topology
       *   high-order/TMOP processing.
       */
      ResultType discretize() const
      {
        const auto& phi = m_phi.get();
        const auto& mesh = phi.getFiniteElementSpace().getMesh();

        if (mesh.getDimension() != 2)
          Alert::MemberFunctionException(*this, __func__)
            << "Expected a mesh of topological dimension 2."
            << Alert::Raise;

        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 2, 1);
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 1, 0);

        const auto& conn = mesh.getConnectivity();

        InterfaceGraph graph = LevelSetInterfaceGraph(phi)
          .setSignTolerance(m_signTolerance)
          .setCrossingSnapTolerance(m_crossingSnapTolerance)
          .setInterfaceAttribute(m_interfaceAttribute)
          .extract();

        ResultType result;
        result.interfaceGraph = graph;
        auto& report = result.report;
        report.nearVertexCrossingCount = graph.nearVertexCrossingCount;
        report.snappedCrossingCount = graph.snappedCrossingCount;
        report.maxSnapDistance = graph.maxSnapDistance;
        report.maxInterfaceDeviation = graph.maxInterfaceDeviation;

        const Index originalVertexCount =
          static_cast<Index>(mesh.getVertexCount());
        std::vector<Math::SpatialPoint> outputVertices;
        outputVertices.reserve(originalVertexCount + graph.vertices.size());

        for (Index v = 0; v < originalVertexCount; ++v)
        {
          report.originalVertexToOutputVertex.emplace(v, v);
          outputVertices.push_back(mesh.getVertexCoordinates(v));

          OutputVertexOrigin origin;
          origin.kind = OutputVertexOriginKind::OriginalMeshVertex;
          origin.originalVertex = v;
          report.vertexOrigins.push_back(std::move(origin));
        }

        std::vector<Index> graphVertexToOutput(graph.vertices.size(), InvalidIndex);
        std::unordered_map<Index, Index> meshEdgeToGraphVertex;
        for (Index gv = 0; gv < graph.vertices.size(); ++gv)
        {
          const auto& vertex = graph.vertices[gv];
          if (vertex.kind == InterfaceVertexKind::OriginalVertex)
          {
            if (vertex.originalVertex)
              graphVertexToOutput[gv] = *vertex.originalVertex;
          }
          else
          {
            const Index out = static_cast<Index>(outputVertices.size());
            graphVertexToOutput[gv] = out;
            outputVertices.push_back(vertex.x);

            OutputVertexOrigin origin;
            origin.kind = OutputVertexOriginKind::InterfaceGraphVertex;
            origin.interfaceGraphVertex = gv;
            report.vertexOrigins.push_back(std::move(origin));

            if (vertex.parentEdge)
              meshEdgeToGraphVertex.emplace(*vertex.parentEdge, gv);
            if (vertex.t &&
                (*vertex.t <= m_diagnosticTolerance ||
                 Real(1) - *vertex.t <= m_diagnosticTolerance))
            {
              report.pathologicalCutCount++;
            }
          }

          if (graphVertexToOutput[gv] != InvalidIndex)
            report.graphVertexToOutputVertex.emplace(gv, graphVertexToOutput[gv]);
        }

        for (Index edge = 0; edge < static_cast<Index>(graph.parentEdgeToVertex.size()); ++edge)
        {
          const Index gv = graph.parentEdgeToVertex[static_cast<size_t>(edge)];
          if (gv != InterfaceGraph::InvalidIndex)
            meshEdgeToGraphVertex[edge] = gv;
        }

        std::vector<std::array<Index, 3>> outputCells;
        std::vector<Optional<Attribute>> outputCellAttributes;

        auto signedArea2 = [&](Index a, Index b, Index c)
        {
          const auto& x0 = outputVertices[a];
          const auto& x1 = outputVertices[b];
          const auto& x2 = outputVertices[c];
          return (x1[0] - x0[0]) * (x2[1] - x0[1]) -
                 (x1[1] - x0[1]) * (x2[0] - x0[0]);
        };

        auto parentBarycentric = [&](Index parentCell, Index outputVertex)
        {
          std::array<Real, 3> bary = {{
            Real(0), Real(0), Real(0) }};
          if (parentCell == InvalidIndex
              || parentCell >= static_cast<Index>(mesh.getCellCount()))
            return bary;

          const auto& cell = conn.getPolytope(2, parentCell);
          const auto x0 = mesh.getVertexCoordinates(cell(0));
          const auto x1 = mesh.getVertexCoordinates(cell(1));
          const auto x2 = mesh.getVertexCoordinates(cell(2));
          const auto& x = outputVertices[outputVertex];

          auto orient = [](
              const Math::SpatialPoint& a,
              const Math::SpatialPoint& b,
              const Math::SpatialPoint& c)
          {
            return (b[0] - a[0]) * (c[1] - a[1])
                 - (b[1] - a[1]) * (c[0] - a[0]);
          };

          const Real denominator = orient(x0, x1, x2);
          if (std::abs(denominator) <= Real(1e-30))
            return bary;

          bary[0] = orient(x,  x1, x2) / denominator;
          bary[1] = orient(x0, x,  x2) / denominator;
          bary[2] = orient(x0, x1, x ) / denominator;
          return bary;
        };

        auto edgeLength = [&](Index a, Index b)
        {
          return (outputVertices[a] - outputVertices[b]).norm();
        };

        auto triangleQuality = [&](Index a, Index b, Index c)
        {
          if (m_qualityMeasure)
            return m_qualityMeasure(
                outputVertices[a], outputVertices[b], outputVertices[c]);
          const Real area = std::abs(signedArea2(a, b, c)) / Real(2);
          const Real l0 = edgeLength(a, b);
          const Real l1 = edgeLength(b, c);
          const Real l2 = edgeLength(c, a);
          const Real denom = l0 * l0 + l1 * l1 + l2 * l2;
          if (denom <= Real(0))
            return Real(0);
          return Real(4) * std::sqrt(Real(3)) * area / denom;
        };

        auto cellAttribute = [&](Index parentCell, LevelSetSide side)
          -> Optional<Attribute>
        {
          if (side == LevelSetSide::Negative && m_negativeCellAttribute)
            return m_negativeCellAttribute;
          if (side == LevelSetSide::Positive && m_positiveCellAttribute)
            return m_positiveCellAttribute;
          return mesh.getCellAttribute(parentCell);
        };

        auto addTriangle = [&](Index a, Index b, Index c,
                               Index parentCell, LevelSetSide side)
        {
          if (a == b || b == c || c == a)
          {
            report.pathologicalCutCount++;
            return;
          }

          if (signedArea2(a, b, c) < Real(0))
            std::swap(b, c);

          const Real area = std::abs(signedArea2(a, b, c)) / Real(2);
          const Real quality = triangleQuality(a, b, c);
          // Never drop a produced element: a small/poor cell is reported but
          // always emitted, so the output mesh stays watertight. Degenerate
          // children are prevented upstream by the conforming snap, not by
          // discarding them here.
          if (area <= m_areaTolerance)
            report.degenerateCellCount++;
          if (quality <= m_qualityTolerance)
            report.pathologicalCutCount++;

          if (outputCells.empty())
          {
            report.minOutputCellArea = area;
            report.minOutputCellQuality = quality;
          }
          else
          {
            report.minOutputCellArea = std::min(report.minOutputCellArea, area);
            report.minOutputCellQuality =
              std::min(report.minOutputCellQuality, quality);
          }

          outputCells.push_back({{ a, b, c }});
          outputCellAttributes.push_back(cellAttribute(parentCell, side));
          report.cellProvenance.push_back({ parentCell, side });

          OutputCellReference reference;
          reference.parentCell = parentCell;
          reference.vertexBarycentric[0] = parentBarycentric(parentCell, a);
          reference.vertexBarycentric[1] = parentBarycentric(parentCell, b);
          reference.vertexBarycentric[2] = parentBarycentric(parentCell, c);
          report.cellReferences.push_back(std::move(reference));
        };

        auto findCellEdge = [&](const Polytope::Key& cell, const IndexVector& edges,
                                std::uint8_t i, std::uint8_t j) -> Index
        {
          const Index vi = cell(i);
          const Index vj = cell(j);
          for (Index edge : edges)
          {
            const auto& e = conn.getPolytope(1, edge);
            if ((e(0) == vi && e(1) == vj) || (e(0) == vj && e(1) == vi))
              return edge;
          }
          return InvalidIndex;
        };

        auto pushUnique = [](std::vector<Index>& polygon, Index value)
        {
          if (value == InvalidIndex)
            return;
          if (polygon.empty() || polygon.back() != value)
            polygon.push_back(value);
        };

        auto inside = [](LevelSetSign sign, LevelSetSide side)
        {
          if (side == LevelSetSide::Negative)
            return sign != LevelSetSign::Positive;
          if (side == LevelSetSide::Positive)
            return sign != LevelSetSign::Negative;
          return false;
        };

        auto crossingVertex =
          [&](const Polytope::Key& cell, const IndexVector& cellEdges,
              const std::array<LevelSetSign, 3>& signs,
              std::uint8_t i, std::uint8_t j) -> Index
        {
          if (signs[i] == LevelSetSign::Zero)
            return report.originalVertexToOutputVertex[cell(i)];
          if (signs[j] == LevelSetSign::Zero)
            return report.originalVertexToOutputVertex[cell(j)];

          const Index parentEdge = findCellEdge(cell, cellEdges, i, j);
          const auto it = meshEdgeToGraphVertex.find(parentEdge);
          if (it == meshEdgeToGraphVertex.end())
          {
            report.pathologicalCutCount++;
            return InvalidIndex;
          }
          return graphVertexToOutput[it->second];
        };

        auto polygonArea2 = [&](const std::vector<Index>& polygon)
        {
          Real res = 0;
          for (Index i = 0; i < polygon.size(); ++i)
          {
            const auto& a = outputVertices[polygon[i]];
            const auto& b = outputVertices[polygon[(i + 1) % polygon.size()]];
            res += a[0] * b[1] - a[1] * b[0];
          }
          return res;
        };

        auto buildSidePolygon =
          [&](const Polytope::Key& cell, const IndexVector& cellEdges,
              const std::array<LevelSetSign, 3>& signs,
              LevelSetSide side) -> std::vector<Index>
        {
          std::vector<Index> polygon;
          polygon.reserve(4);

          for (std::uint8_t k = 0; k < 3; ++k)
          {
            const std::uint8_t next = static_cast<std::uint8_t>((k + 1) % 3);
            const bool in0 = inside(signs[k], side);
            const bool in1 = inside(signs[next], side);

            if (in0)
              pushUnique(polygon, report.originalVertexToOutputVertex[cell(k)]);
            if (in0 != in1)
              pushUnique(polygon, crossingVertex(cell, cellEdges, signs, k, next));
          }

          if (polygon.size() > 1 && polygon.front() == polygon.back())
            polygon.pop_back();

          if (polygon.size() < 3)
            return {};

          if (polygonArea2(polygon) < Real(0))
            std::reverse(polygon.begin(), polygon.end());
          return polygon;
        };

        // Triangulates a side polygon, choosing the better quad diagonal.
        // Pure (no report mutation): improved is set when the alternate
        // diagonal wins, so the commit path can update the diagnostic.
        auto sidePolygonTriangles =
          [&](const std::vector<Index>& polygon, bool& improved)
          -> std::vector<std::array<Index, 3>>
        {
          improved = false;
          std::vector<std::array<Index, 3>> tris;
          if (polygon.size() < 3)
            return tris;
          if (polygon.size() == 4)
          {
            const Real firstDiagonalQuality = std::min(
                triangleQuality(polygon[0], polygon[1], polygon[2]),
                triangleQuality(polygon[0], polygon[2], polygon[3]));
            const Real secondDiagonalQuality = std::min(
                triangleQuality(polygon[1], polygon[2], polygon[3]),
                triangleQuality(polygon[1], polygon[3], polygon[0]));
            if (secondDiagonalQuality > firstDiagonalQuality)
            {
              improved = true;
              tris.push_back({{ polygon[1], polygon[2], polygon[3] }});
              tris.push_back({{ polygon[1], polygon[3], polygon[0] }});
              return tris;
            }
          }
          for (Index i = 1; i + 1 < polygon.size(); ++i)
            tris.push_back({{ polygon[0], polygon[i], polygon[i + 1] }});
          return tris;
        };

        auto trianglesMinQuality =
          [&](const std::vector<std::array<Index, 3>>& tris) -> Real
        {
          if (tris.empty())
            return Real(0);
          Real q = std::numeric_limits<Real>::infinity();
          for (const auto& t : tris)
            q = std::min(q, triangleQuality(t[0], t[1], t[2]));
          return q;
        };

        auto addSidePolygon =
          [&](const Polytope::Key& cell, const IndexVector& cellEdges,
              const std::array<LevelSetSign, 3>& signs,
              Index parentCell, LevelSetSide side)
        {
          const auto polygon = buildSidePolygon(cell, cellEdges, signs, side);
          bool improved = false;
          const auto tris = sidePolygonTriangles(polygon, improved);
          if (improved)
            report.improvedPolygonTriangulationCount++;
          for (const auto& t : tris)
            addTriangle(t[0], t[1], t[2], parentCell, side);
        };

        // Worst child quality the cut of this cell would produce, over both
        // sides. Used by the no-cut quality fallback.
        auto cutMinQuality =
          [&](const Polytope::Key& cell, const IndexVector& cellEdges,
              const std::array<LevelSetSign, 3>& signs) -> Real
        {
          bool dummy = false;
          const Real qn = trianglesMinQuality(sidePolygonTriangles(
              buildSidePolygon(cell, cellEdges, signs, LevelSetSide::Negative),
              dummy));
          const Real qp = trianglesMinQuality(sidePolygonTriangles(
              buildSidePolygon(cell, cellEdges, signs, LevelSetSide::Positive),
              dummy));
          return std::min(qn, qp);
        };

        // A cell may only be kept whole (no-cut fallback) if NONE of its edges
        // carries a real interior crossing vertex: otherwise the neighbouring
        // cut cell would reference a vertex this whole cell does not, creating
        // a hanging node (non-conforming mesh). Crossings that were snapped to
        // an endpoint are not "real" (shared, conforming) so they do not block
        // keeping the cell whole. The conforming knob for avoiding bad cuts is
        // therefore setCrossingSnapTolerance, not the per-cell fallback.
        auto cellHasRealCrossing =
          [&](const Polytope::Key& cell, const IndexVector& cellEdges) -> bool
        {
          static constexpr std::array<std::array<std::uint8_t, 2>, 3> E = {{
            {{0, 1}}, {{1, 2}}, {{2, 0}} }};
          for (const auto& pr : E)
          {
            const Index pe = findCellEdge(cell, cellEdges, pr[0], pr[1]);
            const auto it = meshEdgeToGraphVertex.find(pe);
            if (it == meshEdgeToGraphVertex.end())
              continue;
            if (graph.vertices[it->second].kind
                == InterfaceVertexKind::EdgeCut)
              return true;
          }
          return false;
        };

        std::unordered_set<Index> suppressedCells;

        for (Index c = 0; c < mesh.getCellCount(); ++c)
        {
          const auto& cell = conn.getPolytope(2, c);
          const auto& cellEdges = conn.getIncidence({ 2, 1 }, c);

          if (conn.getGeometry(2, c) != Polytope::Type::Triangle)
            Alert::MemberFunctionException(*this, __func__)
              << "Expected a triangular mesh."
              << Alert::Raise;

          std::array<LevelSetSign, 3> signs;
          Index negative = 0;
          Index zero = 0;
          Index positive = 0;
          bool invalid = false;
          for (std::uint8_t i = 0; i < 3; ++i)
          {
            signs[i] = classify(phi[cell(i)]);
            switch (signs[i])
            {
              case LevelSetSign::Negative:
                negative++;
                break;
              case LevelSetSign::Zero:
                zero++;
                break;
              case LevelSetSign::Positive:
                positive++;
                break;
              case LevelSetSign::Invalid:
                invalid = true;
                break;
            }
          }

          if (invalid || zero == 3)
          {
            report.degenerateCellCount++;
            addTriangle(
              report.originalVertexToOutputVertex[cell(0)],
              report.originalVertexToOutputVertex[cell(1)],
              report.originalVertexToOutputVertex[cell(2)],
              c,
              LevelSetSide::Degenerate);
            continue;
          }

          // No-cut quality fallback: if splitting this genuinely-crossed cell
          // would produce a child below the requested quality, keep it whole
          // on its dominant sign side. The interface is not body-fitted here;
          // a later TargetMatrixOptimization fit recovers it.
          if (m_minCutQuality > Real(0) && negative > 0 && positive > 0
              && !cellHasRealCrossing(cell, cellEdges)
              && cutMinQuality(cell, cellEdges, signs) < m_minCutQuality)
          {
            const LevelSetSide dominant = (negative >= positive)
              ? LevelSetSide::Negative : LevelSetSide::Positive;
            addTriangle(
              report.originalVertexToOutputVertex[cell(0)],
              report.originalVertexToOutputVertex[cell(1)],
              report.originalVertexToOutputVertex[cell(2)],
              c,
              dominant);
            report.uncutCellCount++;
            suppressedCells.insert(c);
            continue;
          }

          if (negative > 0 || zero > 0)
            addSidePolygon(cell, cellEdges, signs, c, LevelSetSide::Negative);
          if (positive > 0 || zero > 0)
            addSidePolygon(cell, cellEdges, signs, c, LevelSetSide::Positive);
        }

        LocalMesh::Builder builder;
        builder
          .initialize(2)
          .nodes(outputVertices.size())
          .reserve(1, graph.edges.size())
          .reserve(2, outputCells.size());

        for (const auto& vertex : outputVertices)
          builder.vertex(vertex);

        for (Index e = 0; e < graph.edges.size(); ++e)
        {
          const auto& edge = graph.edges[e];
          builder
            .polytope(Polytope::Type::Segment, {
              graphVertexToOutput[edge.v0],
              graphVertexToOutput[edge.v1] })
            .attribute({ 1, e }, edge.interfaceAttribute);

          const Real length =
            edgeLength(graphVertexToOutput[edge.v0], graphVertexToOutput[edge.v1]);
          if (length <= m_diagnosticTolerance)
            report.pathologicalCutCount++;
        }

        for (Index c = 0; c < outputCells.size(); ++c)
        {
          const auto& cell = outputCells[c];
          builder
            .polytope(Polytope::Type::Triangle, { cell[0], cell[1], cell[2] })
            .attribute({ 2, c }, outputCellAttributes[c]);
        }

        result.mesh = builder.finalize();
        result.mesh.getConnectivity().compute(2, 1);
        result.mesh.getConnectivity().compute(1, 0);
        result.mesh.getConnectivity().compute(1, 2);

        const auto& outputConn = result.mesh.getConnectivity();
        auto findOutputEdge = [&](Index a, Index b)
        {
          for (Index e = 0; e < outputConn.getCount(1); ++e)
          {
            const auto& edge = outputConn.getPolytope(1, e);
            if ((edge(0) == a && edge(1) == b) ||
                (edge(0) == b && edge(1) == a))
              return e;
          }
          return InvalidIndex;
        };

        std::unordered_set<Index> interfaceOutputEdges;
        for (Index e = 0; e < graph.edges.size(); ++e)
        {
          const auto& edge = graph.edges[e];
          const Index outEdge = findOutputEdge(
            graphVertexToOutput[edge.v0],
            graphVertexToOutput[edge.v1]);
          if (outEdge == InvalidIndex)
          {
            // Expected when the cell that would carry this segment was left
            // whole by the no-cut quality fallback; not a pathology then.
            const bool fromSuppressed = !edge.provenance.empty()
              && suppressedCells.count(edge.provenance.front().parentCell) > 0;
            if (!fromSuppressed)
              report.pathologicalCutCount++;
            continue;
          }

          result.mesh.setAttribute({ 1, outEdge }, edge.interfaceAttribute);
          interfaceOutputEdges.insert(outEdge);

          OutputInterfaceEdgeProvenance provenance;
          provenance.sourceInterfaceGraphEdge = e;
          if (!edge.provenance.empty())
            provenance.parentCell = edge.provenance.front().parentCell;
          report.interfaceEdgeProvenance.emplace(outEdge, std::move(provenance));
        }

        for (Index edge = 0; edge < conn.getCount(1); ++edge)
        {
          const auto attr = mesh.getAttribute(1, edge);
          if (!attr)
            continue;

          const auto& oldEdge = conn.getPolytope(1, edge);
          const Index out0 = report.originalVertexToOutputVertex.at(oldEdge(0));
          const Index out1 = report.originalVertexToOutputVertex.at(oldEdge(1));
          const Index outEdge = findOutputEdge(out0, out1);
          if (outEdge != InvalidIndex &&
              interfaceOutputEdges.find(outEdge) == interfaceOutputEdges.end())
          {
            result.mesh.setAttribute({ 1, outEdge }, attr);
          }
        }

        return result;
      }

    private:
      LevelSetSign classify(Real value) const
      {
        if (!std::isfinite(value))
          return LevelSetSign::Invalid;
        if (value < -m_signTolerance)
          return LevelSetSign::Negative;
        if (value > m_signTolerance)
          return LevelSetSign::Positive;
        return LevelSetSign::Zero;
      }

      std::reference_wrapper<const GridFunctionType> m_phi;
      Real m_signTolerance = 1e-12;
      Optional<Attribute> m_interfaceAttribute;
      Optional<Attribute> m_negativeCellAttribute;
      Optional<Attribute> m_positiveCellAttribute;
      Real m_diagnosticTolerance = 1e-10;
      Real m_crossingSnapTolerance = 0;
      Real m_areaTolerance = 1e-14;
      Real m_qualityTolerance = 1e-8;
      QualityMeasure m_qualityMeasure;
      Real m_minCutQuality = 0;
  };

  template <class Mesh, class Data>
  LevelSetDiscretizerTriangles(
    const Variational::GridFunction<Variational::P1<Real, Mesh>, Data>&)
    -> LevelSetDiscretizerTriangles<
         Variational::GridFunction<Variational::P1<Real, Mesh>, Data>>;
}

#endif
