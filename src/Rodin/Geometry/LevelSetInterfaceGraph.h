/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_LEVELSETINTERFACEGRAPH_H
#define RODIN_GEOMETRY_LEVELSETINTERFACEGRAPH_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <utility>
#include <vector>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Variational/P1/P1.h"

#include "Mesh.h"
#include "Polytope.h"
#include "Types.h"

namespace Rodin::Geometry
{
  enum class LevelSetSign
  {
    Negative,
    Zero,
    Positive,
    Invalid
  };

  enum class InterfaceVertexKind
  {
    OriginalVertex,
    EdgeCut
  };

  struct InterfaceVertex
  {
    Math::SpatialPoint x;
    InterfaceVertexKind kind = InterfaceVertexKind::EdgeCut;
    Optional<Index> originalVertex;
    Optional<Index> parentEdge;
    Optional<Real> t;
  };

  struct InterfaceEdgeProvenance
  {
    Index parentCell = std::numeric_limits<Index>::max();
    std::array<Index, 2> parentEdges =
      {{ std::numeric_limits<Index>::max(), std::numeric_limits<Index>::max() }};
    Optional<Attribute> parentCellAttribute;
  };

  struct InterfaceEdge
  {
    Index v0 = std::numeric_limits<Index>::max();
    Index v1 = std::numeric_limits<Index>::max();

    Optional<Attribute> interfaceAttribute;
    std::vector<InterfaceEdgeProvenance> provenance;
  };

  struct InterfaceChain
  {
    // Open chains satisfy vertices.size() == edges.size() + 1, with
    // edges[k] connecting vertices[k] to vertices[k + 1].
    // Closed chains satisfy vertices.size() == edges.size(), with
    // edges[k] connecting vertices[k] to vertices[(k + 1) % vertices.size()].
    // The orientation is deterministic, but does not yet encode inside/outside
    // or positive/negative side.
    std::vector<Index> vertices;
    std::vector<Index> edges;
    bool closed = false;
  };

  struct InterfaceGraph
  {
    static constexpr Index InvalidIndex = std::numeric_limits<Index>::max();

    std::vector<InterfaceVertex> vertices;
    std::vector<InterfaceEdge> edges;
    std::vector<InterfaceChain> chains;

    /**
     * Maps a background edge carrying a strict sign-changing crossing to the
     * InterfaceGraph vertex used for that crossing. When crossing snap is
     * enabled, the mapped graph vertex may be an original mesh vertex.
     */
    std::vector<Index> parentEdgeToVertex;

    std::vector<Index> degenerateCells;
    std::vector<Index> invalidCells;

    Index nearVertexCrossingCount = 0;
    Index snappedCrossingCount = 0;
    Real maxSnapDistance = 0;
    Real maxInterfaceDeviation = 0;
  };

  struct ProtectedSegment
  {
    Index v0 = std::numeric_limits<Index>::max();
    Index v1 = std::numeric_limits<Index>::max();
    Optional<Attribute> attribute;
    std::vector<Index> sourceInterfaceEdges;
  };

  struct PlanarStraightLineGraph
  {
    std::vector<Math::SpatialPoint> vertices;
    std::vector<ProtectedSegment> segments;
  };

  class InterfaceGraphPSLGBuilder
  {
    public:
      PlanarStraightLineGraph build(const InterfaceGraph& graph) const
      {
        PlanarStraightLineGraph res;
        res.vertices.reserve(graph.vertices.size());
        for (const auto& vertex : graph.vertices)
          res.vertices.push_back(vertex.x);

        res.segments.reserve(graph.edges.size());
        for (Index i = 0; i < graph.edges.size(); ++i)
        {
          const auto& edge = graph.edges[i];

          ProtectedSegment segment;
          segment.v0 = edge.v0;
          segment.v1 = edge.v1;
          segment.attribute = edge.interfaceAttribute;
          segment.sourceInterfaceEdges.push_back(i);

          res.segments.push_back(std::move(segment));
        }
        return res;
      }
  };

  template <class... Params>
  class LevelSetInterfaceGraph;

  template <class Mesh, class Data>
  class LevelSetInterfaceGraph<
    Variational::GridFunction<Variational::P1<Real, Mesh>, Data>>
  {
    public:
      using MeshType = Mesh;
      using FESType = Variational::P1<Real, Mesh>;
      using GridFunctionType = Variational::GridFunction<FESType, Data>;
      using GraphType = InterfaceGraph;

      static constexpr Index InvalidIndex = std::numeric_limits<Index>::max();

      explicit LevelSetInterfaceGraph(const GridFunctionType& phi)
        : m_phi(phi),
          m_signTolerance(1e-12)
      {}

      LevelSetInterfaceGraph& setSignTolerance(Real tol)
      {
        m_signTolerance = std::max(Real(0), tol);
        return *this;
      }

      LevelSetInterfaceGraph& setInterfaceAttribute(const Optional<Attribute>& attr)
      {
        m_interfaceAttribute = attr;
        return *this;
      }

      /**
       * @brief Snaps strict edge crossings close to an endpoint.
       *
       * The tolerance is measured in edge interpolation coordinates. With the
       * default value zero the extracted P1 interface is exact. For @f$\tau>0@f$,
       * crossings with @f$t \le \tau@f$ or @f$1-t \le \tau@f$ are represented by
       * the corresponding existing endpoint and reported in the graph
       * diagnostics. This deliberately trades a bounded local interface
       * deviation for eliminating near-vertex sliver topology.
       */
      LevelSetInterfaceGraph& setCrossingSnapTolerance(Real tol)
      {
        m_crossingSnapTolerance = std::max(Real(0), tol);
        return *this;
      }

      Real getSignTolerance() const
      {
        return m_signTolerance;
      }

      const Optional<Attribute>& getInterfaceAttribute() const
      {
        return m_interfaceAttribute;
      }

      Real getCrossingSnapTolerance() const
      {
        return m_crossingSnapTolerance;
      }

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

      GraphType extract() const
      {
        const auto& phi = getGridFunction();
        const auto& mesh = phi.getFiniteElementSpace().getMesh();

        if (mesh.getDimension() != 2)
          Alert::MemberFunctionException(*this, __func__)
            << "Expected a mesh of topological dimension 2."
            << Alert::Raise;

        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 2, 1);
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 1, 0);

        const auto& conn = mesh.getConnectivity();

        for (Index c = 0; c < mesh.getCellCount(); ++c)
        {
          if (conn.getGeometry(2, c) != Polytope::Type::Triangle)
            Alert::MemberFunctionException(*this, __func__)
              << "Expected a triangular mesh."
              << Alert::Raise;
        }

        GraphType graph;

        const Index nv = static_cast<Index>(mesh.getVertexCount());
        const Index ne = static_cast<Index>(conn.getCount(1));
        graph.parentEdgeToVertex.assign(ne, GraphType::InvalidIndex);

        std::vector<Optional<Index>> originalVertexToGraph(nv);
        std::vector<Optional<Index>> edgeToGraph(ne);
        UnorderedMap<
          Polytope::Key,
          Index,
          Polytope::Key::SymmetricHash,
          Polytope::Key::SymmetricEquality> graphEdgeToIndex;

        auto makeOriginalVertex = [&](Index original) -> Index
        {
          if (originalVertexToGraph[original])
            return *originalVertexToGraph[original];

          InterfaceVertex vertex;
          vertex.x = mesh.getVertexCoordinates(original);
          vertex.kind = InterfaceVertexKind::OriginalVertex;
          vertex.originalVertex = original;

          const Index idx = static_cast<Index>(graph.vertices.size());
          graph.vertices.push_back(std::move(vertex));
          originalVertexToGraph[original] = idx;
          return idx;
        };

        auto makeEdgeVertex = [&](Index edge) -> Index
        {
          const auto& e = conn.getPolytope(1, edge);
          const Index a = e(0);
          const Index b = e(1);
          const Real fa = phi[a];
          const Real fb = phi[b];
          const LevelSetSign sa = classify(fa);
          const LevelSetSign sb = classify(fb);

          if (sa == LevelSetSign::Zero)
            return makeOriginalVertex(a);
          if (sb == LevelSetSign::Zero)
            return makeOriginalVertex(b);

          const Real denom = fa - fb;
          Real t = (denom == Real(0)) ? Real(0.5) : fa / denom;

          if (edgeToGraph[edge])
            return *edgeToGraph[edge];

          t = std::max(Real(0), std::min(Real(1), t));
          const auto exactX = (Real(1) - t) * mesh.getVertexCoordinates(a)
                            + t * mesh.getVertexCoordinates(b);

          if (m_crossingSnapTolerance > Real(0) &&
              (t <= m_crossingSnapTolerance ||
               Real(1) - t <= m_crossingSnapTolerance))
          {
            graph.nearVertexCrossingCount++;
            const Index original = (t <= Real(1) - t) ? a : b;
            const auto snappedX = mesh.getVertexCoordinates(original);
            const Real distance = (exactX - snappedX).norm();
            graph.maxSnapDistance = std::max(graph.maxSnapDistance, distance);
            graph.maxInterfaceDeviation =
              std::max(graph.maxInterfaceDeviation, distance);

            const Index idx = makeOriginalVertex(original);
            edgeToGraph[edge] = idx;
            graph.parentEdgeToVertex[edge] = idx;
            graph.snappedCrossingCount++;
            return idx;
          }

          InterfaceVertex vertex;
          vertex.x = exactX;
          vertex.kind = InterfaceVertexKind::EdgeCut;
          vertex.parentEdge = edge;
          vertex.t = t;

          const Index idx = static_cast<Index>(graph.vertices.size());
          graph.vertices.push_back(std::move(vertex));
          edgeToGraph[edge] = idx;
          graph.parentEdgeToVertex[edge] = idx;
          return idx;
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

        struct Candidate
        {
          Index vertex = InvalidIndex;
          Index parentEdge = InvalidIndex;
        };

        auto addCandidate = [](std::vector<Candidate>& candidates, Candidate candidate)
        {
          for (auto& existing : candidates)
          {
            if (existing.vertex == candidate.vertex)
            {
              if (existing.parentEdge == InvalidIndex)
                existing.parentEdge = candidate.parentEdge;
              return;
            }
          }
          candidates.push_back(candidate);
        };

        auto addEdge = [&](Index a, Index b, InterfaceEdgeProvenance provenance)
        {
          if (a == b)
            return;

          const Polytope::Key key({ a, b });

          auto it = graphEdgeToIndex.find(key);
          if (it != graphEdgeToIndex.end())
          {
            auto& edge = graph.edges[it->second];
            const auto sameProvenance =
              [&provenance](const InterfaceEdgeProvenance& existing)
              {
                const bool sameParentEdges =
                  existing.parentEdges == provenance.parentEdges ||
                  (existing.parentEdges[0] == provenance.parentEdges[1] &&
                   existing.parentEdges[1] == provenance.parentEdges[0]);

                return existing.parentCell == provenance.parentCell &&
                       sameParentEdges;
              };
            if (std::find_if(edge.provenance.begin(), edge.provenance.end(),
                             sameProvenance) == edge.provenance.end())
            {
              edge.provenance.push_back(std::move(provenance));
            }
            return;
          }

          InterfaceEdge edge;
          edge.v0 = a;
          edge.v1 = b;
          edge.interfaceAttribute = m_interfaceAttribute;
          edge.provenance.push_back(std::move(provenance));
          graphEdgeToIndex.emplace(key, static_cast<Index>(graph.edges.size()));
          graph.edges.push_back(std::move(edge));
        };

        static constexpr std::array<std::array<std::uint8_t, 2>, 3> LocalEdges =
          {{{{0, 1}}, {{1, 2}}, {{2, 0}}}};

        for (Index c = 0; c < mesh.getCellCount(); ++c)
        {
          const auto& cell = conn.getPolytope(2, c);
          const auto& cellEdges = conn.getIncidence({2, 1}, c);

          std::array<LevelSetSign, 3> signs;
          size_t negativeCount = 0;
          size_t zeroCount = 0;
          size_t positiveCount = 0;
          bool invalid = false;

          for (std::uint8_t i = 0; i < 3; ++i)
          {
            signs[i] = classify(phi[cell(i)]);
            switch (signs[i])
            {
              case LevelSetSign::Negative:
                negativeCount++;
                break;
              case LevelSetSign::Zero:
                zeroCount++;
                break;
              case LevelSetSign::Positive:
                positiveCount++;
                break;
              case LevelSetSign::Invalid:
                invalid = true;
                break;
            }
          }

          if (invalid)
          {
            graph.invalidCells.push_back(c);
            continue;
          }

          if (zeroCount == 3)
          {
            graph.degenerateCells.push_back(c);
            continue;
          }

          if (negativeCount == 3 || positiveCount == 3)
            continue;

          if (zeroCount == 1 && (negativeCount == 2 || positiveCount == 2))
            continue;

          if (zeroCount == 2)
          {
            std::array<std::uint8_t, 2> zeroLocal = {{0, 0}};
            size_t pos = 0;
            for (std::uint8_t i = 0; i < 3; ++i)
            {
              if (signs[i] == LevelSetSign::Zero)
                zeroLocal[pos++] = i;
            }

            const Index parentEdge =
              findCellEdge(cell, cellEdges, zeroLocal[0], zeroLocal[1]);
            addEdge(
              makeOriginalVertex(cell(zeroLocal[0])),
              makeOriginalVertex(cell(zeroLocal[1])),
              { c, {{ parentEdge, parentEdge }}, mesh.getCellAttribute(c) });
            continue;
          }

          std::vector<Candidate> candidates;
          candidates.reserve(2);

          for (const auto& localEdge : LocalEdges)
          {
            const std::uint8_t i = localEdge[0];
            const std::uint8_t j = localEdge[1];
            const Index parentEdge = findCellEdge(cell, cellEdges, i, j);

            if (signs[i] == LevelSetSign::Zero &&
                signs[j] != LevelSetSign::Zero)
            {
              addCandidate(candidates, { makeOriginalVertex(cell(i)), InvalidIndex });
            }
            else if (signs[j] == LevelSetSign::Zero &&
                     signs[i] != LevelSetSign::Zero)
            {
              addCandidate(candidates, { makeOriginalVertex(cell(j)), InvalidIndex });
            }
            else if ((signs[i] == LevelSetSign::Negative &&
                      signs[j] == LevelSetSign::Positive) ||
                     (signs[i] == LevelSetSign::Positive &&
                      signs[j] == LevelSetSign::Negative))
            {
              addCandidate(candidates, { makeEdgeVertex(parentEdge), parentEdge });
            }
          }

          if (candidates.size() == 2)
          {
            addEdge(
              candidates[0].vertex,
              candidates[1].vertex,
              { c,
                {{ candidates[0].parentEdge, candidates[1].parentEdge }},
                mesh.getCellAttribute(c) });
          }
          else if (candidates.size() > 2)
          {
            graph.degenerateCells.push_back(c);
          }
        }

        extractChains(graph);
        return graph;
      }

    private:
      const GridFunctionType& getGridFunction() const
      {
        return m_phi.get();
      }

      void extractChains(GraphType& graph) const
      {
        std::vector<std::vector<std::array<Index, 2>>> adjacency(graph.vertices.size());
        for (Index e = 0; e < graph.edges.size(); ++e)
        {
          const auto& edge = graph.edges[e];
          adjacency[edge.v0].push_back({{ edge.v1, e }});
          adjacency[edge.v1].push_back({{ edge.v0, e }});
        }

        for (auto& adj : adjacency)
        {
          std::sort(adj.begin(), adj.end(),
            [](const auto& a, const auto& b)
            {
              if (a[0] != b[0])
                return a[0] < b[0];
              return a[1] < b[1];
            });
        }

        std::vector<bool> visited(graph.edges.size(), false);
        for (Index seed = 0; seed < graph.edges.size(); ++seed)
        {
          if (visited[seed])
            continue;

          std::vector<Index> stack{ seed };
          std::vector<Index> componentEdges;
          std::vector<Index> componentVertices;
          std::vector<bool> edgeInComponent(graph.edges.size(), false);
          std::vector<bool> vertexInComponent(graph.vertices.size(), false);

          edgeInComponent[seed] = true;
          while (!stack.empty())
          {
            const Index edgeIndex = stack.back();
            stack.pop_back();
            componentEdges.push_back(edgeIndex);

            const auto& edge = graph.edges[edgeIndex];
            for (Index v : { edge.v0, edge.v1 })
            {
              if (!vertexInComponent[v])
              {
                vertexInComponent[v] = true;
                componentVertices.push_back(v);
              }

              for (const auto& next : adjacency[v])
              {
                const Index nextEdge = next[1];
                if (!edgeInComponent[nextEdge])
                {
                  edgeInComponent[nextEdge] = true;
                  stack.push_back(nextEdge);
                }
              }
            }
          }

          for (Index e : componentEdges)
            visited[e] = true;

          std::vector<Index> endpoints;
          bool ordered = true;
          for (Index v : componentVertices)
          {
            const size_t degree = adjacency[v].size();
            if (degree == 1)
              endpoints.push_back(v);
            else if (degree != 2)
            {
              ordered = false;
              break;
            }
          }

          const bool closed = ordered && endpoints.empty();
          const bool open = ordered && endpoints.size() == 2;
          if (!closed && !open)
            continue;

          const Index startVertex = closed
            ? *std::min_element(componentVertices.begin(), componentVertices.end())
            : *std::min_element(endpoints.begin(), endpoints.end());

          InterfaceChain chain;
          chain.closed = closed;
          std::vector<bool> chainEdgeSeen(graph.edges.size(), false);
          Index currentVertex = startVertex;
          Index previousEdge = InvalidIndex;

          while (true)
          {
            chain.vertices.push_back(currentVertex);

            Index nextEdge = InvalidIndex;
            for (const auto& next : adjacency[currentVertex])
            {
              const Index candidate = next[1];
              if (!edgeInComponent[candidate] ||
                  candidate == previousEdge ||
                  chainEdgeSeen[candidate])
                continue;
              nextEdge = candidate;
              break;
            }

            if (nextEdge == InvalidIndex)
              break;

            chainEdgeSeen[nextEdge] = true;
            chain.edges.push_back(nextEdge);

            const auto& edge = graph.edges[nextEdge];
            currentVertex = edge.v0 == currentVertex ? edge.v1 : edge.v0;
            previousEdge = nextEdge;

            if (currentVertex == startVertex)
              break;
          }

          if (chain.edges.size() == componentEdges.size() &&
              (closed ? currentVertex == startVertex : currentVertex != startVertex))
          {
            graph.chains.push_back(std::move(chain));
          }
        }
      }

      std::reference_wrapper<const GridFunctionType> m_phi;
      Real m_signTolerance;
      Real m_crossingSnapTolerance = 0;
      Optional<Attribute> m_interfaceAttribute;
  };

  template <class Mesh, class Data>
  LevelSetInterfaceGraph(
    const Variational::GridFunction<Variational::P1<Real, Mesh>, Data>&)
    -> LevelSetInterfaceGraph<
         Variational::GridFunction<Variational::P1<Real, Mesh>, Data>>;
}

#endif
