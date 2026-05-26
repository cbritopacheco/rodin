/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_LEVELSETGRAPHCUTCLASSIFIER_H
#define RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_LEVELSETGRAPHCUTCLASSIFIER_H

/**
 * @file
 * @brief Variational cell classifier using seeded minimum s-t graph cut.
 *
 * Converts a continuous level-set phi into discrete cell attributes
 * chi_K in {Inside, Outside} by solving:
 *
 *   min_{chi_K} sum_K [chi_K D_K^In + (1-chi_K) D_K^Out]
 *             + lambda sum_{F=K|L} sigma_F |chi_K - chi_L|
 *
 * with quadrature-based data costs:
 *   D_K^In  = int_K max(phi, 0)^2 dx
 *   D_K^Out = int_K max(-phi, 0)^2 dx
 *
 * and attraction-weighted perimeter:
 *   sigma_F = |F| * (epsilon + |phi_bar_F| / h_F)
 *
 * Seeds force cells far from the interface (|phi_K| > seedFactor * h_K).
 * Min-cut is solved with Edmonds-Karp max-flow (BFS augmenting paths).
 * Interface facets are marked from attribute jumps: F = K+|K- where
 * chi_K+ != chi_K-.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <queue>
#include <vector>

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/Math.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/Types.h"

namespace Rodin::Adaptation::TargetMatrixOptimization
{
  /// Diagnostics returned by LevelSetGraphCutClassifier::classify().
  struct ClassificationResult
  {
    size_t insideCells               = 0;
    size_t outsideCells              = 0;
    size_t interfaceFacets           = 0;
    size_t forcedInsideSeeds         = 0;
    size_t forcedOutsideSeeds        = 0;
    size_t optimizedCells            = 0;
    Real   dataEnergy                = 0;
    Real   perimeterEnergy           = 0;
    Real   totalEnergy               = 0;
    size_t connectedInsideComponents = 0;
    size_t connectedOutsideComponents = 0;
    bool   hasNonmanifoldInterface   = false;
    Real   interfacePhiMax           = 0;
    Real   interfacePhiRms           = 0;
  };

  /**
   * @brief Classifies mesh cells as Inside/Outside via seeded minimum s-t cut.
   *
   * The classifier is stateless; call classify() on any LocalMesh. Connectivity
   * is computed automatically if not already present. All attribute
   * assignments are written directly to the mesh.
   */
  class LevelSetGraphCutClassifier
  {
  public:
    struct Options
    {
      Geometry::Attribute insideAttribute    = 1;
      Geometry::Attribute outsideAttribute   = 10;
      Geometry::Attribute interfaceAttribute = 30;

      /// Regularization weight lambda on the perimeter term.
      Real perimeterWeight        = Real(1);

      /// epsilon in sigma_F = |F| * (epsilon + |phi_bar_F| / h_F).
      /// Prevents zero edge costs far from the interface.
      Real attractionEpsilon      = Real(1e-3);

      /// Cells with |phi_K| > seedDistanceFactor * h_K are hard-seeded.
      Real seedDistanceFactor     = Real(2);

      /// Terminal capacity for hard seeds (acts as infinity).
      Real hardSeedCost           = Real(1e30);

      /// Quadrature order for data costs and facet weights.
      size_t quadratureOrder      = 4;

      /// If false, skip interface facet marking (attributes only).
      bool markInterfaceFacets    = true;
    };

    /**
     * @brief Classify cells using default options.
     *
     * Delegates to the two-parameter overload with a default-constructed Options.
     * This overload exists because Clang does not allow a default argument
     * `opts = {}` for a nested-struct with default member initializers when
     * the declaration appears inside the enclosing class body (CWG 1683).
     *
     * @tparam PhiFn  Callable: Real phi(const Math::SpatialPoint&, Real t).
     * @param mesh    Mesh to classify (modified in-place).
     * @param phi     Level-set function. phi < 0 = Inside, phi > 0 = Outside.
     * @param time    Current time passed to phi.
     * @return        Diagnostics.
     */
    template <class PhiFn>
    ClassificationResult classify(
        Geometry::LocalMesh& mesh,
        const PhiFn&         phi,
        Real                 time) const
    {
      return classify(mesh, phi, time, Options{});
    }

    /**
     * @brief Classify cells and mark interface facets.
     *
     * @tparam PhiFn  Callable: Real phi(const Math::SpatialPoint&, Real t).
     * @param mesh    Mesh to classify (modified in-place).
     * @param phi     Level-set function. phi < 0 = Inside, phi > 0 = Outside.
     * @param time    Current time passed to phi.
     * @param opts    Classifier options.
     * @return        Diagnostics.
     */
    template <class PhiFn>
    ClassificationResult classify(
        Geometry::LocalMesh&  mesh,
        const PhiFn&          phi,
        Real                  time,
        const Options&        opts) const
    {
      // ── Connectivity ─────────────────────────────────────────────────────
      auto& conn = mesh.getConnectivity();
      conn.compute(2, 1);   // cell  → edges
      conn.compute(1, 2);   // edge  → cells
      conn.compute(1, 0);   // edge  → vertices

      const size_t nCells = mesh.getCellCount();
      const size_t nEdges = conn.getCount(1);
      const size_t nVerts = mesh.getVertexCount();

      // ── Build s-t graph ──────────────────────────────────────────────────
      // Nodes: 0=source (Inside), 1=sink (Outside), K+2=cell K
      const int source     = 0;
      const int sink       = 1;
      const int nodeOffset = 2;

      FlowGraph graph(nCells + 2);

      // Per-cell storage for energy post-processing
      std::vector<Real> cellArea(nCells,    Real(0));
      std::vector<Real> costIn(nCells,      Real(0));
      std::vector<Real> costOut(nCells,     Real(0));
      std::vector<Real> phiRep(nCells,      Real(0));

      const auto& qfTri = QF::PolytopeQuadratureFormula::get(
          opts.quadratureOrder, Geometry::Polytope::Type::Triangle);

      size_t forcedInside  = 0;
      size_t forcedOutside = 0;
      size_t optimized     = 0;

      // ── Terminal edges (data costs / seeds) ───────────────────────────────
      for (Index c = 0; c < static_cast<Index>(nCells); ++c)
      {
        if (conn.getGeometry(2, c) != Geometry::Polytope::Type::Triangle)
          continue;

        auto cellIt = mesh.getPolytope(2, c);

        Real area = Real(0), phiW = Real(0);
        Real cIn  = Real(0), cOut = Real(0);

        for (size_t q = 0; q < qfTri.getSize(); ++q)
        {
          const auto& refPt = qfTri.getPoint(q);
          Math::SpatialMatrix<Real> J;
          cellIt->getTransformation().jacobian(J, refPt);
          const Real detJ = std::abs(J.determinant());
          const Real w    = qfTri.getWeight(q) * detJ;

          const auto x  = Geometry::Point(*cellIt, refPt).getPhysicalCoordinates();
          const Real pv = phi(x, time);

          area += w;
          phiW += w * pv;

          const Real posv = std::max(pv, Real(0));
          const Real negv = std::max(-pv, Real(0));
          cIn  += w * posv * posv;
          cOut += w * negv * negv;
        }

        cellArea[c] = area;
        costIn[c]   = cIn;
        costOut[c]  = cOut;
        phiRep[c]   = (area > Real(0)) ? phiW / area : Real(0);

        const Real hK  = std::sqrt(std::max(area, Real(1e-30)));
        const Real tau = opts.seedDistanceFactor * hK;

        Real sCap, tCap; // s→K, K→t
        if (phiRep[c] < -tau)
        {
          // Definitely inside: make labeling Outside impossible
          sCap = opts.hardSeedCost;
          tCap = Real(0);
          ++forcedInside;
        }
        else if (phiRep[c] > tau)
        {
          // Definitely outside: make labeling Inside impossible
          sCap = Real(0);
          tCap = opts.hardSeedCost;
          ++forcedOutside;
        }
        else
        {
          sCap = cOut;
          tCap = cIn;
          ++optimized;
        }

        addEdge(graph, source,       c + nodeOffset, sCap);
        addEdge(graph, c + nodeOffset, sink,          tCap);
      }

      // ── Neighbor edges (perimeter / attraction term) ──────────────────────
      const auto& qfSeg = QF::PolytopeQuadratureFormula::get(
          opts.quadratureOrder, Geometry::Polytope::Type::Segment);

      std::vector<Real> edgeSigma(nEdges, Real(0));

      for (Index e = 0; e < static_cast<Index>(nEdges); ++e)
      {
        const auto& adj = conn.getIncidence({1, 2}, e);
        if (static_cast<size_t>(adj.size()) != 2) continue;

        const Index cK = adj[0];
        const Index cL = adj[1];
        if (conn.getGeometry(2, cK) != Geometry::Polytope::Type::Triangle) continue;
        if (conn.getGeometry(2, cL) != Geometry::Polytope::Type::Triangle) continue;

        auto edgeIt = mesh.getPolytope(1, e);
        Real edgeLen = Real(0), phiAvg = Real(0), phiW = Real(0);

        for (size_t q = 0; q < qfSeg.getSize(); ++q)
        {
          const auto& refPt = qfSeg.getPoint(q);
          Math::SpatialMatrix<Real> J;
          edgeIt->getTransformation().jacobian(J, refPt);
          const Real len = std::hypot(J(0, 0), J(1, 0));
          const Real w   = qfSeg.getWeight(q) * len;

          const auto x  = Geometry::Point(*edgeIt, refPt).getPhysicalCoordinates();
          const Real pv = phi(x, time);

          edgeLen += w;
          phiAvg  += w * pv;
          phiW    += w;
        }
        if (phiW > Real(0)) phiAvg /= phiW;

        // h_F = average sqrt(adjacent cell area)
        const Real hK = std::sqrt(std::max(cellArea[cK], Real(1e-30)));
        const Real hL = std::sqrt(std::max(cellArea[cL], Real(1e-30)));
        const Real hF = Real(0.5) * (hK + hL);

        const Real sigmaF = edgeLen *
            (opts.attractionEpsilon + std::abs(phiAvg) / (hF + Real(1e-14)));
        edgeSigma[e] = sigmaF;

        const Real w = opts.perimeterWeight * sigmaF;
        // Undirected: two directed edges each with capacity w
        addEdge(graph, cK + nodeOffset, cL + nodeOffset, w);
        addEdge(graph, cL + nodeOffset, cK + nodeOffset, w);
      }

      // ── Max-flow → min-cut ────────────────────────────────────────────────
      runMaxFlow(graph, source, sink);
      const std::vector<bool> srcSide = sourceReachable(graph, source);

      // ── Assign cell attributes ────────────────────────────────────────────
      ClassificationResult result;
      result.forcedInsideSeeds  = forcedInside;
      result.forcedOutsideSeeds = forcedOutside;
      result.optimizedCells     = optimized;

      std::vector<bool> insideMask(nCells, false);
      for (Index c = 0; c < static_cast<Index>(nCells); ++c)
      {
        const bool inside = srcSide[static_cast<size_t>(c) + nodeOffset];
        insideMask[c] = inside;
        mesh.setAttribute({2, c}, inside ? opts.insideAttribute : opts.outsideAttribute);
        if (inside) ++result.insideCells;
        else        ++result.outsideCells;
      }

      // ── Data energy ───────────────────────────────────────────────────────
      Real dataEnergy = Real(0);
      for (Index c = 0; c < static_cast<Index>(nCells); ++c)
        dataEnergy += insideMask[c] ? costIn[c] : costOut[c];
      result.dataEnergy = dataEnergy;

      // ── Clear old interface, mark new interface facets, compute energies ──
      // First clear all interface edge attributes
      for (Index e = 0; e < static_cast<Index>(nEdges); ++e)
      {
        const auto attr = mesh.getAttribute(1, e);
        if (attr && *attr == opts.interfaceAttribute)
          mesh.setAttribute({1, e}, Optional<Geometry::Attribute>{});
      }

      Real perimeterEnergy = Real(0);
      for (Index e = 0; e < static_cast<Index>(nEdges); ++e)
      {
        const auto& adj = conn.getIncidence({1, 2}, e);
        if (static_cast<size_t>(adj.size()) != 2) continue;

        const auto a0 = mesh.getAttribute(2, adj[0]);
        const auto a1 = mesh.getAttribute(2, adj[1]);
        if (!a0 || !a1) continue;

        if (*a0 != *a1)
        {
          if (opts.markInterfaceFacets)
            mesh.setAttribute({1, e}, opts.interfaceAttribute);
          ++result.interfaceFacets;
          perimeterEnergy += opts.perimeterWeight * edgeSigma[e];
        }
      }
      result.perimeterEnergy = perimeterEnergy;
      result.totalEnergy     = dataEnergy + perimeterEnergy;

      // ── Connected components (BFS per label) ─────────────────────────────
      {
        // Build cell–cell adjacency through shared edges
        std::vector<std::vector<Index>> nbrs(nCells);
        for (Index e = 0; e < static_cast<Index>(nEdges); ++e)
        {
          const auto& adj = conn.getIncidence({1, 2}, e);
          if (static_cast<size_t>(adj.size()) != 2) continue;
          nbrs[adj[0]].push_back(adj[1]);
          nbrs[adj[1]].push_back(adj[0]);
        }

        std::vector<bool> visited(nCells, false);
        size_t insideComp  = 0;
        size_t outsideComp = 0;

        for (Index c = 0; c < static_cast<Index>(nCells); ++c)
        {
          if (visited[c]) continue;
          visited[c] = true;
          const bool label = insideMask[c];

          std::queue<Index> bfsQ;
          bfsQ.push(c);
          while (!bfsQ.empty())
          {
            const Index cur = bfsQ.front(); bfsQ.pop();
            for (Index nb : nbrs[cur])
            {
              if (!visited[nb] && insideMask[nb] == label)
              {
                visited[nb] = true;
                bfsQ.push(nb);
              }
            }
          }

          if (label) ++insideComp;
          else       ++outsideComp;
        }

        result.connectedInsideComponents  = insideComp;
        result.connectedOutsideComponents = outsideComp;
      }

      // ── Nonmanifold interface check ───────────────────────────────────────
      {
        std::vector<Index> ifaceDeg(nVerts, 0);
        for (Index e = 0; e < static_cast<Index>(nEdges); ++e)
        {
          const auto attr = mesh.getAttribute(1, e);
          if (!attr || *attr != opts.interfaceAttribute) continue;
          const auto& key = conn.getPolytope(1, e);
          ++ifaceDeg[key(0)];
          ++ifaceDeg[key(1)];
        }
        for (Index v = 0; v < static_cast<Index>(nVerts); ++v)
        {
          if (ifaceDeg[v] > 2)
          {
            result.hasNonmanifoldInterface = true;
            break;
          }
        }
      }

      // ── Interface phi diagnostics ─────────────────────────────────────────
      {
        Real phiSqSum  = Real(0);
        Real phiMaxAbs = Real(0);
        size_t cnt     = 0;

        for (Index e = 0; e < static_cast<Index>(nEdges); ++e)
        {
          const auto attr = mesh.getAttribute(1, e);
          if (!attr || *attr != opts.interfaceAttribute) continue;

          const auto& key = conn.getPolytope(1, e);
          for (int vi = 0; vi < 2; ++vi)
          {
            const auto x  = mesh.getVertexCoordinates(key(vi));
            const Real pv = phi(x, time);
            phiSqSum  += pv * pv;
            phiMaxAbs  = std::max(phiMaxAbs, std::abs(pv));
            ++cnt;
          }
        }
        result.interfacePhiRms = (cnt > 0)
            ? std::sqrt(phiSqSum / static_cast<Real>(cnt))
            : Real(0);
        result.interfacePhiMax = phiMaxAbs;
      }

      return result;
    }

  private:
    // ── Max-flow internals (Edmonds-Karp BFS augmenting paths) ───────────────

    struct FlowEdge
    {
      int  to;
      Real cap;
      int  rev; // index of reverse edge in graph[to]
    };

    using FlowGraph = std::vector<std::vector<FlowEdge>>;

    static void addEdge(FlowGraph& g, int u, int v, Real cap)
    {
      g[u].push_back({v, cap, static_cast<int>(g[v].size())});
      g[v].push_back({u, Real(0), static_cast<int>(g[u].size()) - 1});
    }

    static void runMaxFlow(FlowGraph& g, int source, int sink)
    {
      const int n = static_cast<int>(g.size());
      for (;;)
      {
        std::vector<int> parent(n, -1);
        std::vector<int> parentEdge(n, -1);
        parent[source] = source;

        std::queue<int> bfsQ;
        bfsQ.push(source);
        bool found = false;

        while (!bfsQ.empty() && !found)
        {
          const int v = bfsQ.front(); bfsQ.pop();
          for (int i = 0; i < static_cast<int>(g[v].size()) && !found; ++i)
          {
            const FlowEdge& e = g[v][i];
            if (parent[e.to] == -1 && e.cap > Real(1e-15))
            {
              parent[e.to]     = v;
              parentEdge[e.to] = i;
              if (e.to == sink)
                found = true;
              else
                bfsQ.push(e.to);
            }
          }
        }

        if (!found) break;

        // Bottleneck
        Real bottleneck = std::numeric_limits<Real>::max();
        for (int v = sink; v != source; )
        {
          const int u  = parent[v];
          const int ei = parentEdge[v];
          bottleneck = std::min(bottleneck, g[u][ei].cap);
          v = u;
        }
        if (bottleneck < Real(1e-15)) break;

        // Augment
        for (int v = sink; v != source; )
        {
          const int u  = parent[v];
          const int ei = parentEdge[v];
          g[u][ei].cap -= bottleneck;
          g[g[u][ei].to][g[u][ei].rev].cap += bottleneck;
          v = u;
        }
      }
    }

    static std::vector<bool> sourceReachable(const FlowGraph& g, int source)
    {
      const int n = static_cast<int>(g.size());
      std::vector<bool> visited(n, false);
      std::queue<int> q;
      q.push(source);
      visited[source] = true;
      while (!q.empty())
      {
        const int v = q.front(); q.pop();
        for (const FlowEdge& e : g[v])
        {
          if (!visited[e.to] && e.cap > Real(1e-15))
          {
            visited[e.to] = true;
            q.push(e.to);
          }
        }
      }
      return visited;
    }
  };

} // namespace Rodin::Adaptation::TargetMatrixOptimization

#endif // RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_LEVELSETGRAPHCUTCLASSIFIER_H
