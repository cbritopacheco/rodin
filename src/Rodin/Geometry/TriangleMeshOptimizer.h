/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_TRIANGLEMESHOPTIMIZER_H
#define RODIN_GEOMETRY_TRIANGLEMESHOPTIMIZER_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <functional>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Rodin/Types.h"
#include "Mesh.h"
#include "Polytope.h"

namespace Rodin::Geometry
{
  /**
   * @brief Signed twice-area of an oriented triangle (2D).
   */
  inline Real triangleOrient2(const Math::SpatialPoint& a,
                              const Math::SpatialPoint& b,
                              const Math::SpatialPoint& c)
  {
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
  }

  /**
   * @brief Default normalized shape quality @f$4\sqrt3 A/\sum \ell^2@f$.
   *
   * Any callable with this signature can be supplied as the Quality template
   * parameter of TriangleMeshOptimizer.
   */
  struct ShapeQuality
  {
    Real operator()(const Math::SpatialPoint& a,
                     const Math::SpatialPoint& b,
                     const Math::SpatialPoint& c) const
    {
      const Real ar = std::abs(triangleOrient2(a, b, c)) / Real(2);
      const Real d = (b - a).squaredNorm()
                   + (c - b).squaredNorm()
                   + (a - c).squaredNorm();
      if (d <= Real(0))
        return Real(0);
      return Real(4) * std::sqrt(Real(3)) * ar / d;
    }
  };

  struct TriangleMeshOptimizerReport
  {
    std::size_t iterations = 0;
    std::size_t splits = 0;
    std::size_t collapses = 0;
    std::size_t swaps = 0;
    std::size_t smooths = 0;
    std::size_t featureSmooths = 0;
    Real minQualityBefore = 0;
    Real minQualityAfter = 0;
  };

  struct TriangleMeshOptimizerParameters
  {
    Real hMin = 0;
    Real hMax = std::numeric_limits<Real>::infinity();
    Real minQuality = 1e-4;
    Real swapImprovement = 1.01;
    Real areaRelativeTolerance = 1e-12;
    std::size_t maxIterations = 5;
    std::size_t smoothingPasses = 0;
    std::size_t featureSmoothingPasses = 0;

    static TriangleMeshOptimizerParameters levelSetCarryForward(Real h)
    {
      TriangleMeshOptimizerParameters p;
      p.hMin = Real(0.50) * h;
      p.hMax = Real(2.00) * h;
      p.minQuality = Real(0.1);
      p.swapImprovement = Real(1.0001);
      p.maxIterations = 12;
      p.smoothingPasses = 3;
      p.featureSmoothingPasses = 3;
      return p;
    }
  };

  /**
   * @brief Minimal, bulletproof, parallel-ready 2D triangle optimizer.
   *
   * Three MMG-inspired local topological operators driven by edge-length
   * bounds: SPLIT (edge longer than hmax), COLLAPSE (edge shorter than hmin,
   * interior vertex removed), SWAP (quality-improving diagonal flip). It is not
   * an MMG port and has no MMG dependency.
   *
   * Bulletproof: an operation is committed only if every triangle it would
   * produce is strictly valid (orientation preserved, twice-area above a
   * scale-relative floor, quality at least @ref setMinQuality), and a collapse
   * additionally satisfies the manifold link condition. Invalid or topology-
   * breaking operations are rejected, never applied, so the optimizer can
   * never introduce an inverted, degenerate, duplicated or non-manifold cell.
   *
   * Parallel-ready architecture: each pass is split into a read-only EVALUATE
   * phase (a pure function of the current mesh over an edge array — trivially
   * `omp parallel for`) and a deterministic COMMIT phase that accepts a vertex-
   * disjoint independent set of proposals (cavity locking). Disjoint cavities
   * make the commit itself parallelizable; conflicting proposals are deferred
   * to the next outer iteration.
   *
   * The mesh is edited on a detached vertex/triangle copy and a fresh
   * LocalMesh is rebuilt (no in-place connectivity surgery, no MMG).
   */
  template <class Quality = ShapeQuality>
  class TriangleMeshOptimizer
  {
    public:
      TriangleMeshOptimizer() = default;

      explicit TriangleMeshOptimizer(Quality q) : m_quality(std::move(q)) {}

      TriangleMeshOptimizer& setHMin(Real h)
      { m_hmin = std::max(Real(0), h); return *this; }

      TriangleMeshOptimizer& setHMax(Real h)
      { m_hmax = std::max(Real(0), h); return *this; }

      TriangleMeshOptimizer& setMaxIterations(std::size_t n)
      { m_maxIterations = n; return *this; }

      TriangleMeshOptimizer& setMinQuality(Real q)
      { m_minQuality = std::max(Real(0), q); return *this; }

      TriangleMeshOptimizer& setSwapImprovement(Real r)
      { m_swapImprovement = std::max(Real(1), r); return *this; }

      /**
       * @brief Quality-improving vertex relocation sweeps per outer iteration.
       *
       * 0 (default) disables it. When > 0, after the topological operators
       * each free interior vertex (NOT protected, NOT on the domain boundary)
       * is moved toward its one-ring centroid only if the move strictly
       * improves the minimum incident-triangle quality and keeps every
       * incident triangle valid (orientation preserved, no inversion). This is
       * the optimization-based smoothing every production remesher has; it is
       * fit-neutral because the interface/boundary are frozen, and it cannot
       * invert because it is validity-gated.
       */
      TriangleMeshOptimizer& setSmoothingPasses(std::size_t n)
      { m_smoothingPasses = n; return *this; }

      using FeatureProjection =
        std::function<Math::SpatialPoint(const Math::SpatialPoint&)>;

      /**
       * @brief Enables tangential smoothing of protected (interface) vertices.
       *
       * Without this, protected vertices are frozen. With a projection set,
       * each protected non-domain-boundary vertex is smoothed toward its
       * one-ring centroid and then snapped back onto the feature by `proj`
       * (e.g. the analytic phi=0), so the fit stays exact while interface-
       * adjacent element quality improves. Validity-gated and monotone in
       * qmin. `n` = sweeps per outer iteration (0 default = off).
       */
      TriangleMeshOptimizer& setFeatureProjection(FeatureProjection proj)
      { m_featureProjection = std::move(proj); return *this; }

      TriangleMeshOptimizer& setFeatureSmoothingPasses(std::size_t n)
      { m_featureSmoothingPasses = n; return *this; }

      /// Enables collapse of SHORT feature (interface) edges along the
      /// feature, even though both endpoints are protected.
      ///
      /// Tangential smoothing redistributes protected vertices but never
      /// removes one, so a sliver wedged between two interface vertices
      /// survives and the downstream TMOP solve (fixed topology) cannot fix
      /// it. This pass collapses an interface edge shorter than @p length by
      /// merging its two endpoints into one vertex placed on the feature via
      /// setFeatureProjection (so the interface stays on phi=0 and one fewer
      /// node). Requires setFeatureProjection. A feature edge is detected
      /// geometrically: a short edge between two protected vertices whose
      /// midpoint already lies on the feature (projection moves it < 1/4 of
      /// the edge length). Bulletproof: manifold link condition + every
      /// rewritten triangle re-validated with the survivor at its new
      /// position. 0 (default) disables it.
      TriangleMeshOptimizer& setFeatureCollapseLength(Real length)
      { m_featureCollapseLength = std::max(Real(0), length); return *this; }

      /// Scale-relative twice-area floor used by the validity gate.
      TriangleMeshOptimizer& setAreaRelativeTolerance(Real r)
      { m_areaRel = std::max(Real(0), r); return *this; }

      TriangleMeshOptimizer& setQuality(Quality q)
      { m_quality = std::move(q); return *this; }

      TriangleMeshOptimizer& setParameters(
          const TriangleMeshOptimizerParameters& p)
      {
        return setHMin(p.hMin)
          .setHMax(p.hMax)
          .setMinQuality(p.minQuality)
          .setSwapImprovement(p.swapImprovement)
          .setAreaRelativeTolerance(p.areaRelativeTolerance)
          .setMaxIterations(p.maxIterations)
          .setSmoothingPasses(p.smoothingPasses)
          .setFeatureSmoothingPasses(p.featureSmoothingPasses);
      }

      /**
       * @brief Protects vertices from topological operators.
       *
       * No operator touches an edge with a protected endpoint: protected
       * vertices are never removed, relabeled, split across or swapped through.
       * Use this to preserve a fitted level-set interface (or any feature)
       * through coarsening/optimization. Protected vertices remain geometrically
       * fixed unless setFeatureProjection() and setFeatureSmoothingPasses() are
       * enabled, in which case feature smoothing may move them tangentially and
       * project them back onto the supplied feature.
       */
      TriangleMeshOptimizer& setProtectedVertices(std::vector<char> mask)
      { m_protected = std::move(mask); return *this; }

      TriangleMeshOptimizerReport optimize(LocalMesh& mesh) const
      {
        TriangleMeshOptimizerReport report;

        const auto& conn = mesh.getConnectivity();
        m_X.clear();
        m_X.reserve(mesh.getVertexCount());
        for (Index v = 0; v < static_cast<Index>(mesh.getVertexCount()); ++v)
          m_X.push_back(mesh.getVertexCoordinates(v));

        m_T.clear(); m_A.clear(); m_alive.clear();
        for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
        {
          if (conn.getGeometry(2, c) != Polytope::Type::Triangle)
            continue;
          const auto& cell = conn.getPolytope(2, c);
          m_T.push_back({{ cell(0), cell(1), cell(2) }});
          m_A.push_back(mesh.getCellAttribute(c));
          m_alive.push_back(1);
        }

        // Preserve attributed 1-D polytopes (e.g. the level-set interface):
        // their endpoints are expected to be frozen via setProtectedVertices,
        // so the edges survive the optimization and are re-stamped on rebuild.
        m_E.clear(); m_AE.clear();
        for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
        {
          const auto attr = mesh.getAttribute(1, e);
          if (!attr) continue;
          const auto& ed = conn.getPolytope(1, e);
          m_E.push_back({{ ed(0), ed(1) }});
          m_AE.push_back(attr);
        }

        report.minQualityBefore = minQuality();

        for (std::size_t it = 0; it < m_maxIterations; ++it)
        {
          ++report.iterations;
          const std::size_t s = m_hmax > Real(0) ? splitPass() : 0;
          const std::size_t w = swapPass();
          std::size_t sm = 0;
          for (std::size_t p = 0; p < m_smoothingPasses; ++p)
          {
            const std::size_t m = smoothPass();
            sm += m;
            if (m == 0) break;
          }
          std::size_t fsm = 0;
          for (std::size_t p = 0; p < m_featureSmoothingPasses; ++p)
          {
            const std::size_t m = featureSmoothPass();
            fsm += m;
            if (m == 0) break;
          }
          const std::size_t c = m_hmin > Real(0) ? collapsePass() : 0;
          const std::size_t fc =
            (m_featureCollapseLength > Real(0) && m_featureProjection)
              ? featureCollapsePass() : 0;
          report.splits += s; report.collapses += c + fc;
          report.swaps += w; report.smooths += sm;
          report.featureSmooths += fsm;
          if (s == 0 && c == 0 && fc == 0 && w == 0 && sm == 0 && fsm == 0)
            break;
        }

        report.minQualityAfter = minQuality();
        rebuild(mesh);
        return report;
      }

    private:
      using Tri = std::array<Index, 3>;

      // --- geometry / validity ---------------------------------------------

      Real q3(const Tri& t) const
      { return m_quality(m_X[t[0]], m_X[t[1]], m_X[t[2]]); }

      bool prot(Index v) const
      { return v < m_protected.size() && m_protected[v]; }

      // Strictly valid: orientation `sign` preserved, twice-area above a
      // scale-relative floor, quality at least the requested minimum.
      bool validTri(const Tri& t, Real sign) const
      {
        const auto& a = m_X[t[0]];
        const auto& b = m_X[t[1]];
        const auto& c = m_X[t[2]];
        const Real o = triangleOrient2(a, b, c);
        if ((o > 0) != (sign > 0))
          return false;
        const Real scale = std::max({ (b - a).squaredNorm(),
                                      (c - b).squaredNorm(),
                                      (a - c).squaredNorm() });
        if (std::abs(o) <= m_areaRel * scale)
          return false;
        return q3(t) >= m_minQuality;
      }

      bool strictlyOppositeSides(
          Index a,
          Index b,
          Index p,
          Index q) const
      {
        const auto& xa = m_X[a];
        const auto& xb = m_X[b];
        const auto& xp = m_X[p];
        const auto& xq = m_X[q];
        const Real op = triangleOrient2(xa, xb, xp);
        const Real oq = triangleOrient2(xa, xb, xq);
        const Real scale = std::max({
            (xb - xa).squaredNorm(),
            (xp - xa).squaredNorm(),
            (xq - xa).squaredNorm(),
            (xp - xb).squaredNorm(),
            (xq - xb).squaredNorm() });
        const Real eps = m_areaRel * scale;
        return (op > eps && oq < -eps) || (op < -eps && oq > eps);
      }

      bool validSwapCavity(Index u, Index v, Index a, Index b) const
      {
        if (u == v || a == b || u == a || u == b || v == a || v == b)
          return false;
        // The two old triangles are a valid local patch only when a and b
        // straddle the old diagonal uv, and the flipped patch is valid only
        // when u and v straddle the new diagonal ab. Per-triangle orientation
        // alone is weaker: a non-convex cavity can produce two positive
        // triangles that overlap.
        return strictlyOppositeSides(u, v, a, b)
            && strictlyOppositeSides(a, b, u, v);
      }

      Real minQuality() const
      {
        Real q = std::numeric_limits<Real>::infinity();
        bool any = false;
        for (std::size_t i = 0; i < m_T.size(); ++i)
          if (m_alive[i]) { q = std::min(q, q3(m_T[i])); any = true; }
        return any ? q : Real(0);
      }

      // --- adjacency (rebuilt per pass, cache-friendly) --------------------

      // Stable edge key: uses the vertex count captured at buildEdges() time,
      // so keys stay valid even after split appends midpoint vertices within
      // the same commit (edge proposals only ever reference pre-pass vertices).
      Index ekey(Index a, Index b) const
      {
        const Index n = m_keyN;
        return a < b ? a * n + b : b * n + a;
      }

      // Vertex -> incident live triangles, CSR layout.
      void buildVertexRing() const
      {
        const Index nv = m_X.size();
        m_vstart.assign(nv + 1, 0);
        for (std::size_t i = 0; i < m_T.size(); ++i)
          if (m_alive[i])
            for (Index v : m_T[i]) ++m_vstart[v + 1];
        for (Index v = 0; v < nv; ++v) m_vstart[v + 1] += m_vstart[v];
        m_vtri.assign(m_vstart[nv], 0);
        std::vector<Index> cur(m_vstart.begin(), m_vstart.end() - 1);
        for (std::size_t i = 0; i < m_T.size(); ++i)
          if (m_alive[i])
            for (Index v : m_T[i]) m_vtri[cur[v]++] = i;
      }

      template <class F>
      void forVertexTris(Index v, F&& f) const
      {
        for (Index k = m_vstart[v]; k < m_vstart[v + 1]; ++k)
          if (m_alive[m_vtri[k]]) f(m_vtri[k]);
      }

      // Edge -> up to two incident live triangles; >2 marks non-manifold.
      void buildEdges() const
      {
        m_keyN = m_X.size();
        m_edge.clear();
        m_edge.reserve(m_T.size() * 2);
        m_edges.clear();
        for (std::size_t i = 0; i < m_T.size(); ++i)
        {
          if (!m_alive[i]) continue;
          const auto& t = m_T[i];
          for (int e = 0; e < 3; ++e)
          {
            const Index a = t[e], b = t[(e + 1) % 3];
            auto [it, ins] = m_edge.try_emplace(ekey(a, b), EdgeRec{});
            auto& r = it->second;
            if (ins) { r.a = std::min(a, b); r.b = std::max(a, b); }
            if (r.n < 2) r.t[r.n] = i;
            ++r.n;
            if (ins) m_edges.push_back(it->first);
          }
        }
      }

      Index third(std::size_t ti, Index a, Index b) const
      {
        for (Index v : m_T[ti])
          if (v != a && v != b) return v;
        return std::numeric_limits<Index>::max();
      }

      struct EdgeRec { Index a = 0, b = 0; std::array<std::size_t, 2> t{{0,0}}; std::uint8_t n = 0; };
      // cav = the full set of vertices the operation rewrites. The commit
      // accepts a vertex-disjoint independent set, so cav MUST cover every
      // vertex whose incident triangles change — for a collapse that is the
      // entire 1-ring of the removed vertex, not just the 4 edge vertices.
      // Under-specifying cav lets two "independent" collapses with disjoint
      // edge endpoints but overlapping rings both commit and fold the mesh.
      struct Proposal { int kind = 0; Index a = 0, b = 0; Real gain = 0;
                        std::vector<Index> cav; };

      // --- SPLIT ----------------------------------------------------------

      std::size_t splitPass() const
      {
        buildEdges();
        std::vector<Proposal> prop(m_edges.size());
        std::vector<char> ok(m_edges.size(), 0);
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 256)
#endif
        for (std::ptrdiff_t e = 0;
             e < static_cast<std::ptrdiff_t>(m_edges.size()); ++e)
        {
          const auto& r = m_edge.at(m_edges[e]);
          if (prot(r.a) || prot(r.b)) continue;
          const Real len = (m_X[r.a] - m_X[r.b]).norm();
          if (len <= m_hmax || r.n == 0 || r.n > 2) continue;
          Proposal p; p.kind = 1; p.a = r.a; p.b = r.b; p.gain = len;
          p.cav = { r.a, r.b };
          for (std::uint8_t k = 0; k < r.n; ++k)
            p.cav.push_back(third(r.t[k], r.a, r.b));
          prop[e] = p; ok[e] = 1;
        }
        return commit(prop, ok, [&](const Proposal& p)
        {
          const Index w = static_cast<Index>(m_X.size());
          m_X.push_back(Real(0.5) * (m_X[p.a] + m_X[p.b]));
          const auto it = m_edge.find(ekey(p.a, p.b));
          if (it == m_edge.end())
          {
            m_X.pop_back();
            return false;
          }
          std::array<Tri, 4> nw; int cnt = 0;
          std::array<std::size_t, 2> kill; std::uint8_t nk = 0;
          for (std::uint8_t k = 0; k < it->second.n && k < 2; ++k)
          {
            const std::size_t ti = it->second.t[k];
            if (!m_alive[ti])
            {
              m_X.pop_back();
              return false;
            }
            Tri t1 = m_T[ti], t2 = m_T[ti];
            for (int j = 0; j < 3; ++j)
            { if (t1[j] == p.b) t1[j] = w; if (t2[j] == p.a) t2[j] = w; }
            const Real s = triangleOrient2(
                m_X[m_T[ti][0]], m_X[m_T[ti][1]], m_X[m_T[ti][2]]);
            if (!validTri(t1, s) || !validTri(t2, s))
            {
              m_X.pop_back();
              return false;
            }
            nw[cnt++] = t1; nw[cnt++] = t2; kill[nk++] = ti;
          }
          for (std::uint8_t k = 0; k < nk; ++k) m_alive[kill[k]] = 0;
          for (int i = 0; i < cnt; ++i)
          {
            m_T.push_back(nw[i]);
            m_A.push_back(m_A[kill[i / 2]]);
            m_alive.push_back(1);
          }
          return true;
        });
      }

      // --- COLLAPSE -------------------------------------------------------

      std::size_t collapsePass() const
      {
        buildVertexRing();
        buildEdges();
        const Index nv = m_X.size();
        std::vector<char> bnd(nv, 0);
        for (const Index ek : m_edges)
        { const auto& r = m_edge.at(ek);
          if (r.n == 1) { bnd[r.a] = 1; bnd[r.b] = 1; } }

        std::vector<Proposal> prop(m_edges.size());
        std::vector<char> ok(m_edges.size(), 0);
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 256)
#endif
        for (std::ptrdiff_t e = 0;
             e < static_cast<std::ptrdiff_t>(m_edges.size()); ++e)
        {
          const auto& r = m_edge.at(m_edges[e]);
          if (prot(r.a) || prot(r.b)) continue;
          const Real len = (m_X[r.a] - m_X[r.b]).norm();
          if (len >= m_hmin || r.n != 2) continue;     // interior edges only
          Index rem = r.a, keep = r.b;
          if (bnd[rem] && !bnd[keep]) std::swap(rem, keep);
          if (bnd[rem]) continue;                       // never move boundary
          if (!linkCondition(rem, keep, r)) continue;
          if (!collapseValid(rem, keep)) continue;
          Proposal p; p.kind = 2; p.a = rem; p.b = keep;
          p.gain = m_hmin - len;
          // Cavity = the ENTIRE 1-ring of rem (every vertex of every triangle
          // that the collapse rewrites) plus keep. This is what makes the
          // independent-set commit actually correct: two collapses can only
          // both fire if their full rings are vertex-disjoint, so they cannot
          // jointly fold the mesh.
          p.cav = { rem, keep };
          forVertexTris(rem, [&](std::size_t ti)
          {
            for (Index x : m_T[ti]) p.cav.push_back(x);
          });
          prop[e] = p; ok[e] = 1;
        }
        return commit(prop, ok, [&](const Proposal& p)
        {
          const Index rem = p.a, keep = p.b;
          // Re-validate against current state (ring may have changed).
          if (!collapseValid(rem, keep)) return false;
          std::vector<std::size_t> ring;
          forVertexTris(rem, [&](std::size_t ti){ ring.push_back(ti); });
          for (std::size_t ti : ring)
          {
            auto& t = m_T[ti];
            const bool hasK = (t[0]==keep||t[1]==keep||t[2]==keep);
            if (hasK) { m_alive[ti] = 0; continue; }
            for (int j = 0; j < 3; ++j) if (t[j]==rem) t[j]=keep;
          }
          return true;
        });
      }

      // Manifold link condition for edge collapse rem->keep.
      bool linkCondition(Index rem, Index keep, const EdgeRec& r) const
      {
        const Index p0 = third(r.t[0], rem, keep);
        const Index p1 = third(r.t[1], rem, keep);
        // Neighbors of keep.
        bool ok = true;
        forVertexTris(rem, [&](std::size_t ti)
        {
          for (Index v : m_T[ti])
          {
            if (v == rem || v == keep || v == p0 || v == p1) continue;
            // v in link(rem). If v is also adjacent to keep -> fold.
            bool adjKeep = false;
            forVertexTris(keep, [&](std::size_t tj)
            {
              for (Index w : m_T[tj]) if (w == v) adjKeep = true;
            });
            if (adjKeep) ok = false;
          }
        });
        return ok;
      }

      // Every triangle around `rem` (not also on `keep`) must stay valid
      // after rem is moved onto keep.
      bool collapseValid(Index rem, Index keep) const
      {
        bool good = true;
        forVertexTris(rem, [&](std::size_t ti)
        {
          if (!good) return;
          const auto& t = m_T[ti];
          if (t[0]==keep||t[1]==keep||t[2]==keep) return; // dies away
          Tri nt = t;
          const Real s = triangleOrient2(m_X[t[0]], m_X[t[1]], m_X[t[2]]);
          for (int j = 0; j < 3; ++j) if (nt[j]==rem) nt[j]=keep;
          if (!validTri(nt, s)) good = false;
        });
        return good;
      }

      // --- FEATURE COLLAPSE (collapse a short interface edge ALONG phi=0) -

      // validTri but with explicit vertex positions (the survivor moves).
      bool triValidPos(
          const std::array<Math::SpatialPoint, 3>& p, Real sign) const
      {
        const Real o = triangleOrient2(p[0], p[1], p[2]);
        if ((o > 0) != (sign > 0)) return false;
        const Real scale = std::max({ (p[1] - p[0]).squaredNorm(),
                                      (p[2] - p[1]).squaredNorm(),
                                      (p[0] - p[2]).squaredNorm() });
        if (std::abs(o) <= m_areaRel * scale) return false;
        return m_quality(p[0], p[1], p[2]) >= m_minQuality;
      }

      // Collapse rem->keep AND relocate the survivor to `nk`. Every triangle
      // around either endpoint (except the two that die on the shared edge)
      // must stay strictly valid with the new vertex positions.
      bool featureCollapseValid(
          Index rem, Index keep, const Math::SpatialPoint& nk) const
      {
        bool good = true;
        auto checkRing = [&](Index pivot)
        {
          forVertexTris(pivot, [&](std::size_t ti)
          {
            if (!good) return;
            const auto& t = m_T[ti];
            const bool hasRem = (t[0]==rem||t[1]==rem||t[2]==rem);
            const bool hasKeep = (t[0]==keep||t[1]==keep||t[2]==keep);
            if (hasRem && hasKeep) return;       // dies on the collapsed edge
            const Real s =
              triangleOrient2(m_X[t[0]], m_X[t[1]], m_X[t[2]]);
            std::array<Math::SpatialPoint, 3> p;
            for (int j = 0; j < 3; ++j)
            {
              Index v = t[j];
              if (v == rem) v = keep;            // rem maps onto keep
              p[j] = (v == keep) ? nk : m_X[v];  // survivor at new position
            }
            if (!triValidPos(p, s)) good = false;
          });
        };
        checkRing(rem);
        checkRing(keep);
        return good;
      }

      std::size_t featureCollapsePass() const
      {
        buildVertexRing();
        buildEdges();
        const Index nv = m_X.size();
        std::vector<char> bnd(nv, 0);
        for (const Index ek : m_edges)
        { const auto& r = m_edge.at(ek);
          if (r.n == 1) { bnd[r.a] = 1; bnd[r.b] = 1; } }

        std::vector<Proposal> prop(m_edges.size());
        std::vector<char> ok(m_edges.size(), 0);
        for (std::size_t e = 0; e < m_edges.size(); ++e)
        {
          const auto& r = m_edge.at(m_edges[e]);
          if (!prot(r.a) || !prot(r.b)) continue;   // BOTH on the feature
          if (r.n != 2) continue;                   // interior only
          if (bnd[r.a] || bnd[r.b]) continue;       // not domain boundary
          const Real len = (m_X[r.a] - m_X[r.b]).norm();
          if (len >= m_featureCollapseLength || len <= Real(0)) continue;
          // Feature-edge test: the midpoint must already lie on the feature.
          const Math::SpatialPoint mid =
            Real(0.5) * (m_X[r.a] + m_X[r.b]);
          const Math::SpatialPoint proj = m_featureProjection(mid);
          if ((proj - mid).norm() > Real(0.25) * len) continue;
          const Index rem = r.a, keep = r.b;
          if (!linkCondition(rem, keep, r)) continue;
          if (!featureCollapseValid(rem, keep, proj)) continue;
          Proposal p; p.kind = 2; p.a = rem; p.b = keep;
          p.gain = m_featureCollapseLength - len;  // shortest first
          p.cav = { rem, keep };
          forVertexTris(rem, [&](std::size_t ti)
          { for (Index x : m_T[ti]) p.cav.push_back(x); });
          forVertexTris(keep, [&](std::size_t ti)
          { for (Index x : m_T[ti]) p.cav.push_back(x); });
          prop[e] = p; ok[e] = 1;
        }
        return commit(prop, ok, [&](const Proposal& p)
        {
          const Index rem = p.a, keep = p.b;
          const Math::SpatialPoint mid =
            Real(0.5) * (m_X[rem] + m_X[keep]);
          const Math::SpatialPoint proj = m_featureProjection(mid);
          if (!featureCollapseValid(rem, keep, proj)) return false;
          std::vector<std::size_t> ring;
          forVertexTris(rem, [&](std::size_t ti){ ring.push_back(ti); });
          for (std::size_t ti : ring)
          {
            auto& t = m_T[ti];
            const bool hasK = (t[0]==keep||t[1]==keep||t[2]==keep);
            if (hasK) { m_alive[ti] = 0; continue; }
            for (int j = 0; j < 3; ++j) if (t[j]==rem) t[j]=keep;
          }
          m_X[keep] = proj;                        // survivor lands on phi=0
          // Keep the preserved interface-edge list consistent so rebuild()
          // re-stamps the Interface attribute onto the surviving edges
          // (rem -> keep). Entries that become {keep,keep} are the collapsed
          // edge itself; they never match a real edge at re-stamp and are
          // harmlessly ignored.
          for (auto& pr : m_E)
          {
            if (pr[0] == rem) pr[0] = keep;
            if (pr[1] == rem) pr[1] = keep;
          }
          return true;
        });
      }

      // --- SWAP -----------------------------------------------------------

      std::size_t swapPass() const
      {
        buildEdges();
        std::vector<Proposal> prop(m_edges.size());
        std::vector<char> ok(m_edges.size(), 0);
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 256)
#endif
        for (std::ptrdiff_t e = 0;
             e < static_cast<std::ptrdiff_t>(m_edges.size()); ++e)
        {
          const auto& r = m_edge.at(m_edges[e]);
          if (r.n != 2) continue;                       // interior only
          if (m_A[r.t[0]] != m_A[r.t[1]]) continue;      // do not mix regions
          const Index u = r.a, v = r.b;
          if (prot(u) || prot(v)) continue;             // keep feature edges
          const Index a = third(r.t[0], u, v);
          const Index b = third(r.t[1], u, v);
          if (a == b) continue;
          if (!validSwapCavity(u, v, a, b)) continue;
          const Real qOld = std::min(q3(m_T[r.t[0]]), q3(m_T[r.t[1]]));
          Tri n0 = {{ a, u, b }}, n1 = {{ b, v, a }};
          const Real s0 = triangleOrient2(m_X[a], m_X[u], m_X[b]);
          const Real s1 = triangleOrient2(m_X[b], m_X[v], m_X[a]);
          if (s0 == 0 || s1 == 0) continue;
          if (!validTri(n0, s0) || !validTri(n1, s1)) continue;
          const Real qNew = std::min(q3(n0), q3(n1));
          if (qNew <= m_swapImprovement * qOld) continue;
          Proposal p; p.kind = 3; p.a = u; p.b = v; p.gain = qNew - qOld;
          p.cav = { u, v, a, b };
          prop[e] = p; ok[e] = 1;
        }
        return commit(prop, ok, [&](const Proposal& p)
        {
          const auto it = m_edge.find(ekey(p.a, p.b));
          const auto& r = it->second;
          if (r.n != 2 || !m_alive[r.t[0]] || !m_alive[r.t[1]]) return false;
          if (m_A[r.t[0]] != m_A[r.t[1]]) return false;
          const Index u = p.a, v = p.b;
          const Index a = third(r.t[0], u, v), b = third(r.t[1], u, v);
          if (!validSwapCavity(u, v, a, b)) return false;
          Tri n0 = {{ a, u, b }}, n1 = {{ b, v, a }};
          const Real s0 = triangleOrient2(m_X[a], m_X[u], m_X[b]);
          const Real s1 = triangleOrient2(m_X[b], m_X[v], m_X[a]);
          if (!validTri(n0, s0) || !validTri(n1, s1)) return false;
          m_T[r.t[0]] = n0; m_T[r.t[1]] = n1;
          return true;
        });
      }

      // --- COMMIT (deterministic independent set; parallel-ready) ---------

      template <class Apply>
      std::size_t commit(std::vector<Proposal>& prop,
                         const std::vector<char>& ok, Apply&& apply) const
      {
        std::vector<std::size_t> order;
        order.reserve(prop.size());
        for (std::size_t i = 0; i < prop.size(); ++i)
          if (ok[i]) order.push_back(i);
        std::sort(order.begin(), order.end(),
          [&](std::size_t x, std::size_t y)
          { return prop[x].gain > prop[y].gain; });

        std::vector<char> locked(m_X.size(), 0);
        std::size_t n = 0;
        for (std::size_t idx : order)
        {
          const Proposal& p = prop[idx];
          bool free = true;
          for (Index cv : p.cav)
            if (cv < locked.size() && locked[cv]) { free = false; break; }
          if (!free) continue;                  // disjoint-cavity scheduling
          if (!apply(p)) continue;              // bulletproof re-validation
          for (Index cv : p.cav)
            if (cv < locked.size()) locked[cv] = 1;
          ++n;
        }
        return n;
      }

      // --- SMOOTH (quality-improving vertex relocation) -------------------

      // Sequential optimization-based smoothing: each free interior vertex is
      // pulled toward its one-ring centroid, accepted only if it strictly
      // improves the minimum incident-triangle quality and every incident
      // triangle stays valid (orientation preserved, no inversion). Protected
      // (interface) and domain-boundary vertices are frozen, so the fit is
      // untouched. Monotone in qmin by construction.
      std::size_t smoothPass() const
      {
        buildVertexRing();
        buildEdges();
        const Index nv = m_X.size();
        std::vector<char> bnd(nv, 0);
        for (const Index ek : m_edges)
        {
          const auto& r = m_edge.at(ek);
          if (r.n == 1) { bnd[r.a] = 1; bnd[r.b] = 1; }
        }

        std::size_t moved = 0;
        std::vector<std::size_t> ring;
        std::vector<Index> nbr;
        std::vector<Real> oldSign;
        for (Index v = 0; v < nv; ++v)
        {
          if (prot(v) || bnd[v]) continue;
          ring.clear(); nbr.clear(); oldSign.clear();
          forVertexTris(v, [&](std::size_t ti) { ring.push_back(ti); });
          if (ring.empty()) continue;

          Real qOld = std::numeric_limits<Real>::infinity();
          for (std::size_t ti : ring)
          {
            const auto& t = m_T[ti];
            qOld = std::min(qOld, q3(t));
            oldSign.push_back(
                triangleOrient2(m_X[t[0]], m_X[t[1]], m_X[t[2]]));
            for (Index w : t)
              if (w != v
                  && std::find(nbr.begin(), nbr.end(), w) == nbr.end())
                nbr.push_back(w);
          }
          if (nbr.empty()) continue;

          Math::SpatialPoint centroid = m_X[nbr[0]];
          for (std::size_t k = 1; k < nbr.size(); ++k)
            centroid = centroid + m_X[nbr[k]];
          centroid = (Real(1) / static_cast<Real>(nbr.size())) * centroid;

          const Math::SpatialPoint orig = m_X[v];
          bool accepted = false;
          for (const Real f : { Real(1), Real(0.5), Real(0.25) })
          {
            m_X[v] = orig + f * (centroid - orig);
            Real qNew = std::numeric_limits<Real>::infinity();
            bool ok = true;
            for (std::size_t r = 0; r < ring.size() && ok; ++r)
            {
              const auto& t = m_T[ring[r]];
              if (!validTri(t, oldSign[r])) ok = false;
              else qNew = std::min(qNew, q3(t));
            }
            if (ok && qNew > qOld * (Real(1) + Real(1e-12)))
            { accepted = true; ++moved; break; }
          }
          if (!accepted) m_X[v] = orig;
        }
        return moved;
      }

      // Tangential feature smoothing: relocate PROTECTED (interface) vertices
      // toward their one-ring centroid but snap the candidate back onto the
      // feature via m_featureProjection (e.g. the analytic phi=0), so the fit
      // stays exact while interface-adjacent element quality improves.
      // Domain-boundary vertices stay frozen. Inert if no projection is set.
      std::size_t featureSmoothPass() const
      {
        if (!m_featureProjection) return 0;
        buildVertexRing();
        buildEdges();
        const Index nv = m_X.size();
        std::vector<char> bnd(nv, 0);
        for (const Index ek : m_edges)
        {
          const auto& r = m_edge.at(ek);
          if (r.n == 1) { bnd[r.a] = 1; bnd[r.b] = 1; }
        }

        std::size_t moved = 0;
        std::vector<std::size_t> ring;
        std::vector<Index> nbr;
        std::vector<Real> oldSign;
        for (Index v = 0; v < nv; ++v)
        {
          if (!prot(v) || bnd[v]) continue;       // interface vertices only
          ring.clear(); nbr.clear(); oldSign.clear();
          forVertexTris(v, [&](std::size_t ti) { ring.push_back(ti); });
          if (ring.empty()) continue;

          Real qOld = std::numeric_limits<Real>::infinity();
          for (std::size_t ti : ring)
          {
            const auto& t = m_T[ti];
            qOld = std::min(qOld, q3(t));
            oldSign.push_back(
                triangleOrient2(m_X[t[0]], m_X[t[1]], m_X[t[2]]));
            for (Index w : t)
              if (w != v
                  && std::find(nbr.begin(), nbr.end(), w) == nbr.end())
                nbr.push_back(w);
          }
          if (nbr.empty()) continue;

          Math::SpatialPoint centroid = m_X[nbr[0]];
          for (std::size_t k = 1; k < nbr.size(); ++k)
            centroid = centroid + m_X[nbr[k]];
          centroid = (Real(1) / static_cast<Real>(nbr.size())) * centroid;

          const Math::SpatialPoint orig = m_X[v];
          bool accepted = false;
          for (const Real f : { Real(1), Real(0.5), Real(0.25) })
          {
            const Math::SpatialPoint cand = orig + f * (centroid - orig);
            m_X[v] = m_featureProjection(cand);   // snap back onto phi=0
            Real qNew = std::numeric_limits<Real>::infinity();
            bool ok = true;
            for (std::size_t r = 0; r < ring.size() && ok; ++r)
            {
              const auto& t = m_T[ring[r]];
              if (!validTri(t, oldSign[r])) ok = false;
              else qNew = std::min(qNew, q3(t));
            }
            if (ok && qNew > qOld * (Real(1) + Real(1e-12)))
            { accepted = true; ++moved; break; }
          }
          if (!accepted) m_X[v] = orig;
        }
        return moved;
      }

      // --- rebuild --------------------------------------------------------

      void rebuild(LocalMesh& mesh) const
      {
        std::vector<Index> remap(m_X.size(),
                                 std::numeric_limits<Index>::max());
        std::vector<Math::SpatialPoint> verts;
        std::size_t cells = 0;
        for (std::size_t i = 0; i < m_T.size(); ++i)
        {
          if (!m_alive[i]) continue;
          ++cells;
          for (Index v : m_T[i])
            if (remap[v] == std::numeric_limits<Index>::max())
            { remap[v] = verts.size(); verts.push_back(m_X[v]); }
        }

        LocalMesh::Builder builder;
        builder.initialize(2).nodes(verts.size()).reserve(2, cells);
        for (const auto& pt : verts) builder.vertex(pt);
        Index ci = 0;
        for (std::size_t i = 0; i < m_T.size(); ++i)
        {
          if (!m_alive[i]) continue;
          Index x0 = remap[m_T[i][0]], x1 = remap[m_T[i][1]],
                x2 = remap[m_T[i][2]];
          if (triangleOrient2(verts[x0], verts[x1], verts[x2]) < Real(0))
            std::swap(x1, x2);
          builder.polytope(Polytope::Type::Triangle, { x0, x1, x2 })
                 .attribute({ 2, ci }, m_A[i]);
          ++ci;
        }
        mesh = builder.finalize();
        mesh.getConnectivity().compute(2, 1);
        mesh.getConnectivity().compute(1, 0);

        // Re-stamp preserved attributed edges (the interface). Their endpoints
        // are frozen, so the edge still exists in the rebuilt triangulation;
        // we only need to locate it and restore its attribute.
        if (!m_E.empty())
        {
          const auto& oc = mesh.getConnectivity();
          UnorderedMap<Index, Index> pair2edge;
          pair2edge.reserve(oc.getCount(1) * 2);
          const Index n = mesh.getVertexCount();
          for (Index e = 0; e < static_cast<Index>(oc.getCount(1)); ++e)
          {
            const auto& ed = oc.getPolytope(1, e);
            const Index a = ed(0), b = ed(1);
            pair2edge[a < b ? a * n + b : b * n + a] = e;
          }
          for (std::size_t k = 0; k < m_E.size(); ++k)
          {
            const Index v0 = m_E[k][0], v1 = m_E[k][1];
            if (remap[v0] == std::numeric_limits<Index>::max()
                || remap[v1] == std::numeric_limits<Index>::max())
              continue;                       // region legitimately coarsened
            const Index a = remap[v0], b = remap[v1];
            const auto it =
              pair2edge.find(a < b ? a * n + b : b * n + a);
            if (it != pair2edge.end())
              mesh.setAttribute({ 1, it->second }, m_AE[k]);
          }
        }
      }

      Quality m_quality{};
      Real m_hmin = 0;
      Real m_hmax = std::numeric_limits<Real>::infinity();
      Real m_minQuality = 1e-4;
      Real m_swapImprovement = 1.01;
      Real m_areaRel = 1e-12;
      std::size_t m_maxIterations = 5;
      std::size_t m_smoothingPasses = 0;
      std::size_t m_featureSmoothingPasses = 0;
      Real m_featureCollapseLength = 0;
      FeatureProjection m_featureProjection;

      mutable std::vector<Math::SpatialPoint> m_X;
      mutable std::vector<Tri> m_T;
      mutable std::vector<Optional<Attribute>> m_A;
      mutable std::vector<std::array<Index, 2>> m_E;
      mutable std::vector<Optional<Attribute>> m_AE;
      mutable std::vector<char> m_alive;
      mutable UnorderedMap<Index, EdgeRec> m_edge;
      mutable std::vector<Index> m_edges;
      mutable std::vector<Index> m_vstart;
      mutable std::vector<Index> m_vtri;
      mutable Index m_keyN = 1;
      std::vector<char> m_protected;
  };

  TriangleMeshOptimizer() -> TriangleMeshOptimizer<ShapeQuality>;
}

#endif
