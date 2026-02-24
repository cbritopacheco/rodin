/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_MARCHINGTETRAHEDRA_H
#define RODIN_GEOMETRY_MARCHINGTETRAHEDRA_H

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <tuple>
#include <type_traits>
#include <variant>
#include <vector>

#include "Rodin/Alert/MemberFunctionException.h"

#include "MarchingBase.h"
#include "Mesh.h"
#include "Rodin/Geometry/Polytope.h"

namespace Rodin::Geometry
{
  /**
   * @brief Marching tetrahedra discretizer for a level-set on a 3D tetrahedral mesh.
   *
   * @tparam Params Specialization parameters.
   */
  template <class ... Params>
  class MarchingTetrahedra;

  /**
   * @brief Specialization for P1 grid functions.
   *
   * This class implements a robust, connectivity-preserving variant of “marching tetrahedra”
   * for a scalar P1 grid function @f$\phi@f$ on a tetrahedral volume mesh.
   *
   * The algorithm:
   * - iterates over input tetrahedra,
   * - classifies vertices by the sign of @f$\phi@f$ (with tolerance),
   * - if a strict sign change is detected, splits the tetrahedron into a small set of tetrahedra
   *   whose union matches the original tetrahedron, and whose new vertices lie on edges where
   *   @f$\phi@f=0@f$ (linear interpolation),
   * - assigns output cell attributes using the SplitMap semantics inherited from MarchingBase,
   * - optionally marks interface faces/edges with a dedicated interface attribute.
   *
   * Robustness measures included (relative to a naive marching implementation):
   * - sign tolerance and snap tolerance,
   * - edge intersection caching by input edge id (topologically consistent),
   * - deterministic tie-breaking for near-zero values (to avoid random “zero” patterns),
   * - best-of-two decomposition for the 2-neg/2-pos case based on a volume quality proxy,
   * - conservative attribute transfer: only original faces/edges inherit non-interface attributes.
   *
   * @tparam Mesh Mesh type (must be 3D tetrahedral volume mesh).
   * @tparam Data Underlying storage/backend of the GridFunction.
   */
  template <class Mesh, class Data>
  class MarchingTetrahedra<Variational::GridFunction<Variational::P1<Real, Mesh>, Data>>
    : public MarchingBase<Variational::GridFunction<Variational::P1<Real, Mesh>, Data>>
  {
    public:
      using Parent = MarchingBase<Variational::GridFunction<Variational::P1<Real, Mesh>, Data>>;

      using FESType          = typename Parent::FESType;
      using MeshType         = typename Parent::MeshType;
      using GridFunctionType = typename Parent::GridFunctionType;

      using SplitMap = typename Parent::SplitMap;
      using NoSplitT = typename Parent::NoSplitT;
      using Split    = typename Parent::Split;

      /**
       * @brief Construct from a P1 level-set grid function.
       *
       * @param ls Input P1 scalar field (level set).
       */
      MarchingTetrahedra(const GridFunctionType& ls)
        : Parent(ls),
          m_sign_tolerance(1e-12),
          m_snap_tolerance(1e-12)
      {}

      /**
       * @brief Discretize the implicit geometry described by the level set.
       *
       * Produces a tetrahedral mesh where input tetrahedra may be split along the
       * @f$\phi=0@f$ interface. Output cell attributes are inferred from the input cell
       * attribute using the SplitMap configuration (see MarchingBase).
       *
       * Interface marking:
       * - If an interface attribute is configured (setInterface), output faces/edges that
       *   separate negative and positive cells (or are “fitted” to zero) may be marked with it.
       *
       * Attribute transfer policy for faces/edges:
       * - Only *original* faces/edges (composed solely of original vertices) may inherit
       *   non-interface attributes from the input mesh.
       * - Newly created faces/edges remain un-attributed unless they are detected as interface.
       *
       * @return Output mesh resulting from the marching tetrahedra discretization.
       */
      MeshType discretize() const override
      {
        const auto& ls   = this->getGridFunction();
        const auto& mesh = ls.getFiniteElementSpace().getMesh();

        // We require adjacency/incidence needed for:
        // - reading boundary constraints (noSplit on edges/faces),
        // - identifying original boundary entities by vertices,
        // - rebuilding interface tagging by adjacency in the output.
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 3, 2);
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 3, 1);
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 2, 0);
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 1, 0);

        if (mesh.getDimension() != 3)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "MarchingTetrahedra expects a 3D volume mesh."
            << Alert::Raise;
        }

        auto& conn = mesh.getConnectivity();

        const Index nv = static_cast<Index>(mesh.getVertexCount());
        const Index ne = static_cast<Index>(mesh.getPolytopeCount(1));

        using Point = Math::SpatialPoint;
        static constexpr Index InvalidIndex = std::numeric_limits<Index>::max();

        // User-configured tolerances:
        // - eps_sign_user: classification threshold for deciding sign/zero.
        // - eps_snap_user: snapping threshold for edge intersections (avoid tiny slivers).
        const Real eps_sign_user = this->getSignTolerance();
        const Real eps_snap_user = this->getSnapTolerance();

        auto finite = [](Real x) -> bool { return std::isfinite(x); };

        // “cut” predicate consistent with legacy behavior:
        // - if both ends are classified as zero, no cut
        // - otherwise, any sign difference (including zero vs nonzero) is a cut.
        auto isCut = [&](int sa, int sb) -> bool
        {
          if (sa == 0 && sb == 0) return false;
          return sa != sb;
        };

        // Per-output-cell side (used for interface tagging in the output):
        // - Unknown: not confidently classified (or blocked from splitting)
        // - Positive: phi > 0 side
        // - Negative: phi < 0 side
        static constexpr int8_t SideUnknown  = -1;
        static constexpr int8_t SidePositive = 0;
        static constexpr int8_t SideNegative = 1;

        // Determine cell-side label from vertex signs (ignoring zero only through classification).
        auto sideFromSigns = [&](const int ss[4]) -> int8_t
        {
          bool hasNeg = false, hasPos = false;
          for (int i = 0; i < 4; ++i)
          {
            hasNeg = hasNeg || (ss[i] < 0);
            hasPos = hasPos || (ss[i] > 0);
          }
          if (hasNeg && hasPos) return SideUnknown;
          if (hasNeg)           return SideNegative;
          if (hasPos)           return SidePositive;
          return SideUnknown;
        };

        // --------------------------
        // SplitMap helpers
        // --------------------------

        // Return true if attribute a is explicitly configured as NoSplit at dimension d.
        auto isNoSplit = [&](size_t d, const Optional<Attribute>& a) -> bool
        {
          if (!a) return false;
          const auto& sm = this->getSplitMap(d);
          if (sm.empty()) return false;
          auto it = sm.find(*a);
          if (it == sm.end()) return false;
          return std::holds_alternative<NoSplitT>(it->second);
        };

        // Return true if we are allowed to split elements with attribute a at dimension d.
        // Policy:
        // - empty map => split allowed
        // - missing attribute => split allowed (default)
        // - present:
        //   * Split => allowed
        //   * NoSplit => blocked
        auto shouldSplit = [&](size_t d, const Optional<Attribute>& a) -> bool
        {
          const auto& sm = this->getSplitMap(d);
          if (sm.empty()) return true;
          if (!a) return true; // no provenance => still split geometry
          auto it = sm.find(*a);
          if (it == sm.end()) return true; // default: split
          return std::holds_alternative<Split>(it->second);
        };

        // Infer output label from an input attribute and a side.
        // If base is not in SplitMap, the attribute is kept unchanged.
        auto inferLabel = [&](size_t d, const Optional<Attribute>& base, bool negativeSide)
          -> Optional<Attribute>
        {
          if (!base) return {}; // no provenance => cannot label
          const auto& sm = this->getSplitMap(d);
          if (!sm.empty())
          {
            auto it = sm.find(*base);
            if (it != sm.end())
            {
              return std::visit(
                [&](const auto& v) -> Optional<Attribute>
                {
                  using T = std::decay_t<decltype(v)>;
                  if constexpr (std::is_same_v<T, NoSplitT>)
                    return base;
                  else
                    return negativeSide ? Optional<Attribute>(v.negative)
                                        : Optional<Attribute>(v.positive);
                },
                it->second);
            }
          }
          return base;
        };

        // Signed volume (6x actual tetra volume). Used to enforce orientation and quality choice.
        auto signedVolume = [&](Index a, Index b, Index c, Index d) -> Real
        {
          const auto A = mesh.getVertexCoordinates(a);
          const auto B = mesh.getVertexCoordinates(b);
          const auto C = mesh.getVertexCoordinates(c);
          const auto D = mesh.getVertexCoordinates(d);
          const auto u = B - A;
          const auto v = C - A;
          const auto w = D - A;
          return u.dot(v.cross(w));
        };

        // ------------------------------------------------------------
        // Attribute provenance maps for original edges/faces
        // ------------------------------------------------------------
        struct EdgeKey
        {
          Index a, b;
          bool operator<(const EdgeKey& o) const { return std::tie(a, b) < std::tie(o.a, o.b); }
        };

        struct FaceKey
        {
          Index a, b, c; // sorted triple
          bool operator<(const FaceKey& o) const { return std::tie(a, b, c) < std::tie(o.a, o.b, o.c); }
        };

        auto makeEdgeKey = [](Index i, Index j) -> EdgeKey
        {
          if (i > j) std::swap(i, j);
          return { i, j };
        };

        auto makeFaceKey = [](Index i, Index j, Index k) -> FaceKey
        {
          std::array<Index, 3> s{ i, j, k };
          std::sort(s.begin(), s.end());
          return { s[0], s[1], s[2] };
        };

        FlatMap<EdgeKey, Attribute> inEdgeAttrByVerts;
        FlatMap<FaceKey, Attribute> inFaceAttrByVerts;

        // Record attributes only for ORIGINAL edges/faces (in the input mesh).
        for (auto eit = mesh.getPolytope(1); !eit.end(); ++eit)
        {
          const auto& ev = eit->getVertices();
          const auto a = eit->getAttribute();
          if (a) inEdgeAttrByVerts[makeEdgeKey(ev(0), ev(1))] = *a;
        }

        for (auto fit = mesh.getPolytope(2); !fit.end(); ++fit)
        {
          if (fit->getGeometry() != Polytope::Type::Triangle) continue;
          const auto& fv = fit->getVertices();
          const auto a = fit->getAttribute();
          if (a) inFaceAttrByVerts[makeFaceKey(fv(0), fv(1), fv(2))] = *a;
        }

        // ------------------------------------------------------------
        // Output vertex set + edge intersection cache
        // ------------------------------------------------------------
        PointCloud outVerts;
        outVerts.setDimension(3);
        outVerts.reserve(static_cast<size_t>(nv) + static_cast<size_t>(ne / 4));

        // First, copy original vertices 1:1.
        for (Index i = 0; i < nv; ++i)
          outVerts.push_back(mesh.getVertexCoordinates(i));

        auto addVertex = [&](const Point& p) -> Index
        {
          outVerts.push_back(p);
          return static_cast<Index>(outVerts.getCount() - 1);
        };

        // Cache intersection vertex per input edge id. This enforces topological consistency:
        // the same input edge is intersected once and reused across all incident cells.
        std::vector<Index> isect(static_cast<size_t>(ne), InvalidIndex);

        /**
         * @brief Compute (or reuse) the intersection vertex on an input edge.
         *
         * The intersection is computed by linear interpolation of the P1 level-set:
         * @f$\phi((1-t)\,a + t\,b) = 0@f$.
         *
         * Robustness:
         * - snap to endpoints if |phi| <= epsS,
         * - clamp and snap the interpolation parameter,
         * - if phi(a) and phi(b) are nearly symmetric around zero, snap to mid-edge.
         *
         * @param edgeIdx Input edge index (for caching).
         * @param va,vb Input vertex indices (original vertices).
         * @param fa,fb Level-set values at va and vb.
         * @param eps0 Per-cell classification tolerance.
         * @param eps_snap Per-cell snapping tolerance (>= eps0).
         * @return Index of the intersection vertex in the output vertex array.
         */
        auto getIntersectionOnEdge = [&](Index edgeIdx,
                                         Index va, Index vb,
                                         Real fa, Real fb,
                                         Real eps0, Real eps_snap) -> Index
        {
          if (!finite(fa) || !finite(fb))
            return va;

          const Real epsS = std::max(eps_snap, eps0);

          if (std::abs(fa) <= epsS) return va;
          if (std::abs(fb) <= epsS) return vb;

          Index& cached = isect[static_cast<size_t>(edgeIdx)];
          if (cached != InvalidIndex)
            return cached;

          const Real denom = (fa - fb);
          if (std::abs(denom) <= epsS)
            return (std::abs(fa) < std::abs(fb)) ? va : vb;

          // (1-t)*fa + t*fb = 0 => t = fa/(fa-fb)
          Real t = fa / denom;
          t = std::min<Real>(Real(1), std::max<Real>(Real(0), t));

          // snap parameter close to endpoints
          if (t <= epsS)           return va;
          if (t >= Real(1) - epsS) return vb;

          // symmetric near-zero case => stabilize to mid-edge
          if (std::abs(fa + fb) <= Real(4) * eps0)
            t = Real(0.5);

          const auto pa = outVerts[va];
          const auto pb = outVerts[vb];
          const auto p  = pa + t * (pb - pa);

          cached = addVertex(p);
          return cached;
        };

        // ------------------------------------------------------------
        // Output cells + metadata
        // ------------------------------------------------------------
        std::vector<std::array<Index, 4>> outCells;
        std::vector<Optional<Attribute>>  outCellAttr;
        std::vector<int8_t>               outCellSide;

        outCells.reserve(static_cast<size_t>(mesh.getCellCount()) * 2);
        outCellAttr.reserve(outCells.capacity());
        outCellSide.reserve(outCells.capacity());

        // Emit a tetrahedron with consistent orientation (positive signed volume).
        auto emitTet = [&](Index a, Index b, Index c, Index d,
                           const Optional<Attribute>& outAttr,
                           int8_t side)
        {
          std::array<Index, 4> t{{a, b, c, d}};
          if (signedVolume(t[0], t[1], t[2], t[3]) < 0)
            std::swap(t[1], t[2]);
          outCells.push_back(t);
          outCellAttr.push_back(outAttr);
          outCellSide.push_back(side);
        };

        // Split templates (standard marching tetrahedra patterns).
        auto split_1neg = [&](Index vN,
                              Index p0, Index p1, Index p2,
                              Index i0, Index i1, Index i2,
                              const Optional<Attribute>& aNeg,
                              const Optional<Attribute>& aPos)
        {
          emitTet(vN, i0, i1, i2, aNeg, SideNegative);
          emitTet(p0, p1, p2, i0, aPos, SidePositive);
          emitTet(p1, p2, i0, i1, aPos, SidePositive);
          emitTet(p2, i0, i1, i2, aPos, SidePositive);
        };

        auto split_1pos = [&](Index vP,
                              Index n0, Index n1, Index n2,
                              Index i0, Index i1, Index i2,
                              const Optional<Attribute>& aNeg,
                              const Optional<Attribute>& aPos)
        {
          emitTet(vP, i0, i1, i2, aPos, SidePositive);
          emitTet(n0, n1, n2, i0, aNeg, SideNegative);
          emitTet(n1, n2, i0, i1, aNeg, SideNegative);
          emitTet(n2, i0, i1, i2, aNeg, SideNegative);
        };

        // Quality proxy: maximize the minimum absolute signed volume among generated tets.
        auto minAbsVol6 = [&](const std::array<std::array<Index, 4>, 6>& tets) -> Real
        {
          Real m = std::numeric_limits<Real>::infinity();
          for (const auto& T : tets)
            m = std::min(m, std::abs(signedVolume(T[0], T[1], T[2], T[3])));
          return m;
        };

        // 2-neg/2-pos split: choose between two valid decompositions based on minAbsVol6.
        auto split_2neg2pos_best = [&](Index n0, Index n1,
                                       Index p0, Index p1,
                                       Index a, Index b, Index c, Index d,
                                       const Optional<Attribute>& aNeg,
                                       const Optional<Attribute>& aPos) -> void
        {
          std::array<std::array<Index, 4>, 6> A = {{
            {{n0, a, b, n1}},
            {{a,  b, n1, c}},
            {{b,  n1, c,  d}},
            {{p0, a, c,  p1}},
            {{a,  c, p1, b}},
            {{c,  p1, b, d}}
          }};

          std::array<std::array<Index, 4>, 6> B = {{
            {{n0, a, b, n1}},
            {{a,  b, n1, d}},
            {{a,  n1, c,  d}},
            {{p0, a, c,  p1}},
            {{a,  c, p1, d}},
            {{a,  p1, b, d}}
          }};

          const Real qa = minAbsVol6(A);
          const Real qb = minAbsVol6(B);
          const auto& best = (qb > qa) ? B : A;

          emitTet(best[0][0], best[0][1], best[0][2], best[0][3], aNeg, SideNegative);
          emitTet(best[1][0], best[1][1], best[1][2], best[1][3], aNeg, SideNegative);
          emitTet(best[2][0], best[2][1], best[2][2], best[2][3], aNeg, SideNegative);

          emitTet(best[3][0], best[3][1], best[3][2], best[3][3], aPos, SidePositive);
          emitTet(best[4][0], best[4][1], best[4][2], best[4][3], aPos, SidePositive);
          emitTet(best[5][0], best[5][1], best[5][2], best[5][3], aPos, SidePositive);
        };

        static constexpr std::array<std::pair<int,int>, 6> TetEdges = {{
          {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
        }};

        // Detect a strict sign change ignoring “zero” vertices (zeros are handled by tie-break).
        auto hasStrictCutNonZero = [&](const int ss[4]) -> bool
        {
          for (const auto& [a,b] : TetEdges)
          {
            if (ss[a] == 0 || ss[b] == 0) continue;
            if (ss[a] * ss[b] == -1) return true;
          }
          return false;
        };

        auto hasAnyZero = [&](const int ss[4]) -> bool
        {
          return (ss[0] == 0) || (ss[1] == 0) || (ss[2] == 0) || (ss[3] == 0);
        };

        // ------------------------------------------------------------
        // Main loop over input tetrahedra
        // ------------------------------------------------------------
        for (auto cit = mesh.getCell(); !cit.end(); ++cit)
        {
          const Index ci = cit->getIndex();
          if (cit->getGeometry() != Polytope::Type::Tetrahedron)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "MarchingTetrahedra only supports tetrahedral meshes."
              << Alert::Raise;
          }

          const auto& cv = cit->getVertices();
          const std::array<Index, 4> v{{ cv(0), cv(1), cv(2), cv(3) }};

          const std::array<Real, 4> phi{{ ls[v[0]], ls[v[1]], ls[v[2]], ls[v[3]] }};
          const auto cellAttr = cit->getAttribute();

          // Non-finite values => do not split, keep cell as unknown.
          if (!finite(phi[0]) || !finite(phi[1]) || !finite(phi[2]) || !finite(phi[3]))
          {
            emitTet(v[0], v[1], v[2], v[3], cellAttr, SideUnknown);
            continue;
          }

          // Adaptive classification tolerance without needing h or grad(phi):
          // eps0 = max(user_eps_sign, rel0 * max(|phi_i|)).
          const Real phi_scale = std::max(std::max(std::abs(phi[0]), std::abs(phi[1])),
                                std::max(std::abs(phi[2]), std::abs(phi[3])));
          const Real rel0 = Real(1e-10);
          const Real eps0 = std::max(eps_sign_user, rel0 * phi_scale);

          // Per-cell snap tolerance is at least eps0.
          const Real eps_snap = std::max(eps_snap_user, eps0);

          // Initial strict sign classification using eps0.
          int ss[4] = {
            (phi[0] < -eps0) ? -1 : (phi[0] > eps0 ? +1 : 0),
            (phi[1] < -eps0) ? -1 : (phi[1] > eps0 ? +1 : 0),
            (phi[2] < -eps0) ? -1 : (phi[2] > eps0 ? +1 : 0),
            (phi[3] < -eps0) ? -1 : (phi[3] > eps0 ? +1 : 0)
          };

          // Deterministic tie-break for near-zero vertices: assign zeros to a stable sign.
          // This prevents “random” patterns of zeros causing topological noise.
          const bool anyZero = hasAnyZero(ss);
          if (anyZero)
          {
            const Real phi_avg = (phi[0] + phi[1] + phi[2] + phi[3]) / Real(4);
            int fallbackSign = 0;

            if (phi_avg > eps0) fallbackSign = +1;
            else if (phi_avg < -eps0) fallbackSign = -1;
            else
            {
              int imax = 0;
              Real amax = std::abs(phi[0]);
              for (int i = 1; i < 4; ++i)
              {
                const Real ai = std::abs(phi[static_cast<size_t>(i)]);
                if (ai > amax) { amax = ai; imax = i; }
              }
              fallbackSign = (phi[static_cast<size_t>(imax)] >= Real(0)) ? +1 : -1;
            }

            for (int i = 0; i < 4; ++i)
              if (ss[i] == 0) ss[i] = fallbackSign;
          }

          const int8_t side0 = sideFromSigns(ss);

          // Decide whether to split geometrically:
          // if no strict cut remains => do not split.
          const bool doSplitGeom = hasStrictCutNonZero(ss);

          if (!doSplitGeom)
          {
            if (side0 == SideNegative)
              emitTet(v[0], v[1], v[2], v[3], inferLabel(3, cellAttr, true), SideNegative);
            else if (side0 == SidePositive)
              emitTet(v[0], v[1], v[2], v[3], inferLabel(3, cellAttr, false), SidePositive);
            else
              emitTet(v[0], v[1], v[2], v[3], cellAttr, SideUnknown);
            continue;
          }

          // Split allowed by SplitMap?
          if (!shouldSplit(3, cellAttr))
          {
            emitTet(v[0], v[1], v[2], v[3], cellAttr, SideUnknown);
            continue;
          }

          // Block splitting if a noSplit edge/face is actually crossed by the interface.
          bool blocked = false;

          const auto& eids = conn.getIncidence({3, 1}, ci);

          // Use the same eps0 as classification to decide cuts for noSplit constraints.
          auto sstrict_cell = [&](Real x) -> int
          {
            if (x < -eps0) return -1;
            if (x >  eps0) return +1;
            return 0;
          };

          for (Index ei : eids)
          {
            const auto ea = mesh.getAttribute(1, ei);
            if (!isNoSplit(1, ea)) continue;

            const auto& ev = conn.getIncidence({1, 0}, ei);
            const int sa = sstrict_cell(ls[ev[0]]);
            const int sb = sstrict_cell(ls[ev[1]]);
            if (isCut(sa, sb)) { blocked = true; break; }
          }

          if (!blocked)
          {
            const auto& fids = conn.getIncidence({3, 2}, ci);
            for (Index fi : fids)
            {
              const auto fa = mesh.getAttribute(2, fi);
              if (!isNoSplit(2, fa)) continue;

              const auto& fv = conn.getIncidence({2, 0}, fi);
              const int s0 = sstrict_cell(ls[fv[0]]);
              const int s1 = sstrict_cell(ls[fv[1]]);
              const int s2 = sstrict_cell(ls[fv[2]]);
              const bool uniform = (s0 == s1 && s1 == s2);
              if (!uniform) { blocked = true; break; }
            }
          }

          if (blocked)
          {
            emitTet(v[0], v[1], v[2], v[3], cellAttr, SideUnknown);
            continue;
          }

          // Determine output labels on each side.
          const auto aNeg = inferLabel(3, cellAttr, true);
          const auto aPos = inferLabel(3, cellAttr, false);

          // Map local edge (ia,ib) to local edge index [0..5] in the incidence list.
          auto edgeIdx = [&](int a, int b) -> int
          {
            if (a > b) std::swap(a, b);
            if (a == 0 && b == 1) return 0;
            if (a == 0 && b == 2) return 1;
            if (a == 0 && b == 3) return 2;
            if (a == 1 && b == 2) return 3;
            if (a == 1 && b == 3) return 4;
            if (a == 2 && b == 3) return 5;
            return -1;
          };

          // Retrieve/create intersection vertex on the corresponding input edge.
          auto I = [&](int ia, int ib) -> Index
          {
            const int e = edgeIdx(ia, ib);
            const Index eid = eids[static_cast<size_t>(e)];
            const Index va = v[static_cast<size_t>(ia)];
            const Index vb = v[static_cast<size_t>(ib)];
            const Real  fa = phi[static_cast<size_t>(ia)];
            const Real  fb = phi[static_cast<size_t>(ib)];
            return getIntersectionOnEdge(eid, va, vb, fa, fb, eps0, eps_snap);
          };

          int nneg = 0;
          for (int i = 0; i < 4; ++i) nneg += (ss[i] < 0);

          // Fully classified after tie-break: handle degenerate cases defensively.
          if (nneg == 0)
          {
            emitTet(v[0], v[1], v[2], v[3], aPos, SidePositive);
            continue;
          }
          if (nneg == 4)
          {
            emitTet(v[0], v[1], v[2], v[3], aNeg, SideNegative);
            continue;
          }

          if (nneg == 1)
          {
            int in = -1; int ip[3]; int k = 0;
            for (int i = 0; i < 4; ++i) (ss[i] < 0) ? (in = i) : (ip[k++] = i);

            const Index i0 = I(in, ip[0]);
            const Index i1 = I(in, ip[1]);
            const Index i2 = I(in, ip[2]);

            split_1neg(v[static_cast<size_t>(in)],
                       v[static_cast<size_t>(ip[0])],
                       v[static_cast<size_t>(ip[1])],
                       v[static_cast<size_t>(ip[2])],
                       i0, i1, i2, aNeg, aPos);
            continue;
          }

          if (nneg == 3)
          {
            int ipos = -1; int ineg[3]; int k = 0;
            for (int i = 0; i < 4; ++i) (ss[i] > 0) ? (ipos = i) : (ineg[k++] = i);

            const Index i0 = I(ipos, ineg[0]);
            const Index i1 = I(ipos, ineg[1]);
            const Index i2 = I(ipos, ineg[2]);

            split_1pos(v[static_cast<size_t>(ipos)],
                       v[static_cast<size_t>(ineg[0])],
                       v[static_cast<size_t>(ineg[1])],
                       v[static_cast<size_t>(ineg[2])],
                       i0, i1, i2, aNeg, aPos);
            continue;
          }

          // nneg == 2
          {
            int ineg[2]; int ipos[2]; int kn = 0, kp = 0;
            for (int i = 0; i < 4; ++i) (ss[i] > 0) ? (ipos[kp++] = i) : (ineg[kn++] = i);

            const Index a = I(ineg[0], ipos[0]);
            const Index b = I(ineg[0], ipos[1]);
            const Index c = I(ineg[1], ipos[0]);
            const Index d = I(ineg[1], ipos[1]);

            split_2neg2pos_best(v[static_cast<size_t>(ineg[0])],
                                v[static_cast<size_t>(ineg[1])],
                                v[static_cast<size_t>(ipos[0])],
                                v[static_cast<size_t>(ipos[1])],
                                a, b, c, d, aNeg, aPos);
            continue;
          }
        }

        // ------------------------------------------------------------
        // Build output mesh (cells)
        // ------------------------------------------------------------
        typename MeshType::Builder builder;
        builder.initialize(3).nodes(outVerts.getCount());
        builder.setVertices(std::move(outVerts));
        builder.reserve(3, outCells.size());

        for (size_t i = 0; i < outCells.size(); ++i)
        {
          const auto& t = outCells[i];
          builder.tetrahedron(IndexArray{{ t[0], t[1], t[2], t[3] }});
          if (outCellAttr[i])
            builder.attribute({3, static_cast<Index>(i)}, *outCellAttr[i]);
        }

        auto out = builder.finalize();

        // ------------------------------------------------------------
        // Attribute transfer / interface marking
        // ------------------------------------------------------------
        out.getConnectivity().compute(2, 3);
        out.getConnectivity().compute(1, 3);
        const auto& oconn = out.getConnectivity();

        const auto& ifaceOpt = this->getInterface();

        // Fitted interface detection uses the *global* sign tolerance to preserve user expectation.
        auto isFittedZeroV = [&](Index vi) -> bool
        {
          return (vi < nv) && (std::abs(ls[vi]) <= eps_sign_user);
        };

        // Faces: mark interface if adjacent to both sides OR fitted to zero;
        // otherwise, transfer attribute only if face is original.
        for (auto fit = out.getPolytope(2); !fit.end(); ++fit)
        {
          if (fit->getGeometry() != Polytope::Type::Triangle)
            continue;

          const Index fidx = fit->getIndex();
          const auto& fv = fit->getVertices();

          const bool allOriginal = (fv(0) < nv) && (fv(1) < nv) && (fv(2) < nv);

          const auto& inc = oconn.getIncidence({2, 3}, fidx);
          bool hasNeg = false, hasPos = false;
          for (const auto ci : inc)
          {
            hasNeg = hasNeg || (outCellSide[ci] == SideNegative);
            hasPos = hasPos || (outCellSide[ci] == SidePositive);
          }

          const bool fittedInterface =
            allOriginal && isFittedZeroV(fv(0)) && isFittedZeroV(fv(1)) && isFittedZeroV(fv(2));

          if ((hasNeg && hasPos) || fittedInterface)
          {
            if (ifaceOpt) out.setAttribute({2, fidx}, *ifaceOpt);
            continue;
          }

          if (!allOriginal)
            continue;

          auto itIn = inFaceAttrByVerts.find(makeFaceKey(fv(0), fv(1), fv(2)));
          if (itIn == inFaceAttrByVerts.end())
            continue;

          const Attribute base = itIn->second;

          // Never carry the old interface tag as a “material” label.
          if (ifaceOpt && base == *ifaceOpt)
            continue;

          const Optional<Attribute> mapped = inferLabel(2, Optional<Attribute>(base), /*neg*/ hasNeg);
          if (mapped)
            out.setAttribute({2, fidx}, *mapped);
        }

        // Edges: same policy as faces.
        for (auto eit = out.getPolytope(1); !eit.end(); ++eit)
        {
          const Index eidx = eit->getIndex();
          const auto& ev = eit->getVertices();

          const bool allOriginal = (ev(0) < nv) && (ev(1) < nv);

          const auto& inc = oconn.getIncidence({1, 3}, eidx);
          bool hasNeg = false, hasPos = false;
          for (const auto ci : inc)
          {
            hasNeg = hasNeg || (outCellSide[ci] == SideNegative);
            hasPos = hasPos || (outCellSide[ci] == SidePositive);
          }

          const bool fittedInterface =
            allOriginal && isFittedZeroV(ev(0)) && isFittedZeroV(ev(1));

          if ((hasNeg && hasPos) || fittedInterface)
          {
            if (ifaceOpt) out.setAttribute({1, eidx}, *ifaceOpt);
            continue;
          }

          if (!allOriginal)
            continue;

          auto itIn = inEdgeAttrByVerts.find(makeEdgeKey(ev(0), ev(1)));
          if (itIn == inEdgeAttrByVerts.end())
            continue;

          const Attribute base = itIn->second;

          if (ifaceOpt && base == *ifaceOpt)
            continue;

          const Optional<Attribute> mapped = inferLabel(1, Optional<Attribute>(base), /*neg*/ hasNeg);
          if (mapped)
            out.setAttribute({1, eidx}, *mapped);
        }

        return out;
      }

      /**
       * @brief Set the user-defined sign tolerance for classification.
       *
       * Values with |phi| <= tol are treated as “near-zero” during sign tests.
       * The discretizer may still apply per-cell adaptation (see discretize()).
       *
       * @param tol Sign tolerance.
       * @return *this
       */
      MarchingTetrahedra& setSignTolerance(Real tol)
      {
        m_sign_tolerance = tol;
        return *this;
      }

      /**
       * @brief Set the user-defined snapping tolerance for edge intersections.
       *
       * This tolerance is used to snap intersections to endpoints / avoid generating
       * tiny sliver cuts when the interface passes extremely close to a vertex.
       *
       * @param tol Snap tolerance.
       * @return *this
       */
      MarchingTetrahedra& setSnapTolerance(Real tol)
      {
        m_snap_tolerance = tol;
        return *this;
      }

      /// @return Current sign tolerance.
      Real getSignTolerance() const
      {
        return m_sign_tolerance;
      }

      /// @return Current snap tolerance.
      Real getSnapTolerance() const
      {
        return m_snap_tolerance;
      }

    private:
      /// User sign tolerance (see setSignTolerance()).
      Real m_sign_tolerance;

      /// User snap tolerance (see setSnapTolerance()).
      Real m_snap_tolerance;
  };

  template <class Mesh, class Data>
  MarchingTetrahedra(const Variational::GridFunction<Variational::P1<Real, Mesh>, Data>& ls)
    -> MarchingTetrahedra<Variational::GridFunction<Variational::P1<Real, Mesh>, Data>>;
}

#endif
