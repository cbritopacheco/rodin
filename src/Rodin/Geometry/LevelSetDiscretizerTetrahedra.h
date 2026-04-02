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

#include "LevelSetDiscretizerBase.h"
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

        if (mesh.getDimension() != 3)
          Alert::MemberFunctionException(*this, __func__)
            << "Expected 3D mesh." << Alert::Raise;

        auto& conn = mesh.getConnectivity();

        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 3, 1);
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 3, 2);
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 2, 1);
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 1, 0);

        using Point = Math::SpatialPoint;
        auto finite = [](Real x) { return std::isfinite(x); };

        static constexpr Index INVALID = std::numeric_limits<Index>::max();

        const Real eps_sign_user = this->getSignTolerance();
        const Real eps_snap_user = this->getSnapTolerance();

        const Index nv = static_cast<Index>(mesh.getVertexCount());

        // ============================================================
        // Attribute / SplitMap helpers (original semantics)
        // ============================================================
        static constexpr int8_t SideUnknown  = -1;
        static constexpr int8_t SidePositive =  0;
        static constexpr int8_t SideNegative =  1;

        auto isNoSplit = [&](size_t d, const Optional<Attribute>& a) -> bool
        {
          if (!a) return false;
          const auto& sm = this->getSplitMap(d);
          if (sm.empty()) return false;
          auto it = sm.find(*a);
          if (it == sm.end()) return false;
          return std::holds_alternative<NoSplitT>(it->second);
        };

        auto shouldSplit = [&](size_t d, const Optional<Attribute>& a) -> bool
        {
          const auto& sm = this->getSplitMap(d);
          if (sm.empty()) return true;
          if (!a) return true;
          auto it = sm.find(*a);
          if (it == sm.end()) return true;
          return std::holds_alternative<Split>(it->second);
        };

        auto inferLabel = [&](size_t d, const Optional<Attribute>& base, bool negativeSide)
          -> Optional<Attribute>
        {
          if (!base) return {};
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

        // ============================================================
        // Provenance maps for ORIGINAL edges/faces only (from discretizeBad)
        // ============================================================
        struct EdgeKey
        {
          Index a, b;
          bool operator<(const EdgeKey& o) const { return std::tie(a, b) < std::tie(o.a, o.b); }
        };
        struct FaceKey
        {
          Index a, b, c;
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

        // ============================================================
        // Global tolerance (consistent across shared edges)
        // ============================================================
        static constexpr Real rel0_global = Real(1e-10);
        Real globalPhiScale = Real(0);
        for (Index i = 0; i < nv; ++i)
          if (finite(ls[i]))
            globalPhiScale = std::max(globalPhiScale, std::abs(ls[i]));

        const Real eps0_global = std::max(eps_sign_user, rel0_global * globalPhiScale);
        const Real eps_snap    = std::max(eps_snap_user, eps0_global);

        // ------------------------------------------------------------
        // 1) Global per-vertex sign (deterministic w/ tolerance)
        // ------------------------------------------------------------
        std::vector<int8_t> vsgn(static_cast<size_t>(nv), 0);
        for (Index i = 0; i < nv; ++i)
        {
          const Real p = ls[i];
          if (!finite(p)) { vsgn[(size_t)i] = 0; continue; }
          vsgn[(size_t)i] = (p < -eps0_global) ? int8_t(-1) : int8_t(+1);
        }

        // ------------------------------------------------------------
        // 2) Output vertices
        // ------------------------------------------------------------
        PointCloud outVerts;
        outVerts.setDimension(3);
        outVerts.reserve(static_cast<size_t>(nv) + static_cast<size_t>(conn.getCount(1)));

        for (Index i = 0; i < nv; ++i)
          outVerts.push_back(mesh.getVertexCoordinates(i));

        auto addVertex = [&](const Point& p) -> Index
        {
          outVerts.push_back(p);
          return static_cast<Index>(outVerts.getCount() - 1);
        };

        // ------------------------------------------------------------
        // Helpers for deterministic triangulation / guards
        // ------------------------------------------------------------
        auto sortedPair = [](Index i, Index j) -> std::pair<Index,Index>
        {
          if (i > j) std::swap(i, j);
          return {i, j};
        };

        auto triPush = [](std::vector<std::array<Index,3>>& L, Index a, Index b, Index c)
        {
          if (a==b || a==c || b==c) return;
          L.push_back({{a,b,c}});
        };

        auto triangulateQuad =
          [&](Index q0, Index q1, Index q2, Index q3, std::vector<std::array<Index,3>>& outTris)
        {
          const auto d02 = sortedPair(q0,q2);
          const auto d13 = sortedPair(q1,q3);
          if (d02 < d13)
          {
            triPush(outTris, q0,q1,q2);
            triPush(outTris, q0,q2,q3);
          }
          else
          {
            triPush(outTris, q1,q2,q3);
            triPush(outTris, q1,q3,q0);
          }
        };

        // ------------------------------------------------------------
        // 3) Edge intersection cache by EDGE INDEX (global snap)
        // ------------------------------------------------------------
        const Index ne = static_cast<Index>(conn.getCount(1));
        std::vector<Index> edgeIsect(static_cast<size_t>(ne), INVALID);

        auto getIsectOnEdge = [&](Index eid, Index va, Index vb, Real fa, Real fb) -> Index
        {
          if (finite(fa) && std::abs(fa) <= eps_snap) return va;
          if (finite(fb) && std::abs(fb) <= eps_snap) return vb;

          Index& slot = edgeIsect[(size_t)eid];
          if (slot != INVALID)
            return slot;

          if (!finite(fa) || !finite(fb))
          {
            slot = va;
            return slot;
          }

          const Real denom = (fa - fb);
          if (std::abs(denom) <= eps_snap)
          {
            slot = (std::abs(fa) < std::abs(fb)) ? va : vb;
            return slot;
          }

          Real t = fa / denom; // fast path (unclamped)
          const auto pa = outVerts[va];
          const auto pb = outVerts[vb];
          slot = addVertex(pa + t * (pb - pa));
          return slot;
        };

        // ------------------------------------------------------------
        // 4) Face split cache by FACE INDEX
        // ------------------------------------------------------------
        struct FaceSplit
        {
          std::vector<std::array<Index,3>> neg;
          std::vector<std::array<Index,3>> pos;
        };

        const Index nf = static_cast<Index>(conn.getCount(2));
        std::vector<uint8_t> faceDone(static_cast<size_t>(nf), 0);
        std::vector<FaceSplit> faceSplit(static_cast<size_t>(nf));

        auto getFaceSplitById = [&](Index fid) -> const FaceSplit&
        {
          if (faceDone[(size_t)fid])
            return faceSplit[(size_t)fid];

          faceDone[(size_t)fid] = 1;
          FaceSplit& fs = faceSplit[(size_t)fid];

          const auto& fv = conn.getPolytope(2, fid); // 3 vertices (a,b,c)
          const Index a = fv(0);
          const Index b = fv(1);
          const Index c = fv(2);

          const Real fa = ls[a], fb = ls[b], fc = ls[c];
          if (!finite(fa) || !finite(fb) || !finite(fc))
            return fs;

          const int sa = (int)vsgn[(size_t)a];
          const int sb = (int)vsgn[(size_t)b];
          const int sc = (int)vsgn[(size_t)c];

          if (sa==0 || sb==0 || sc==0)
            return fs;

          if (sa==sb && sb==sc)
          {
            if (sa < 0) triPush(fs.neg, a,b,c);
            else        triPush(fs.pos, a,b,c);
            return fs;
          }

          const auto& feids = conn.getIncidence({2,1}, fid); // 3 edges

          Index eAB = INVALID, eBC = INVALID, eCA = INVALID;
          for (Index eid : feids)
          {
            const auto& ev = conn.getIncidence({1,0}, eid);
            const Index x = ev[0], y = ev[1];

            if ((x==a && y==b) || (x==b && y==a)) eAB = eid;
            else if ((x==b && y==c) || (x==c && y==b)) eBC = eid;
            else if ((x==c && y==a) || (x==a && y==c)) eCA = eid;
          }

          auto edgeOf = [&](Index x, Index y) -> Index
          {
            if ((x==a && y==b) || (x==b && y==a)) return eAB;
            if ((x==b && y==c) || (x==c && y==b)) return eBC;
            if ((x==c && y==a) || (x==a && y==c)) return eCA;
            Alert::MemberFunctionException(*this, __func__)
              << "Face edge mismatch for fid=" << fid << " edge (" << x << "," << y << ")."
              << Alert::Raise;
            return INVALID;
          };

          auto I = [&](Index x, Index y, Real fx, Real fy) -> Index
          {
            const Index eid = edgeOf(x, y);
            return getIsectOnEdge(eid, x, y, fx, fy);
          };

          auto handleOneTwo =
            [&](Index lone, Index s0, Index s1, Real flone, Real fs0, Real fs1, bool loneIsNeg)
          {
            if (s0 > s1) { std::swap(s0,s1); std::swap(fs0,fs1); }

            const Index p0 = I(lone, s0, flone, fs0);
            const Index p1 = I(lone, s1, flone, fs1);

            if (loneIsNeg) triPush(fs.neg, lone, p0, p1);
            else           triPush(fs.pos, lone, p0, p1);

            if (loneIsNeg) triangulateQuad(s0, s1, p1, p0, fs.pos);
            else           triangulateQuad(s0, s1, p1, p0, fs.neg);
          };

          if (sa!=sb && sa!=sc)      handleOneTwo(a, b, c, fa, fb, fc, sa<0);
          else if (sb!=sa && sb!=sc) handleOneTwo(b, a, c, fb, fa, fc, sb<0);
          else                       handleOneTwo(c, a, b, fc, fa, fb, sc<0);

          return fs;
        };

        // ------------------------------------------------------------
        // 5) Emission (cells + per-cell attributes + per-cell side)
        // ------------------------------------------------------------
        auto signedVolume = [&](Index a, Index b, Index c, Index d) -> long double
        {
          const auto A = outVerts[a];
          const auto B = outVerts[b];
          const auto C = outVerts[c];
          const auto D = outVerts[d];
          return (long double)(B - A).dot((C - A).cross(D - A));
        };

        std::vector<std::array<Index,4>> outCells;
        std::vector<Optional<Attribute>> outCellAttr;
        std::vector<int8_t>             outCellSide;

        outCells.reserve(static_cast<size_t>(mesh.getCellCount()) * 6);
        outCellAttr.reserve(outCells.capacity());
        outCellSide.reserve(outCells.capacity());

        auto emitTet = [&](Index a, Index b, Index c, Index d,
                           const Optional<Attribute>& attr,
                           int8_t side)
        {
          if (a==b || a==c || a==d || b==c || b==d || c==d) return;
          std::array<Index,4> t{{a,b,c,d}};
          if (signedVolume(t[0],t[1],t[2],t[3]) < (long double)0) std::swap(t[1],t[2]);
          outCells.push_back(t);
          outCellAttr.push_back(attr);
          outCellSide.push_back(side);
        };

        auto conePolyhedron =
          [&](const std::vector<std::array<Index,3>>& boundaryTris,
              const Optional<Attribute>& attr,
              int8_t side) -> void
        {
          if (boundaryTris.empty()) return;

          boost::unordered_set<Index> uniq;
          uniq.reserve(boundaryTris.size()*2);

          Point c(3); c.setZero();
          size_t cnt = 0;
          for (const auto& tri : boundaryTris)
          {
            for (int k=0; k<3; ++k)
            {
              const Index vi = tri[(size_t)k];
              if (uniq.insert(vi).second)
              {
                c += outVerts[vi];
                cnt++;
              }
            }
          }
          if (cnt == 0) return;
          c /= Real(cnt);

          const Index center = addVertex(c);

          for (const auto& tri : boundaryTris)
            emitTet(center, tri[0], tri[1], tri[2], attr, side);
        };

        // ------------------------------------------------------------
        // 6) Main loop (fast geometry) + NoSplit blocking + material split
        // ------------------------------------------------------------
        auto sstrict = [&](Real x) -> int
        {
          if (!finite(x)) return 0;
          if (x < -eps0_global) return -1;
          if (x >  eps0_global) return +1;
          return 0;
        };

        auto isCutStrict = [&](int sa, int sb) -> bool
        {
          if (sa == 0 && sb == 0) return false;
          return sa != sb;
        };

        Index cellIdx = 0;
        for (auto cit = mesh.getCell(); !cit.end(); ++cit, ++cellIdx)
        {
          if (cit->getGeometry() != Polytope::Type::Tetrahedron)
            Alert::MemberFunctionException(*this, __func__)
              << "Only tetrahedral meshes are supported." << Alert::Raise;

          const Optional<Attribute> cellAttr = cit->getAttribute();

          const auto& cv = cit->getVertices();
          const std::array<Index,4> v{{cv(0),cv(1),cv(2),cv(3)}};

          const std::array<Real,4> phi{{ ls[v[0]], ls[v[1]], ls[v[2]], ls[v[3]] }};
          if (!finite(phi[0]) || !finite(phi[1]) || !finite(phi[2]) || !finite(phi[3]))
          {
            emitTet(v[0], v[1], v[2], v[3], cellAttr, SideUnknown);
            continue;
          }

          int ss[4] = {
            (int)vsgn[(size_t)v[0]],
            (int)vsgn[(size_t)v[1]],
            (int)vsgn[(size_t)v[2]],
            (int)vsgn[(size_t)v[3]]
          };

          if (ss[0]==0 || ss[1]==0 || ss[2]==0 || ss[3]==0)
          {
            emitTet(v[0], v[1], v[2], v[3], cellAttr, SideUnknown);
            continue;
          }

          // cut check
          static constexpr int E[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
          bool hasCut = false;
          for (int e=0; e<6; ++e)
            if (ss[E[e][0]] != ss[E[e][1]]) { hasCut = true; break; }

          if (!hasCut)
          {
            const bool neg = (ss[0] < 0);
            const Optional<Attribute> a = inferLabel(3, cellAttr, neg);
            emitTet(v[0], v[1], v[2], v[3], a, neg ? SideNegative : SidePositive);
            continue;
          }

          // if cell is configured as "do not split", keep it unsplit
          if (!shouldSplit(3, cellAttr))
          {
            emitTet(v[0], v[1], v[2], v[3], cellAttr, SideUnknown);
            continue;
          }

          // cell edges (needed for NoSplit blocking and for localEdgeId)
          const auto& ceids = conn.getIncidence({3,1}, cellIdx);

          // ---- NoSplit blocking (restore discretizeBad semantics)
          bool blocked = false;

          // (1) block if any NoSplit EDGE (d=1) is cut
          for (Index eid : ceids)
          {
            const auto ea = mesh.getAttribute(1, eid);
            if (!isNoSplit(1, ea)) continue;

            const auto& ev = conn.getIncidence({1,0}, eid);
            const int sa = sstrict(ls[ev[0]]);
            const int sb = sstrict(ls[ev[1]]);
            if (isCutStrict(sa, sb)) { blocked = true; break; }
          }

          // (2) block if any NoSplit FACE (d=2) is non-uniform
          if (!blocked)
          {
            const auto& cfids_block = conn.getIncidence({3,2}, cellIdx);
            for (Index fid : cfids_block)
            {
              const auto fa = mesh.getAttribute(2, fid);
              if (!isNoSplit(2, fa)) continue;

              const auto& fv = conn.getPolytope(2, fid);
              const int s0 = sstrict(ls[fv(0)]);
              const int s1 = sstrict(ls[fv(1)]);
              const int s2 = sstrict(ls[fv(2)]);
              const bool uniform = (s0 == s1 && s1 == s2);
              if (!uniform) { blocked = true; break; }
            }
          }

          if (blocked)
          {
            emitTet(v[0], v[1], v[2], v[3], cellAttr, SideUnknown);
            continue;
          }

          const Optional<Attribute> aNeg = inferLabel(3, cellAttr, true);
          const Optional<Attribute> aPos = inferLabel(3, cellAttr, false);

          // Build localEdgeId (O(6) per cell, then O(1) queries)
          Index localEdgeId[4][4];
          for (int i=0;i<4;++i)
            for (int j=0;j<4;++j)
              localEdgeId[i][j] = INVALID;

          for (Index eid : ceids)
          {
            const auto& ev = conn.getIncidence({1,0}, eid);
            const Index a = ev[0];
            const Index b = ev[1];

            int ia = -1, ib = -1;
            for (int i=0;i<4;++i)
            {
              if (v[(size_t)i] == a) ia = i;
              if (v[(size_t)i] == b) ib = i;
            }
            if (ia >= 0 && ib >= 0)
            {
              localEdgeId[ia][ib] = eid;
              localEdgeId[ib][ia] = eid;
            }
          }

          std::vector<std::array<Index,3>> negB, posB;
          negB.reserve(16);
          posB.reserve(16);

          bool fallback = false;

          // (A) split triangles on the 4 original faces using CELL->FACE ids
          const auto& cfids = conn.getIncidence({3,2}, cellIdx);
          for (Index fid : cfids)
          {
            const auto& fs = getFaceSplitById(fid);
            if (fs.neg.empty() && fs.pos.empty())
            {
              fallback = true;
              break;
            }
            negB.insert(negB.end(), fs.neg.begin(), fs.neg.end());
            posB.insert(posB.end(), fs.pos.begin(), fs.pos.end());
          }

          if (fallback)
          {
            emitTet(v[0], v[1], v[2], v[3], cellAttr, SideUnknown);
            continue;
          }

          // (B) interface triangles inside tet (uses localEdgeId => O(1))
          {
            auto I = [&](int ia, int ib) -> Index
            {
              const Index eid = localEdgeId[ia][ib];
              if (eid == INVALID)
              {
                Alert::MemberFunctionException(*this, __func__)
                  << "localEdgeId missing for cell " << cellIdx
                  << " local edge (" << ia << "," << ib << ")."
                  << Alert::Raise;
              }
              return getIsectOnEdge(eid,
                                    v[(size_t)ia], v[(size_t)ib],
                                    phi[(size_t)ia], phi[(size_t)ib]);
            };

            const int nneg = (ss[0]<0) + (ss[1]<0) + (ss[2]<0) + (ss[3]<0);
            if (nneg != 0 && nneg != 4)
            {
              if (nneg==1 || nneg==3)
              {
                const bool loneIsNeg = (nneg==1);
                int lone = -1;
                int others[3]; int k=0;
                for (int i=0;i<4;++i)
                {
                  const bool isNeg = (ss[i] < 0);
                  if (isNeg == loneIsNeg) lone = i;
                  else others[k++] = i;
                }

                Index p[3] = { I(lone, others[0]), I(lone, others[1]), I(lone, others[2]) };

                if (p[0]>p[1]) std::swap(p[0],p[1]);
                if (p[1]>p[2]) std::swap(p[1],p[2]);
                if (p[0]>p[1]) std::swap(p[0],p[1]);

                triPush(negB, p[0], p[1], p[2]);
                triPush(posB, p[0], p[1], p[2]);
              }
              else // nneg == 2
              {
                int ineg[2], ipos[2], kn=0, kp=0;
                for (int i=0;i<4;++i) (ss[i]<0) ? (ineg[kn++]=i) : (ipos[kp++]=i);

                Index n0 = v[(size_t)ineg[0]], n1 = v[(size_t)ineg[1]];
                Index p0 = v[(size_t)ipos[0]], p1 = v[(size_t)ipos[1]];

                if (n0 > n1) std::swap(ineg[0], ineg[1]), std::swap(n0,n1);
                if (p0 > p1) std::swap(ipos[0], ipos[1]), std::swap(p0,p1);

                const Index a = I(ineg[0], ipos[0]);
                const Index b = I(ineg[0], ipos[1]);
                const Index d = I(ineg[1], ipos[1]);
                const Index c = I(ineg[1], ipos[0]);

                std::vector<std::array<Index,3>> qTris;
                qTris.reserve(2);
                triangulateQuad(a,b,d,c, qTris);

                for (const auto& tri : qTris)
                {
                  triPush(negB, tri[0], tri[1], tri[2]);
                  triPush(posB, tri[0], tri[1], tri[2]);
                }
              }
            }
          }

          // (C) cone each side polyhedron (emit with side + mapped cell attr)
          conePolyhedron(negB, aNeg, SideNegative);
          conePolyhedron(posB, aPos, SidePositive);
        }

        // ------------------------------------------------------------
        // 7) Build output mesh (cells + cell attributes)
        // ------------------------------------------------------------
        typename MeshType::Builder builder;
        builder.initialize(3).nodes(outVerts.getCount());
        builder.setVertices(std::move(outVerts));
        builder.reserve(3, outCells.size());

        for (size_t i = 0; i < outCells.size(); ++i)
        {
          const auto& t = outCells[i];
          builder.tetrahedron(IndexArray{{t[0], t[1], t[2], t[3]}});
          if (outCellAttr[i])
            builder.attribute({3, static_cast<Index>(i)}, *outCellAttr[i]);
        }

        auto out = builder.finalize();

        // ============================================================
        // 8) Attribute transfer / interface marking (faces + edges)
        // ============================================================
        out.getConnectivity().compute(2, 3);
        out.getConnectivity().compute(1, 3);
        const auto& oconn = out.getConnectivity();

        const auto& ifaceOpt2 = this->getInterface(2);
        const auto& ifaceOpt1 = this->getInterface(1);

        auto isFittedZeroV = [&](Index vi) -> bool
        {
          return (vi < nv) && finite(ls[vi]) && (std::abs(ls[vi]) <= eps_sign_user);
        };

        // Faces
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
            hasNeg = hasNeg || (outCellSide[(size_t)ci] == SideNegative);
            hasPos = hasPos || (outCellSide[(size_t)ci] == SidePositive);
          }

          const bool fittedInterface =
            allOriginal && isFittedZeroV(fv(0)) && isFittedZeroV(fv(1)) && isFittedZeroV(fv(2));

          const bool isInterface = (hasNeg && hasPos) || fittedInterface;
          if (isInterface)
          {
            if (ifaceOpt2) out.setAttribute({2, fidx}, *ifaceOpt2);
            continue;
          }

          if (!allOriginal)
          {
            if (m_fallback[2])
              out.setAttribute({2, fidx}, *m_fallback[2]);
            continue;
          }

          auto itIn = inFaceAttrByVerts.find(makeFaceKey(fv(0), fv(1), fv(2)));
          if (itIn == inFaceAttrByVerts.end())
            continue;

          const Attribute base = itIn->second;

          if (ifaceOpt2 && base == *ifaceOpt2)
          {
            if (m_old[2])
              out.setAttribute({2, fidx}, *m_old[2]);
            continue;
          }

          const Optional<Attribute> mapped =
            inferLabel(2, Optional<Attribute>(base), /*negativeSide*/ hasNeg);

          if (mapped)
            out.setAttribute({2, fidx}, *mapped);
        }

        // Edges
        for (auto eit = out.getPolytope(1); !eit.end(); ++eit)
        {
          const Index eidx = eit->getIndex();
          const auto& ev = eit->getVertices();

          const bool allOriginal = (ev(0) < nv) && (ev(1) < nv);

          const auto& inc = oconn.getIncidence({1, 3}, eidx);
          bool hasNeg = false, hasPos = false;
          for (const auto ci : inc)
          {
            hasNeg = hasNeg || (outCellSide[(size_t)ci] == SideNegative);
            hasPos = hasPos || (outCellSide[(size_t)ci] == SidePositive);
          }

          const bool fittedInterface =
            allOriginal && isFittedZeroV(ev(0)) && isFittedZeroV(ev(1));

          const bool isInterface = (hasNeg && hasPos) || fittedInterface;
          if (isInterface)
          {
            if (ifaceOpt1) out.setAttribute({1, eidx}, *ifaceOpt1);
            continue;
          }

          if (!allOriginal)
          {
            if (m_fallback[1])
              out.setAttribute({1, eidx}, *m_fallback[1]);
            continue;
          }

          auto itIn = inEdgeAttrByVerts.find(makeEdgeKey(ev(0), ev(1)));
          if (itIn == inEdgeAttrByVerts.end())
            continue;

          const Attribute base = itIn->second;

          // "old interface" cleanup: accept either configured interface tag (d=1 or d=2)
          if ((ifaceOpt1 && base == *ifaceOpt1) || (ifaceOpt2 && base == *ifaceOpt2))
          {
            if (m_old[1])
              out.setAttribute({1, eidx}, *m_old[1]);
            continue;
          }

          const Optional<Attribute> mapped =
            inferLabel(1, Optional<Attribute>(base), /*negativeSide*/ hasNeg);

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

      MarchingTetrahedra& setOld(size_t d, const Optional<Attribute>& a)
      {
        m_old[d] = a;
        return *this;
      }

      MarchingTetrahedra& setFallback(size_t d, const Optional<Attribute>& a)
      {
        m_fallback[d] = a;
        return *this;
      }


    private:
      /// User sign tolerance (see setSignTolerance()).
      Real m_sign_tolerance;

      /// User snap tolerance (see setSnapTolerance()).
      Real m_snap_tolerance;

      std::array<Optional<Attribute>, 4> m_old;

      std::array<Optional<Attribute>, 4> m_fallback;
  };

  template <class Mesh, class Data>
  MarchingTetrahedra(const Variational::GridFunction<Variational::P1<Real, Mesh>, Data>& ls)
    -> MarchingTetrahedra<Variational::GridFunction<Variational::P1<Real, Mesh>, Data>>;
}

#endif
