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
  template <class ... Params>
  class MarchingTetrahedra;

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

      MarchingTetrahedra(const GridFunctionType& ls)
        : Parent(ls)
      {}

      MeshType discretize() const override
      {
        const auto& ls   = this->getGridFunction();
        const auto& mesh = ls.getFiniteElementSpace().getMesh();

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

        const Real eps = this->getTolerance();

        enum class Sign : uint8_t { Negative, Positive, Zero };

        auto sgn = [&](Real x) -> Sign
        {
          if (x < -eps) return Sign::Negative;
          if (x >  eps) return Sign::Positive;
          return Sign::Zero;
        };

        auto finite = [](Real x) -> bool { return std::isfinite(x); };

        // ---- SplitMap helpers
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
          if (!a) return true; // no attribute => default behavior: split
          const auto& sm = this->getSplitMap(d);
          if (sm.empty()) return true;      // default: split everything
          auto it = sm.find(*a);
          if (it == sm.end()) return false; // rule missing => do not split for this attr
          return std::holds_alternative<Split>(it->second);
        };

        // Labeling policy:
        // - If base attribute is missing => use global +/- if provided.
        // - If base attribute exists and there is an explicit split rule => apply it.
        // - Otherwise KEEP base attribute (do not override with global +/-).
        auto splitLabel = [&](size_t d, const Optional<Attribute>& base, bool negativeSide)
          -> Optional<Attribute>
        {
          if (!base)
          {
            const auto& want = negativeSide ? this->getNegative() : this->getPositive();
            return want; // may be empty
          }

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

        // ---- Signed volume (6x tetra volume)
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

        // ---------------------------------------------------------------------
        // Input lookup maps (by vertex sets) to RETAIN original EDGE/FACE attrs.
        // ---------------------------------------------------------------------
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

        for (auto eit = mesh.getPolytope(1); !eit.end(); ++eit)
        {
          const auto& ev = eit->getVertices();
          const auto a = eit->getAttribute();
          if (a)
            inEdgeAttrByVerts[makeEdgeKey(ev(0), ev(1))] = *a;
        }

        for (auto fit = mesh.getPolytope(2); !fit.end(); ++fit)
        {
          const auto g = fit->getGeometry();
          if (g != Polytope::Type::Triangle)
            continue;
          const auto& fv = fit->getVertices();
          const auto a = fit->getAttribute();
          if (a)
            inFaceAttrByVerts[makeFaceKey(fv(0), fv(1), fv(2))] = *a;
        }

        // ---- Output vertices
        PointCloud outVerts;
        outVerts.setDimension(3);
        outVerts.reserve(static_cast<size_t>(nv) + static_cast<size_t>(ne / 4));

        for (Index i = 0; i < nv; ++i)
          outVerts.push_back(mesh.getVertexCoordinates(i));

        auto addVertex = [&](const Point& p) -> Index
        {
          outVerts.push_back(p);
          return static_cast<Index>(outVerts.getCount() - 1);
        };

        // ---- Edge intersection cache by global edge id
        std::vector<Index> isect(static_cast<size_t>(ne), InvalidIndex);

        auto getIntersectionOnEdge = [&](Index edgeIdx,
                                         Index va, Index vb,
                                         Real fa, Real fb) -> Index
        {
          // Prevent NaN propagation to geometry
          if (!finite(fa) || !finite(fb))
            return va;

          if (std::abs(fa) <= eps) return va;
          if (std::abs(fb) <= eps) return vb;

          Index& cached = isect[static_cast<size_t>(edgeIdx)];
          if (cached != InvalidIndex)
            return cached;

          const Real denom = (fa - fb);

          if (std::abs(denom) <= eps)
            return (std::abs(fa) < std::abs(fb)) ? va : vb;

          Real t = fa / denom;
          t = std::min<Real>(Real(1), std::max<Real>(Real(0), t));

          if (t <= eps)            return va;
          if (t >= Real(1) - eps)  return vb;

          const auto pa = outVerts[va];
          const auto pb = outVerts[vb];
          const auto p  = pa + t * (pb - pa);

          cached = addVertex(p);
          return cached;
        };

        // ---- Output cells + metadata (needed to label faces/edges)
        static constexpr int8_t SideUnknown  = -1;
        static constexpr int8_t SidePositive = 0;
        static constexpr int8_t SideNegative = 1;

        std::vector<std::array<Index, 4>> outCells;
        std::vector<Optional<Attribute>> outCellAttr;
        std::vector<Optional<Attribute>> outCellSourceAttr;
        std::vector<int8_t> outCellSide;

        outCells.reserve(static_cast<size_t>(mesh.getCellCount()) * 2);
        outCellAttr.reserve(outCells.capacity());
        outCellSourceAttr.reserve(outCells.capacity());
        outCellSide.reserve(outCells.capacity());

        auto emitTet = [&](Index a, Index b, Index c, Index d,
                           const Optional<Attribute>& outAttr,
                           const Optional<Attribute>& sourceAttr,
                           int8_t side)
        {
          std::array<Index, 4> t{{a, b, c, d}};
          if (signedVolume(t[0], t[1], t[2], t[3]) < 0)
            std::swap(t[1], t[2]);

          outCells.push_back(t);
          outCellAttr.push_back(outAttr);
          outCellSourceAttr.push_back(sourceAttr);
          outCellSide.push_back(side);
        };

        // ---- Templates
        auto split_1neg = [&](Index vN,
                              Index p0, Index p1, Index p2,
                              Index i0, Index i1, Index i2,
                              const Optional<Attribute>& aNeg,
                              const Optional<Attribute>& aPos,
                              const Optional<Attribute>& sourceAttr)
        {
          emitTet(vN, i0, i1, i2, aNeg, sourceAttr, SideNegative);

          emitTet(p0, p1, p2, i0, aPos, sourceAttr, SidePositive);
          emitTet(p1, p2, i0, i1, aPos, sourceAttr, SidePositive);
          emitTet(p2, i0, i1, i2, aPos, sourceAttr, SidePositive);
        };

        auto split_1pos = [&](Index vP,
                              Index n0, Index n1, Index n2,
                              Index i0, Index i1, Index i2,
                              const Optional<Attribute>& aNeg,
                              const Optional<Attribute>& aPos,
                              const Optional<Attribute>& sourceAttr)
        {
          emitTet(vP, i0, i1, i2, aPos, sourceAttr, SidePositive);

          emitTet(n0, n1, n2, i0, aNeg, sourceAttr, SideNegative);
          emitTet(n1, n2, i0, i1, aNeg, sourceAttr, SideNegative);
          emitTet(n2, i0, i1, i2, aNeg, sourceAttr, SideNegative);
        };

        auto minAbsVol6 = [&](const std::array<std::array<Index, 4>, 6>& tets) -> Real
        {
          Real m = std::numeric_limits<Real>::infinity();
          for (const auto& T : tets)
            m = std::min(m, std::abs(signedVolume(T[0], T[1], T[2], T[3])));
          return m;
        };

        auto split_2neg2pos_best = [&](Index n0, Index n1,
                                       Index p0, Index p1,
                                       Index a, Index b, Index c, Index d,
                                       const Optional<Attribute>& aNeg,
                                       const Optional<Attribute>& aPos,
                                       const Optional<Attribute>& sourceAttr)
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

          emitTet(best[0][0], best[0][1], best[0][2], best[0][3], aNeg, sourceAttr, SideNegative);
          emitTet(best[1][0], best[1][1], best[1][2], best[1][3], aNeg, sourceAttr, SideNegative);
          emitTet(best[2][0], best[2][1], best[2][2], best[2][3], aNeg, sourceAttr, SideNegative);

          emitTet(best[3][0], best[3][1], best[3][2], best[3][3], aPos, sourceAttr, SidePositive);
          emitTet(best[4][0], best[4][1], best[4][2], best[4][3], aPos, sourceAttr, SidePositive);
          emitTet(best[5][0], best[5][1], best[5][2], best[5][3], aPos, sourceAttr, SidePositive);
        };

        // ---- Main loop
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
          if (!finite(phi[0]) || !finite(phi[1]) || !finite(phi[2]) || !finite(phi[3]))
          {
            const auto cellAttr = cit->getAttribute();
            emitTet(v[0], v[1], v[2], v[3], cellAttr, cellAttr, SideUnknown);
            continue;
          }

          const std::array<Sign, 4> s{{ sgn(phi[0]), sgn(phi[1]), sgn(phi[2]), sgn(phi[3]) }};
          const auto cellAttr = cit->getAttribute();

          int nneg = 0;
          for (int i = 0; i < 4; ++i)
            nneg += (s[i] != Sign::Positive); // Negative OR Zero treated as "neg side"

          if (nneg == 0)
          {
            emitTet(v[0], v[1], v[2], v[3],
                    splitLabel(3, cellAttr, /*negative*/ false),
                    cellAttr, SidePositive);
            continue;
          }

          if (nneg == 4)
          {
            emitTet(v[0], v[1], v[2], v[3],
                    splitLabel(3, cellAttr, /*negative*/ true),
                    cellAttr, SideNegative);
            continue;
          }

          if (!shouldSplit(3, cellAttr))
          {
            emitTet(v[0], v[1], v[2], v[3], cellAttr, cellAttr, SideUnknown);
            continue;
          }

          // Block by noSplit edges/faces
          bool blocked = false;

          // RELIES ON DOCUMENTED INVARIANT:
          // conn.getIncidence({3,1}, ci) ordering for tets is:
          // 0:(0,1) 1:(0,2) 2:(0,3) 3:(1,2) 4:(1,3) 5:(2,3)
          const auto& eids = conn.getIncidence({3, 1}, ci);

          for (Index ei : eids)
          {
            const auto ea = mesh.getAttribute(1, ei);
            if (!isNoSplit(1, ea))
              continue;

            const auto& ev = conn.getIncidence({1, 0}, ei);
            const Sign sa = sgn(ls[ev[0]]);
            const Sign sb = sgn(ls[ev[1]]);

            const bool cut =
              (sa == Sign::Negative && sb == Sign::Positive) || (sa == Sign::Positive && sb == Sign::Negative) ||
              (sa == Sign::Zero && sb != Sign::Zero) || (sb == Sign::Zero && sa != Sign::Zero);

            if (cut) { blocked = true; break; }
          }

          if (!blocked)
          {
            const auto& fids = conn.getIncidence({3, 2}, ci);
            for (Index fi : fids)
            {
              const auto fa = mesh.getAttribute(2, fi);
              if (!isNoSplit(2, fa))
                continue;

              const auto& fv = conn.getIncidence({2, 0}, fi);
              const Sign s0 = sgn(ls[fv[0]]);
              const Sign s1 = sgn(ls[fv[1]]);
              const Sign s2 = sgn(ls[fv[2]]);

              const bool uniform = (s0 == s1 && s1 == s2);
              if (!uniform) { blocked = true; break; }
            }
          }

          if (blocked)
          {
            emitTet(v[0], v[1], v[2], v[3], cellAttr, cellAttr, SideUnknown);
            continue;
          }

          const auto aNeg = splitLabel(3, cellAttr, /*negative*/ true);
          const auto aPos = splitLabel(3, cellAttr, /*negative*/ false);

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

          auto I = [&](int ia, int ib) -> Index
          {
            const int e = edgeIdx(ia, ib);
            const Index eid = eids[static_cast<size_t>(e)];
            const Index va = v[static_cast<size_t>(ia)];
            const Index vb = v[static_cast<size_t>(ib)];
            const Real  fa = phi[static_cast<size_t>(ia)];
            const Real  fb = phi[static_cast<size_t>(ib)];
            return getIntersectionOnEdge(eid, va, vb, fa, fb);
          };

          if (nneg == 1)
          {
            int in = -1; int ip[3]; int k = 0;
            for (int i = 0; i < 4; ++i)
              (s[i] != Sign::Positive) ? (in = i) : (ip[k++] = i); // Zero with negative

            const Index i0 = I(in, ip[0]);
            const Index i1 = I(in, ip[1]);
            const Index i2 = I(in, ip[2]);

            split_1neg(v[static_cast<size_t>(in)],
                       v[static_cast<size_t>(ip[0])],
                       v[static_cast<size_t>(ip[1])],
                       v[static_cast<size_t>(ip[2])],
                       i0, i1, i2, aNeg, aPos, cellAttr);
            continue;
          }

          if (nneg == 3)
          {
            int ipos = -1; int ineg[3]; int k = 0;
            for (int i = 0; i < 4; ++i)
              (s[i] == Sign::Positive) ? (ipos = i) : (ineg[k++] = i); // Zero with negative

            const Index i0 = I(ipos, ineg[0]);
            const Index i1 = I(ipos, ineg[1]);
            const Index i2 = I(ipos, ineg[2]);

            split_1pos(v[static_cast<size_t>(ipos)],
                       v[static_cast<size_t>(ineg[0])],
                       v[static_cast<size_t>(ineg[1])],
                       v[static_cast<size_t>(ineg[2])],
                       i0, i1, i2, aNeg, aPos, cellAttr);
            continue;
          }

          // nneg == 2
          {
            int ineg[2]; int ipos[2]; int kn = 0, kp = 0;
            for (int i = 0; i < 4; ++i)
              (s[i] == Sign::Positive) ? (ipos[kp++] = i) : (ineg[kn++] = i); // Zero with negative

            const Index a = I(ineg[0], ipos[0]);
            const Index b = I(ineg[0], ipos[1]);
            const Index c = I(ineg[1], ipos[0]);
            const Index d = I(ineg[1], ipos[1]);

            split_2neg2pos_best(v[static_cast<size_t>(ineg[0])],
                                v[static_cast<size_t>(ineg[1])],
                                v[static_cast<size_t>(ipos[0])],
                                v[static_cast<size_t>(ipos[1])],
                                a, b, c, d, aNeg, aPos, cellAttr);
            continue;
          }
        }

        // ---- Build output mesh (cells)
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

        // ---------------------------------------------------------------------
        // Transfer FACE/EDGE attributes using input maps.
        // IMPORTANT: do NOT inherit old interface tags (Gamma) from input.
        // Interface is recomputed from adjacency (hasNeg && hasPos) each time.
        // ---------------------------------------------------------------------
        out.getConnectivity().compute(2, 3);
        out.getConnectivity().compute(1, 3);
        const auto& oconn = out.getConnectivity();

        const auto& ifaceOpt = this->getInterface();

        // Faces (triangles only) + interface marking
        for (auto fit = out.getPolytope(2); !fit.end(); ++fit)
        {
          if (fit->getGeometry() != Polytope::Type::Triangle)
            continue;

          const Index fidx = fit->getIndex();
          const auto& inc = oconn.getIncidence({2, 3}, fidx);
          if (inc.empty())
            continue;

          bool hasNeg = false, hasPos = false, hasUnknown = false;
          for (const auto ci : inc)
          {
            hasNeg     = hasNeg     || (outCellSide[ci] == SideNegative);
            hasPos     = hasPos     || (outCellSide[ci] == SidePositive);
            hasUnknown = hasUnknown || (outCellSide[ci] == SideUnknown);
          }

          if (hasNeg && hasPos)
          {
            if (ifaceOpt)
              out.setAttribute({2, fidx}, *ifaceOpt);
            continue;
          }

          const auto& fv = fit->getVertices();
          const bool allOriginal = (fv(0) < nv) && (fv(1) < nv) && (fv(2) < nv);
          if (!allOriginal)
            continue;

          auto itIn = inFaceAttrByVerts.find(makeFaceKey(fv(0), fv(1), fv(2)));
          if (itIn == inFaceAttrByVerts.end())
            continue;

          const Attribute baseFaceAttr = itIn->second;

          // Do not carry over old interface tag
          if (ifaceOpt && baseFaceAttr == *ifaceOpt)
            continue;

          if (hasUnknown)
            out.setAttribute({2, fidx}, baseFaceAttr);
          else
            out.setAttribute({2, fidx}, splitLabel(2, baseFaceAttr, /*negative*/ hasNeg));
        }

        // Edges + interface marking
        for (auto eit = out.getPolytope(1); !eit.end(); ++eit)
        {
          const Index eidx = eit->getIndex();
          const auto& inc = oconn.getIncidence({1, 3}, eidx);
          if (inc.empty())
            continue;

          bool hasNeg = false, hasPos = false, hasUnknown = false;
          for (const auto ci : inc)
          {
            hasNeg     = hasNeg     || (outCellSide[ci] == SideNegative);
            hasPos     = hasPos     || (outCellSide[ci] == SidePositive);
            hasUnknown = hasUnknown || (outCellSide[ci] == SideUnknown);
          }

          if (hasNeg && hasPos)
          {
            if (ifaceOpt)
              out.setAttribute({1, eidx}, *ifaceOpt);
            continue;
          }

          const auto& ev = eit->getVertices();
          const bool allOriginal = (ev(0) < nv) && (ev(1) < nv);
          if (!allOriginal)
            continue;

          auto itIn = inEdgeAttrByVerts.find(makeEdgeKey(ev(0), ev(1)));
          if (itIn == inEdgeAttrByVerts.end())
            continue;

          const Attribute baseEdgeAttr = itIn->second;

          // Do not carry over old interface tag
          if (ifaceOpt && baseEdgeAttr == *ifaceOpt)
            continue;

          if (hasUnknown)
            out.setAttribute({1, eidx}, baseEdgeAttr);
          else
            out.setAttribute({1, eidx}, splitLabel(1, baseEdgeAttr, /*negative*/ hasNeg));
        }

        return out;
      }
  };

  template <class Mesh, class Data>
  MarchingTetrahedra(const Variational::GridFunction<Variational::P1<Real, Mesh>, Data>& ls)
    -> MarchingTetrahedra<Variational::GridFunction<Variational::P1<Real, Mesh>, Data>>;
}

#endif
