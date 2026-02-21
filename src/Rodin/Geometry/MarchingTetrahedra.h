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
#include <map>
#include <optional>
#include <tuple>
#include <type_traits>
#include <variant>
#include <vector>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Alert/Warning.h"

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

      using FESType         = typename Parent::FESType;
      using MeshType        = typename Parent::MeshType;
      using GridFunctionType= typename Parent::GridFunctionType;

      using SplitMap        = typename Parent::SplitMap;
      using NoSplitT        = typename Parent::NoSplitT;
      using Split           = typename Parent::Split;

      MarchingTetrahedra(const GridFunctionType& ls)
        : Parent(ls)
      {}

      MeshType discretize() const override
      {
        const auto& ls   = this->getGridFunction();
        const auto& mesh = ls.getFiniteElementSpace().getMesh();

        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 3, 2);
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 3, 1);

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

        // ---- Robust sign classification
        const Real eps = this->getTolerance();
        enum class Sign : uint8_t { Negative, Positive, Zero };

        auto sgn = [&](Real x) -> Sign
        {
          if (x < -eps) return Sign::Negative;
          if (x >  eps) return Sign::Positive;
          return Sign::Zero;
        };

        // ---- SplitMap helpers
        auto isNoSplit = [&](size_t d, Attribute a) -> bool
        {
          const auto& sm = this->getSplitMap(d);
          if (sm.empty()) return false;
          auto it = sm.find(a);
          if (it == sm.end()) return false;
          return std::holds_alternative<NoSplitT>(it->second);
        };

        auto shouldSplit = [&](size_t d, Attribute a) -> bool
        {
          const auto& sm = this->getSplitMap(d);
          if (sm.empty()) return true;      // default: split everything
          auto it = sm.find(a);
          if (it == sm.end()) return false; // default: do not split if rule missing
          return std::holds_alternative<Split>(it->second);
        };

        auto splitLabel = [&](size_t d, Attribute base, bool negativeSide) -> Attribute
        {
          const auto def = [&]()
          {
            return negativeSide ? this->getNegative() : this->getPositive();
          };

          const auto& sm = this->getSplitMap(d);
          if (sm.empty())
            return def();

          auto it = sm.find(base);
          if (it == sm.end())
            return base;

          return std::visit(
            [&](const auto& v) -> Attribute
            {
              using T = std::decay_t<decltype(v)>;
              if constexpr (std::is_same_v<T, NoSplitT>)
                return base;
              else
                return negativeSide ? v.negative : v.positive;
            },
            it->second);
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

        // ---- Characteristic edge length (mesh scale)
        // Uses eps everywhere small thresholds are needed (no hard-coded 1e-12/1e-14).
        auto characteristicEdgeLength = [&]() -> Real
        {
          if (ne <= 0)
            return Real(1);

          Real sum = 0;
          Index cnt = 0;

          // assumes you can iterate over 1D polytopes (edges)
          for (auto eit = mesh.getPolytope(1); !eit.end(); ++eit)
          {
            const Index ei = eit->getIndex();
            const auto& ev = conn.getIncidence({1, 0}, ei);
            const auto p0 = mesh.getVertexCoordinates(ev[0]);
            const auto p1 = mesh.getVertexCoordinates(ev[1]);
            sum += (p1 - p0).norm();
            ++cnt;
          }

          return (cnt > 0) ? (sum / cnt) : Real(1);
        };

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
          // Snap to existing mesh vertices when near interface
          if (std::abs(fa) <= eps) return va;
          if (std::abs(fb) <= eps) return vb;

          Index& cached = isect[static_cast<size_t>(edgeIdx)];
          if (cached != InvalidIndex)
            return cached;

          const Real denom = (fa - fb);

          // Near-constant along edge: choose the closer endpoint
          if (std::abs(denom) <= eps)
            return (std::abs(fa) < std::abs(fb)) ? va : vb;

          // Linear interpolation
          Real t = fa / denom;

          // Clamp (roundoff safety)
          t = std::min<Real>(Real(1), std::max<Real>(Real(0), t));

          // Snap to endpoint if extremely close (uses eps; t is unitless, but eps is your only knob here)
          if (t <= eps)            return va;
          if (t >= Real(1) - eps)  return vb;

          const auto pa = outVerts[va];
          const auto pb = outVerts[vb];
          const auto p  = pa + t * (pb - pa);

          cached = addVertex(p);
          return cached;
        };

        // ---- Output cells + attributes
        std::vector<std::array<Index, 4>> outCells;
        std::vector<Attribute> outCellAttr;
        outCells.reserve(static_cast<size_t>(mesh.getCellCount()) * 2);
        outCellAttr.reserve(outCells.capacity());

        auto emitTet = [&](Index a, Index b, Index c, Index d, Attribute attr)
        {
          std::array<Index, 4> t{{a, b, c, d}};
          if (signedVolume(t[0], t[1], t[2], t[3]) < 0)
            std::swap(t[1], t[2]);

          outCells.push_back(t);
          outCellAttr.push_back(attr);
        };

        // ---- Templates
        auto split_1neg = [&](Index vN,
                              Index p0, Index p1, Index p2,
                              Index i0, Index i1, Index i2,
                              Attribute aNeg, Attribute aPos)
        {
          emitTet(vN, i0, i1, i2, aNeg);

          emitTet(p0, p1, p2, i0, aPos);
          emitTet(p1, p2, i0, i1, aPos);
          emitTet(p2, i0, i1, i2, aPos);
        };

        auto split_1pos = [&](Index vP,
                              Index n0, Index n1, Index n2,
                              Index i0, Index i1, Index i2,
                              Attribute aNeg, Attribute aPos)
        {
          emitTet(vP, i0, i1, i2, aPos);

          emitTet(n0, n1, n2, i0, aNeg);
          emitTet(n1, n2, i0, i1, aNeg);
          emitTet(n2, i0, i1, i2, aNeg);
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
                                       Attribute aNeg, Attribute aPos)
        {
          // Candidate A (your original)
          std::array<std::array<Index, 4>, 6> A = {{
            {{n0, a, b, n1}},
            {{a,  b, n1, c}},
            {{b,  n1, c,  d}},
            {{p0, a, c,  p1}},
            {{a,  c, p1, b}},
            {{c,  p1, b, d}}
          }};

          // Candidate B (alternate)
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

          emitTet(best[0][0], best[0][1], best[0][2], best[0][3], aNeg);
          emitTet(best[1][0], best[1][1], best[1][2], best[1][3], aNeg);
          emitTet(best[2][0], best[2][1], best[2][2], best[2][3], aNeg);

          emitTet(best[3][0], best[3][1], best[3][2], best[3][3], aPos);
          emitTet(best[4][0], best[4][1], best[4][2], best[4][3], aPos);
          emitTet(best[5][0], best[5][1], best[5][2], best[5][3], aPos);
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
          const std::array<Sign, 4> s{{ sgn(phi[0]), sgn(phi[1]), sgn(phi[2]), sgn(phi[3]) }};

          const Attribute cellAttr = cit->getAttribute();

          int nneg = 0;
          for (int i = 0; i < 4; ++i)
            nneg += (s[i] == Sign::Negative);

          if (nneg == 0)
          {
            emitTet(v[0], v[1], v[2], v[3], splitLabel(3, cellAttr, /*negative*/ false));
            continue;
          }
          if (nneg == 4)
          {
            emitTet(v[0], v[1], v[2], v[3], splitLabel(3, cellAttr, /*negative*/ true));
            continue;
          }

          if (!shouldSplit(3, cellAttr))
          {
            emitTet(v[0], v[1], v[2], v[3], cellAttr);
            continue;
          }

          // Block by noSplit edges/faces
          bool blocked = false;

          const auto& eids = conn.getIncidence({3, 1}, ci);
          for (Index ei : eids)
          {
            const Attribute ea = mesh.getAttribute(1, ei);
            if (!isNoSplit(1, ea))
              continue;

            const auto& ev = conn.getIncidence({1, 0}, ei);
            const Sign sa = sgn(ls[ev[0]]);
            const Sign sb = sgn(ls[ev[1]]);

            const bool cut =
              (sa == Sign::Negative && sb == Sign::Positive) || (sa == Sign::Positive && sb == Sign::Negative) ||
              (sa == Sign::Zero && sb != Sign::Zero) || (sb == Sign::Zero && sa != Sign::Zero);

            if (cut)
            {
              blocked = true;
              break;
            }
          }

          if (!blocked)
          {
            const auto& fids = conn.getIncidence({3, 2}, ci);
            for (Index fi : fids)
            {
              const Attribute fa = mesh.getAttribute(2, fi);
              if (!isNoSplit(2, fa))
                continue;

              const auto& fv = conn.getIncidence({2, 0}, fi);
              const Sign s0 = sgn(ls[fv[0]]);
              const Sign s1 = sgn(ls[fv[1]]);
              const Sign s2 = sgn(ls[fv[2]]);

              const bool uniform = (s0 == s1 && s1 == s2);
              if (!uniform)
              {
                blocked = true;
                break;
              }
            }
          }

          if (blocked)
          {
            emitTet(v[0], v[1], v[2], v[3], cellAttr);
            continue;
          }

          const Attribute aNeg = splitLabel(3, cellAttr, /*negative*/ true);
          const Attribute aPos = splitLabel(3, cellAttr, /*negative*/ false);

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
            int e = edgeIdx(ia, ib);
            assert(e >= 0);
            const Index eid = eids[static_cast<size_t>(e)];
            const Index va = v[static_cast<size_t>(ia)];
            const Index vb = v[static_cast<size_t>(ib)];
            const Real fa = phi[static_cast<size_t>(ia)];
            const Real fb = phi[static_cast<size_t>(ib)];
            return getIntersectionOnEdge(eid, va, vb, fa, fb);
          };

          if (nneg == 1)
          {
            int in = -1; int ip[3]; int k = 0;
            for (int i = 0; i < 4; ++i)
            {
              if (s[i] == Sign::Negative) in = i;
              else ip[k++] = i;
            }

            const Index i0 = I(in, ip[0]);
            const Index i1 = I(in, ip[1]);
            const Index i2 = I(in, ip[2]);

            split_1neg(
              v[static_cast<size_t>(in)],
              v[static_cast<size_t>(ip[0])],
              v[static_cast<size_t>(ip[1])],
              v[static_cast<size_t>(ip[2])],
              i0, i1, i2,
              aNeg, aPos);
            continue;
          }

          if (nneg == 3)
          {
            int ipos = -1; int ineg[3]; int k = 0;
            for (int i = 0; i < 4; ++i)
            {
              if (s[i] == Sign::Positive) ipos = i;
              else ineg[k++] = i;
            }

            const Index i0 = I(ipos, ineg[0]);
            const Index i1 = I(ipos, ineg[1]);
            const Index i2 = I(ipos, ineg[2]);

            split_1pos(
              v[static_cast<size_t>(ipos)],
              v[static_cast<size_t>(ineg[0])],
              v[static_cast<size_t>(ineg[1])],
              v[static_cast<size_t>(ineg[2])],
              i0, i1, i2,
              aNeg, aPos);
            continue;
          }

          // nneg == 2
          {
            int ineg[2]; int ipos[2]; int kn = 0, kp = 0;
            for (int i = 0; i < 4; ++i)
            {
              if (s[i] == Sign::Negative) ineg[kn++] = i;
              else ipos[kp++] = i;
            }

            const Index a = I(ineg[0], ipos[0]); // t_n0p0
            const Index b = I(ineg[0], ipos[1]); // t_n0p1
            const Index c = I(ineg[1], ipos[0]); // t_n1p0
            const Index d = I(ineg[1], ipos[1]); // t_n1p1

            split_2neg2pos_best(
              v[static_cast<size_t>(ineg[0])],
              v[static_cast<size_t>(ineg[1])],
              v[static_cast<size_t>(ipos[0])],
              v[static_cast<size_t>(ipos[1])],
              a, b, c, d,
              aNeg, aPos);
            continue;
          }
        }

        // ---- Build output mesh
        typename MeshType::Builder builder;
        builder.initialize(3).nodes(outVerts.getCount());
        builder.setVertices(std::move(outVerts));
        builder.reserve(3, outCells.size());
        for (size_t i = 0; i < outCells.size(); ++i)
        {
          const auto& t = outCells[i];
          builder.tetrahedron(IndexArray{{ t[0], t[1], t[2], t[3] }});
          builder.attribute({3, static_cast<Index>(i)}, outCellAttr[i]);
        }
        return builder.finalize();
      }
  };

  template <class Mesh, class Data>
  MarchingTetrahedra(const Variational::GridFunction<Variational::P1<Real, Mesh>, Data>& ls)
    -> MarchingTetrahedra<Variational::GridFunction<Variational::P1<Real, Mesh>, Data>>;
}

#endif
