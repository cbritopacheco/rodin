/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_MARCHINGTRIANGLES_H
#define RODIN_GEOMETRY_MARCHINGTRIANGLES_H

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
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
  class MarchingTriangles;

  template <class Mesh, class Data>
  class MarchingTriangles<Variational::GridFunction<Variational::P1<Real, Mesh>, Data>>
    : public MarchingBase<Variational::GridFunction<Variational::P1<Real, Mesh>, Data>>
  {
    public:
      using Parent = MarchingBase<Variational::GridFunction<Variational::P1<Real, Mesh>, Data>>;
      using MeshType = typename Parent::MeshType;
      using GridFunctionType = typename Parent::GridFunctionType;
      using SplitMap = typename Parent::SplitMap;
      using NoSplitT = typename Parent::NoSplitT;
      using Split = typename Parent::Split;

      MarchingTriangles(const GridFunctionType& ls)
        : Parent(ls),
          // Defaults mirror MarchingTetrahedra and can be tuned by callers.
          // - sign tolerance: classify near-zero level-set values.
          // - snap tolerance: snap edge intersections close to endpoints.
          m_sign_tolerance(1e-12),
          m_snap_tolerance(1e-12)
      {}

      MeshType discretize() const override
      {
        const auto& ls = this->getGridFunction();
        MeshType mesh(ls.getFiniteElementSpace().getMesh());

        if (mesh.getDimension() != 2)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "MarchingTriangles expects a 2D volume mesh."
            << Alert::Raise;
        }

        auto& conn = mesh.getConnectivity();
        conn.compute(2, 1);
        conn.compute(1, 0);
        const Index nv = static_cast<Index>(mesh.getVertexCount());
        const Index ne = static_cast<Index>(mesh.getPolytopeCount(1));
        static constexpr Index InvalidIndex = std::numeric_limits<Index>::max();

        const Real eps_sign = this->getSignTolerance();
        const Real eps_snap = this->getSnapTolerance();

        auto signOf = [&](Real x) -> int
        {
          if (x < -eps_sign) return -1;
          if (x >  eps_sign) return 1;
          return 0;
        };

        auto isCut = [](int sa, int sb)
        {
          if (sa == 0 && sb == 0) return false;
          return sa != sb;
        };

        static constexpr int8_t SideUnknown  = -1;
        static constexpr int8_t SidePositive = 0;
        static constexpr int8_t SideNegative = 1;

        auto isNoSplit = [&](size_t d, const Optional<Attribute>& a) -> bool
        {
          if (!a) return false;
          const auto& sm = this->getSplitMap(d);
          if (sm.empty()) return false;
          const auto it = sm.find(*a);
          return it != sm.end() && std::holds_alternative<NoSplitT>(it->second);
        };

        auto shouldSplit = [&](size_t d, const Optional<Attribute>& a) -> bool
        {
          const auto& sm = this->getSplitMap(d);
          if (sm.empty()) return true;
          if (!a) return true;
          const auto it = sm.find(*a);
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
            const auto it = sm.find(*base);
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

        struct EdgeKey
        {
          Index a;
          Index b;
          bool operator<(const EdgeKey& o) const { return std::tie(a, b) < std::tie(o.a, o.b); }
        };

        auto makeEdgeKey = [](Index a, Index b) -> EdgeKey
        {
          if (a > b) std::swap(a, b);
          return { a, b };
        };

        std::map<EdgeKey, Attribute> inEdgeAttrByVerts;
        for (auto eit = mesh.getPolytope(1); !eit.end(); ++eit)
        {
          const auto a = eit->getAttribute();
          if (!a) continue;
          const auto& ev = eit->getVertices();
          inEdgeAttrByVerts.emplace(makeEdgeKey(ev(0), ev(1)), *a);
        }

        PointCloud outVerts;
        outVerts.setDimension(mesh.getSpaceDimension());
        // Typical cut patterns create O(number_of_edges) intersections;
        // reserving half the edge count keeps growth amortized without over-allocation.
        outVerts.reserve(static_cast<size_t>(nv) + static_cast<size_t>(ne / 2));
        for (Index i = 0; i < nv; ++i)
          outVerts.push_back(mesh.getVertexCoordinates(i));

        auto addVertex = [&](const Math::SpatialPoint& p) -> Index
        {
          outVerts.push_back(p);
          return static_cast<Index>(outVerts.getCount() - 1);
        };

        std::vector<Index> isect(static_cast<size_t>(ne), InvalidIndex);
        auto getIntersection = [&](Index edgeIdx, Index va, Index vb, Real fa, Real fb) -> Index
        {
          if (std::abs(fa) <= eps_snap) return va;
          if (std::abs(fb) <= eps_snap) return vb;

          Index& cached = isect[static_cast<size_t>(edgeIdx)];
          if (cached != InvalidIndex) return cached;

          const Real denom = (fa - fb);
          if (std::abs(denom) <= eps_sign)
            return (std::abs(fa) < std::abs(fb)) ? va : vb;

          Real t = fa / denom;
          t = std::min<Real>(Real(1), std::max<Real>(Real(0), t));
          if (t <= eps_snap) return va;
          if (t >= Real(1) - eps_snap) return vb;

          const auto p = outVerts[va] + t * (outVerts[vb] - outVerts[va]);
          cached = addVertex(p);
          return cached;
        };

        std::vector<std::array<Index, 3>> outCells;
        std::vector<Optional<Attribute>> outCellAttr;
        std::vector<int8_t> outCellSide;
        outCells.reserve(static_cast<size_t>(mesh.getCellCount()) * 3);
        outCellAttr.reserve(outCells.capacity());
        outCellSide.reserve(outCells.capacity());

        auto emitTri = [&](Index a, Index b, Index c, const Optional<Attribute>& attr, int8_t side)
        {
          outCells.push_back({{a, b, c}});
          outCellAttr.push_back(attr);
          outCellSide.push_back(side);
        };

        for (auto cit = mesh.getCell(); !cit.end(); ++cit)
        {
          if (cit->getGeometry() != Polytope::Type::Triangle)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "MarchingTriangles only supports triangular meshes."
              << Alert::Raise;
          }

          const Index ci = cit->getIndex();
          const auto& cv = cit->getVertices();
          const std::array<Index, 3> v{{ cv(0), cv(1), cv(2) }};
          const std::array<Real, 3> phi{{ ls[v[0]], ls[v[1]], ls[v[2]] }};
          int ss[3]{ signOf(phi[0]), signOf(phi[1]), signOf(phi[2]) };

          int nneg = 0;
          for (int i = 0; i < 3; ++i) nneg += (ss[i] < 0);

          const auto cellAttr = cit->getAttribute();
          if (nneg == 0 || nneg == 3)
          {
            const int8_t side = (nneg == 3) ? SideNegative : SidePositive;
            emitTri(v[0], v[1], v[2], inferLabel(2, cellAttr, side == SideNegative), side);
            continue;
          }

          if (!shouldSplit(2, cellAttr))
          {
            emitTri(v[0], v[1], v[2], cellAttr, SideUnknown);
            continue;
          }

          bool blocked = false;
          const auto& eids = conn.getIncidence({2, 1}, ci);
          for (const auto ei : eids)
          {
            if (!isNoSplit(1, mesh.getAttribute(1, ei)))
              continue;
            const auto& ev = conn.getIncidence({1, 0}, ei);
            if (isCut(signOf(ls[ev[0]]), signOf(ls[ev[1]])))
            {
              blocked = true;
              break;
            }
          }

          if (blocked)
          {
            emitTri(v[0], v[1], v[2], cellAttr, SideUnknown);
            continue;
          }

          std::map<EdgeKey, Index> cellEdgeIdByVerts;
          for (const auto ei : eids)
          {
            const auto& ev = conn.getIncidence({1, 0}, ei);
            cellEdgeIdByVerts.emplace(makeEdgeKey(ev[0], ev[1]), ei);
          }
          auto I = [&](int a, int b) -> Index
          {
            const auto va = v[static_cast<size_t>(a)];
            const auto vb = v[static_cast<size_t>(b)];
            const auto it = cellEdgeIdByVerts.find(makeEdgeKey(va, vb));
            assert(it != cellEdgeIdByVerts.end());
            return getIntersection(
              it->second,
              va, vb,
              phi[static_cast<size_t>(a)],
              phi[static_cast<size_t>(b)]);
          };

          const auto aNeg = inferLabel(2, cellAttr, true);
          const auto aPos = inferLabel(2, cellAttr, false);

          if (nneg == 1)
          {
            int in = -1;
            int ip[2];
            int k = 0;
            for (int i = 0; i < 3; ++i)
              (ss[i] < 0) ? in = i : ip[k++] = i;

            const Index i0 = I(in, ip[0]);
            const Index i1 = I(in, ip[1]);

            emitTri(v[static_cast<size_t>(in)], i0, i1, aNeg, SideNegative);
            emitTri(v[static_cast<size_t>(ip[0])], v[static_cast<size_t>(ip[1])], i0, aPos, SidePositive);
            emitTri(v[static_cast<size_t>(ip[1])], i0, i1, aPos, SidePositive);
          }
          else // nneg == 2
          {
            int ip = -1;
            int in[2];
            int k = 0;
            for (int i = 0; i < 3; ++i)
              (ss[i] < 0) ? in[k++] = i : ip = i;

            const Index i0 = I(ip, in[0]);
            const Index i1 = I(ip, in[1]);

            emitTri(v[static_cast<size_t>(ip)], i0, i1, aPos, SidePositive);
            emitTri(v[static_cast<size_t>(in[0])], v[static_cast<size_t>(in[1])], i0, aNeg, SideNegative);
            emitTri(v[static_cast<size_t>(in[1])], i0, i1, aNeg, SideNegative);
          }
        }

        typename MeshType::Builder builder;
        builder.initialize(2).nodes(outVerts.getCount());
        builder.setVertices(std::move(outVerts));
        builder.reserve(2, outCells.size());
        for (size_t i = 0; i < outCells.size(); ++i)
        {
          builder.triangle(IndexArray{{ outCells[i][0], outCells[i][1], outCells[i][2] }});
          if (outCellAttr[i])
            builder.attribute({2, static_cast<Index>(i)}, *outCellAttr[i]);
        }

        auto out = builder.finalize();
        out.getConnectivity().compute(1, 2);
        const auto& oconn = out.getConnectivity();
        const auto& ifaceOpt = this->getInterface(2);

        auto isFittedZeroV = [&](Index vi) -> bool
        {
          return (vi < nv) && (std::abs(ls[vi]) <= eps_sign);
        };

        for (auto eit = out.getPolytope(1); !eit.end(); ++eit)
        {
          const Index eidx = eit->getIndex();
          const auto& ev = eit->getVertices();
          const bool allOriginal = (ev(0) < nv) && (ev(1) < nv);
          const auto& inc = oconn.getIncidence({1, 2}, eidx);
          bool hasNeg = false;
          bool hasPos = false;
          for (const auto ci : inc)
          {
            hasNeg = hasNeg || (outCellSide[ci] == SideNegative);
            hasPos = hasPos || (outCellSide[ci] == SidePositive);
          }
          const bool fittedInterface = allOriginal && isFittedZeroV(ev(0)) && isFittedZeroV(ev(1));
          if ((hasNeg && hasPos) || fittedInterface)
          {
            if (this->getInterface(1))
              out.setAttribute({1, eidx}, *this->getInterface(1));
            continue;
          }

          if (!allOriginal)
            continue;
          const auto in = inEdgeAttrByVerts.find(makeEdgeKey(ev(0), ev(1)));
          if (in == inEdgeAttrByVerts.end())
            continue;
          if (ifaceOpt && in->second == *ifaceOpt)
            continue;
          const auto mapped = inferLabel(1, Optional<Attribute>(in->second), hasNeg);
          if (mapped)
            out.setAttribute({1, eidx}, *mapped);
        }

        return out;
      }

      MarchingTriangles& setSignTolerance(Real tol)
      {
        m_sign_tolerance = std::max<Real>(tol, 0);
        return *this;
      }

      Real getSignTolerance() const
      {
        return m_sign_tolerance;
      }

      MarchingTriangles& setSnapTolerance(Real tol)
      {
        m_snap_tolerance = std::max<Real>(tol, 0);
        return *this;
      }

      Real getSnapTolerance() const
      {
        return m_snap_tolerance;
      }

    private:
      Real m_sign_tolerance;
      Real m_snap_tolerance;
  };

  template <class Mesh, class Data>
  MarchingTriangles(const Variational::GridFunction<Variational::P1<Real, Mesh>, Data>& ls)
    -> MarchingTriangles<Variational::GridFunction<Variational::P1<Real, Mesh>, Data>>;
}

#endif
