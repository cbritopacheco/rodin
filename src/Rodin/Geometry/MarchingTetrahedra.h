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

        const Real eps_sign_user = this->getSignTolerance();
        const Real eps_snap_user = this->getSnapTolerance();

        auto finite = [](Real x) -> bool { return std::isfinite(x); };

        auto isCut = [&](int sa, int sb) -> bool
        {
          if (sa == 0 && sb == 0) return false;
          return sa != sb;
        };

        static constexpr int8_t SideUnknown  = -1;
        static constexpr int8_t SidePositive = 0;
        static constexpr int8_t SideNegative = 1;

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

        // SplitMap helpers
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

        // ---- provenance maps for ORIGINAL edges/faces only
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

        // ---- output vertices + edge intersection cache
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

        std::vector<Index> isect(static_cast<size_t>(ne), InvalidIndex);

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

          static constexpr Real tminFrac = Real(1e-2);
          Real t = fa / denom;

          if (std::abs(fa + fb) <= Real(4) * eps0)
            t = Real(0.5);
          else
            t = std::min<Real>(Real(1) - tminFrac, std::max<Real>(tminFrac, t));

          const auto pa = outVerts[va];
          const auto pb = outVerts[vb];
          const auto p  = pa + t * (pb - pa);

          cached = addVertex(p);
          return cached;
        };

        // ---- output cells + metadata
        std::vector<std::array<Index, 4>> outCells;
        std::vector<Optional<Attribute>>  outCellAttr;
        std::vector<int8_t>               outCellSide;

        outCells.reserve(static_cast<size_t>(mesh.getCellCount()) * 2);
        outCellAttr.reserve(outCells.capacity());
        outCellSide.reserve(outCells.capacity());

        std::vector<Index> edgeCellDegree(static_cast<size_t>(mesh.getPolytopeCount(1)), 0);
        for (auto cit = mesh.getCell(); !cit.end(); ++cit)
        {
          for (const auto ei : conn.getIncidence({3, 1}, cit->getIndex()))
            edgeCellDegree[static_cast<size_t>(ei)]++;
        }

        auto signedVolume = [&](Index a, Index b, Index c, Index d) -> Real
        {
          const auto A = outVerts[a];
          const auto B = outVerts[b];
          const auto C = outVerts[c];
          const auto D = outVerts[d];
          const auto u = B - A;
          const auto v = C - A;
          const auto w = D - A;
          return u.dot(v.cross(w));
        };
        auto tetQuality = [&](Index a, Index b, Index c, Index d) -> Real
        {
          // Lower bound for denominator in the scale-normalized quality metric.
          // Kept very small to avoid perturbing quality ordering while preventing
          // division by zero in highly collapsed configurations.
          // Safety factor chosen to keep denominator strictly positive while
          // remaining close to machine precision for normalized coordinates.
          static constexpr Real epsilonSafetyFactor = Real(100);
          const Real minQualityDivisor =
            epsilonSafetyFactor * std::numeric_limits<Real>::epsilon();
          const Real absoluteVolume = std::abs(signedVolume(a, b, c, d));
          const auto A = outVerts[a];
          const auto B = outVerts[b];
          const auto C = outVerts[c];
          const auto D = outVerts[d];
          const Real l01 = (A - B).norm();
          const Real l02 = (A - C).norm();
          const Real l03 = (A - D).norm();
          const Real l12 = (B - C).norm();
          const Real l13 = (B - D).norm();
          const Real l23 = (C - D).norm();
          Real h = std::max(l01, l02);
          h = std::max(h, l03);
          h = std::max(h, l12);
          h = std::max(h, l13);
          h = std::max(h, l23);
          h = std::max(h, minQualityDivisor);
          return absoluteVolume / (h * h * h);
        };
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
        auto minQuality4 = [&](const std::array<std::array<Index, 4>, 4>& tets) -> Real
        {
          Real minQuality = std::numeric_limits<Real>::infinity();
          for (const auto& tet : tets)
            minQuality = std::min(minQuality, tetQuality(tet[0], tet[1], tet[2], tet[3]));
          return minQuality;
        };
        auto bestAndQualityOf3ByMinQuality4 = [&](const std::array<std::array<Index, 4>, 4>& A,
                                                  const std::array<std::array<Index, 4>, 4>& B,
                                                  const std::array<std::array<Index, 4>, 4>& C)
          -> std::pair<const std::array<std::array<Index, 4>, 4>*, Real>
        {
          const Real qa = minQuality4(A);
          const Real qb = minQuality4(B);
          const Real qc = minQuality4(C);
          const std::array<std::array<Index, 4>, 4>* best = &A;
          Real qbest = qa;
          if (qb > qbest) { best = &B; qbest = qb; }
          if (qc > qbest) { best = &C; qbest = qc; }
          return {best, qbest};
        };
        auto split_1neg = [&](Index vN,
                              Index p0, Index p1, Index p2,
                              Index i0, Index i1, Index i2,
                              const Optional<Attribute>& aNeg,
                              const Optional<Attribute>& aPos,
                              bool allowLocalNoSplit)
          -> bool
        {
          std::array<std::array<Index, 4>, 4> A = {{
            {{vN, i0, i1, i2}},
            {{p0, p1, p2, i0}},
            {{p1, p2, i0, i1}},
            {{p2, i0, i1, i2}}
          }};
          std::array<std::array<Index, 4>, 4> B = {{
            {{vN, i0, i1, i2}},
            {{p0, p1, p2, i1}},
            {{p0, p2, i1, i2}},
            {{p0, i0, i1, i2}}
          }};
          std::array<std::array<Index, 4>, 4> C = {{
            {{vN, i0, i1, i2}},
            {{p0, p1, p2, i2}},
            {{p0, p1, i0, i2}},
            {{p1, i0, i1, i2}}
          }};
          const auto [best, bestQuality] = bestAndQualityOf3ByMinQuality4(A, B, C);
          if (allowLocalNoSplit && !finite(bestQuality))
            return false;
          emitTet((*best)[0][0], (*best)[0][1], (*best)[0][2], (*best)[0][3], aNeg, SideNegative);
          emitTet((*best)[1][0], (*best)[1][1], (*best)[1][2], (*best)[1][3], aPos, SidePositive);
          emitTet((*best)[2][0], (*best)[2][1], (*best)[2][2], (*best)[2][3], aPos, SidePositive);
          emitTet((*best)[3][0], (*best)[3][1], (*best)[3][2], (*best)[3][3], aPos, SidePositive);
          return true;
        };

        auto split_1pos = [&](Index vP,
                              Index n0, Index n1, Index n2,
                              Index i0, Index i1, Index i2,
                              const Optional<Attribute>& aNeg,
                              const Optional<Attribute>& aPos,
                              bool allowLocalNoSplit)
          -> bool
        {
          std::array<std::array<Index, 4>, 4> A = {{
            {{vP, i0, i1, i2}},
            {{n0, n1, n2, i0}},
            {{n1, n2, i0, i1}},
            {{n2, i0, i1, i2}}
          }};
          std::array<std::array<Index, 4>, 4> B = {{
            {{vP, i0, i1, i2}},
            {{n0, n1, n2, i1}},
            {{n0, n2, i1, i2}},
            {{n0, i0, i1, i2}}
          }};
          std::array<std::array<Index, 4>, 4> C = {{
            {{vP, i0, i1, i2}},
            {{n0, n1, n2, i2}},
            {{n0, n1, i0, i2}},
            {{n1, i0, i1, i2}}
          }};
          const auto [best, bestQuality] = bestAndQualityOf3ByMinQuality4(A, B, C);
          if (allowLocalNoSplit && !finite(bestQuality))
            return false;
          emitTet((*best)[0][0], (*best)[0][1], (*best)[0][2], (*best)[0][3], aPos, SidePositive);
          emitTet((*best)[1][0], (*best)[1][1], (*best)[1][2], (*best)[1][3], aNeg, SideNegative);
          emitTet((*best)[2][0], (*best)[2][1], (*best)[2][2], (*best)[2][3], aNeg, SideNegative);
          emitTet((*best)[3][0], (*best)[3][1], (*best)[3][2], (*best)[3][3], aNeg, SideNegative);
          return true;
        };

        auto minQuality6 = [&](const std::array<std::array<Index, 4>, 6>& tets) -> Real
        {
          Real minQuality = std::numeric_limits<Real>::infinity();
          for (const auto& tet : tets)
            minQuality = std::min(minQuality, tetQuality(tet[0], tet[1], tet[2], tet[3]));
          return minQuality;
        };
        auto split_2neg2pos_best = [&](Index n0, Index n1,
                                       Index p0, Index p1,
                                       Index a, Index b, Index c, Index d,
                                       const Optional<Attribute>& aNeg,
                                       const Optional<Attribute>& aPos,
                                       bool allowLocalNoSplit)
          -> bool
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

          const Real qa = minQuality6(A);
          const Real qb = minQuality6(B);
          const auto& best = (qb > qa) ? B : A;
          if (allowLocalNoSplit && !finite(std::max(qa, qb)))
            return false;

          emitTet(best[0][0], best[0][1], best[0][2], best[0][3], aNeg, SideNegative);
          emitTet(best[1][0], best[1][1], best[1][2], best[1][3], aNeg, SideNegative);
          emitTet(best[2][0], best[2][1], best[2][2], best[2][3], aNeg, SideNegative);

          emitTet(best[3][0], best[3][1], best[3][2], best[3][3], aPos, SidePositive);
          emitTet(best[4][0], best[4][1], best[4][2], best[4][3], aPos, SidePositive);
          emitTet(best[5][0], best[5][1], best[5][2], best[5][3], aPos, SidePositive);
          return true;
        };

        static constexpr std::array<std::pair<int,int>, 6> TetEdges = {{
          {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
        }};

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

        // ---- main loop
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

          if (!finite(phi[0]) || !finite(phi[1]) || !finite(phi[2]) || !finite(phi[3]))
          {
            emitTet(v[0], v[1], v[2], v[3], cellAttr, SideUnknown);
            continue;
          }

          const Real phi_scale = std::max(std::max(std::abs(phi[0]), std::abs(phi[1])),
                                std::max(std::abs(phi[2]), std::abs(phi[3])));
          const Real rel0 = Real(1e-10);
          const Real eps0 = std::max(eps_sign_user, rel0 * phi_scale);
          const Real eps_snap = std::max(eps_snap_user, eps0);

          int ss[4] = {
            (phi[0] < -eps0) ? -1 : (phi[0] > eps0 ? +1 : 0),
            (phi[1] < -eps0) ? -1 : (phi[1] > eps0 ? +1 : 0),
            (phi[2] < -eps0) ? -1 : (phi[2] > eps0 ? +1 : 0),
            (phi[3] < -eps0) ? -1 : (phi[3] > eps0 ? +1 : 0)
          };

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

          if (!shouldSplit(3, cellAttr))
          {
            emitTet(v[0], v[1], v[2], v[3], cellAttr, SideUnknown);
            continue;
          }

          bool blocked = false;
          const auto& eids = conn.getIncidence({3, 1}, ci);

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

          const auto aNeg = inferLabel(3, cellAttr, true);
          const auto aPos = inferLabel(3, cellAttr, false);

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
            return getIntersectionOnEdge(eid, va, vb, fa, fb, eps0, eps_snap);
          };

          int nneg = 0;
          for (int i = 0; i < 4; ++i) nneg += (ss[i] < 0);
          bool allowLocalNoSplit = true;
          for (size_t ei = 0; ei < TetEdges.size(); ++ei)
          {
            const auto [ea, eb] = TetEdges[ei];
            if (ss[ea] * ss[eb] != -1)
              continue;
            if (edgeCellDegree[static_cast<size_t>(eids[ei])] > 1)
            {
              allowLocalNoSplit = false;
              break;
            }
          }

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
            if (!split_1neg(v[static_cast<size_t>(in)],
                            v[static_cast<size_t>(ip[0])],
                            v[static_cast<size_t>(ip[1])],
                            v[static_cast<size_t>(ip[2])],
                            i0, i1, i2, aNeg, aPos, allowLocalNoSplit))
            {
              emitTet(v[0], v[1], v[2], v[3], cellAttr, SideUnknown);
            }
            continue;
          }

          if (nneg == 3)
          {
            int ipos = -1; int ineg[3]; int k = 0;
            for (int i = 0; i < 4; ++i) (ss[i] > 0) ? (ipos = i) : (ineg[k++] = i);

            const Index i0 = I(ipos, ineg[0]);
            const Index i1 = I(ipos, ineg[1]);
            const Index i2 = I(ipos, ineg[2]);
            if (!split_1pos(v[static_cast<size_t>(ipos)],
                            v[static_cast<size_t>(ineg[0])],
                            v[static_cast<size_t>(ineg[1])],
                            v[static_cast<size_t>(ineg[2])],
                            i0, i1, i2, aNeg, aPos, allowLocalNoSplit))
            {
              emitTet(v[0], v[1], v[2], v[3], cellAttr, SideUnknown);
            }
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
            if (!split_2neg2pos_best(v[static_cast<size_t>(ineg[0])],
                                     v[static_cast<size_t>(ineg[1])],
                                     v[static_cast<size_t>(ipos[0])],
                                     v[static_cast<size_t>(ipos[1])],
                                     a, b, c, d, aNeg, aPos, allowLocalNoSplit))
            {
              emitTet(v[0], v[1], v[2], v[3], cellAttr, SideUnknown);
            }
            continue;
          }
        }

        // ---- build output mesh (cells)
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
        // Attribute transfer / interface marking (+ m_old, m_fallback policy)
        //
        // Cases handled:
        //  (A) Interface face/edge: mark with ifaceOpt (if configured).
        //  (B) Original (allOriginal) non-interface:
        //      - if input base == ifaceOpt:
        //          * if m_old set -> set to *m_old
        //          * else -> leave nullopt (legacy “drop old interface tag”)
        //      - else: map with inferLabel(...) and set if not nullopt.
        //  (C) Created (NOT allOriginal) non-interface:
        //      - if m_fallback set -> set to *m_fallback
        //      - else -> leave nullopt (legacy behavior)
        // ------------------------------------------------------------
        out.getConnectivity().compute(2, 3);
        out.getConnectivity().compute(1, 3);
        const auto& oconn = out.getConnectivity();

        const auto& ifaceOpt = this->getInterface(2);

        auto isFittedZeroV = [&](Index vi) -> bool
        {
          return (vi < nv) && (std::abs(ls[vi]) <= eps_sign_user);
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
            hasNeg = hasNeg || (outCellSide[ci] == SideNegative);
            hasPos = hasPos || (outCellSide[ci] == SidePositive);
          }

          const bool fittedInterface =
            allOriginal && isFittedZeroV(fv(0)) && isFittedZeroV(fv(1)) && isFittedZeroV(fv(2));

          const bool isInterface = (hasNeg && hasPos) || fittedInterface;
          if (isInterface)
          {
            if (ifaceOpt) out.setAttribute({2, fidx}, *ifaceOpt);
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

          if (ifaceOpt && base == *ifaceOpt)
          {
            if (m_old[2])
              out.setAttribute({2, fidx}, *m_old[2]);
            continue;
          }

          const Optional<Attribute> mapped = inferLabel(2, Optional<Attribute>(base), /*neg*/ hasNeg);
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
            hasNeg = hasNeg || (outCellSide[ci] == SideNegative);
            hasPos = hasPos || (outCellSide[ci] == SidePositive);
          }

          const bool fittedInterface =
            allOriginal && isFittedZeroV(ev(0)) && isFittedZeroV(ev(1));

          const bool isInterface = (hasNeg && hasPos) || fittedInterface;
          if (isInterface)
          {
            if (this->getInterface(1))
              out.setAttribute({1, eidx}, *this->getInterface(1));
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

          if (ifaceOpt && base == *ifaceOpt)
          {
            if (m_old[1])
              out.setAttribute({1, eidx}, *m_old[1]);
            continue;
          }

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
