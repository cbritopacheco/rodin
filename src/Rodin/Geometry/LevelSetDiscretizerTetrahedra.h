/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_LEVELSETDISCRETIZERTETRAHEDRA_H
#define RODIN_GEOMETRY_LEVELSETDISCRETIZERTETRAHEDRA_H

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <variant>
#include <vector>

#include "Rodin/Alert/MemberFunctionException.h"

#include "LevelSetDiscretizerBase.h"
#include "Mesh.h"
#include "Rodin/Geometry/Polytope.h"

namespace Rodin::Geometry
{
  template <class... Params>
  class LevelSetDiscretizerTetrahedra;

  template <class Mesh, class Data>
  class LevelSetDiscretizerTetrahedra<
    Variational::GridFunction<Variational::P1<Real, Mesh>, Data>>
    : public LevelSetDiscretizerBase<
        Variational::GridFunction<Variational::P1<Real, Mesh>, Data>>
  {
    public:
      using Parent = LevelSetDiscretizerBase<
        Variational::GridFunction<Variational::P1<Real, Mesh>, Data>>;

      using FESType          = typename Parent::FESType;
      using MeshType         = typename Parent::MeshType;
      using GridFunctionType = typename Parent::GridFunctionType;

      using Split            = typename Parent::Split;
      using PreserveAttrT    = typename Parent::PreserveAttributeT;
      using SplitMap         = typename Parent::SplitMap;

      using Point = Math::SpatialPoint;
      using Tet   = std::array<Index, 4>;
      using Tri   = std::array<Index, 3>;

      using QualityMetric = std::function<long double(const PointCloud&, const Tet&)>;

      LevelSetDiscretizerTetrahedra(const GridFunctionType& ls,
                                    QualityMetric q = QualityMetric{})
        : Parent(ls),
          m_sign_tolerance(1e-12),
          m_snap_tolerance(1e-12),
          m_min_quality(1e-10L),
          m_quality_metric(q ? std::move(q) : defaultQualityMetric())
      {}

      LevelSetDiscretizerTetrahedra& setSignTolerance(Real tol)
      {
        m_sign_tolerance = tol;
        return *this;
      }

      LevelSetDiscretizerTetrahedra& setSnapTolerance(Real tol)
      {
        m_snap_tolerance = tol;
        return *this;
      }

      LevelSetDiscretizerTetrahedra& setQualityMetric(QualityMetric q)
      {
        m_quality_metric = q ? std::move(q) : defaultQualityMetric();
        return *this;
      }

      LevelSetDiscretizerTetrahedra& setMinimumQuality(long double qmin)
      {
        m_min_quality = std::max(0.0L, qmin);
        return *this;
      }

      Real getSignTolerance() const
      {
        return m_sign_tolerance;
      }

      Real getSnapTolerance() const
      {
        return m_snap_tolerance;
      }

      const QualityMetric& getQualityMetric() const
      {
        return m_quality_metric;
      }

      long double getMinimumQuality() const
      {
        return m_min_quality;
      }

      LevelSetDiscretizerTetrahedra& setOld(size_t d, const Optional<Attribute>& a)
      {
        m_old[d] = a;
        return *this;
      }

      LevelSetDiscretizerTetrahedra& setFallback(size_t d, const Optional<Attribute>& a)
      {
        m_fallback[d] = a;
        return *this;
      }

      MeshType discretize() const override
      {
        const auto& ls   = this->getGridFunction();
        const auto& mesh = ls.getFiniteElementSpace().getMesh();

        if (mesh.getDimension() != 3)
          Alert::MemberFunctionException(*this, __func__)
            << "Expected 3D mesh."
            << Alert::Raise;

        auto& conn = mesh.getConnectivity();

        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 3, 2);
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 3, 1);
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 2, 1);
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 1, 0);

        static constexpr Index INVALID = std::numeric_limits<Index>::max();

        const Index nv = static_cast<Index>(mesh.getVertexCount());
        const Index ne = static_cast<Index>(conn.getCount(1));
        const Index nf = static_cast<Index>(conn.getCount(2));

        const Real eps_sign_user = this->getSignTolerance();
        const Real eps_snap_user = this->getSnapTolerance();

        auto finite = [](Real x) { return std::isfinite(x); };

        // --------------------------------------------------------------------
        // Attribute transfer helpers
        // --------------------------------------------------------------------
        auto inferLabel = [&](size_t d, const Optional<Attribute>& base, bool negativeSide)
          -> Optional<Attribute>
        {
          if (!base)
            return {};

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
                  if constexpr (std::is_same_v<T, PreserveAttrT>)
                  {
                    return base;
                  }
                  else
                  {
                    return negativeSide ? Optional<Attribute>(v.negative)
                                        : Optional<Attribute>(v.positive);
                  }
                },
                it->second);
            }
          }

          return base;
        };

        // --------------------------------------------------------------------
        // Provenance maps for original lower-dimensional attributes
        // --------------------------------------------------------------------
        struct EdgeKey
        {
          Index a, b;
          bool operator<(const EdgeKey& o) const
          {
            return std::tie(a, b) < std::tie(o.a, o.b);
          }
        };

        struct FaceKey
        {
          Index a, b, c;
          bool operator<(const FaceKey& o) const
          {
            return std::tie(a, b, c) < std::tie(o.a, o.b, o.c);
          }
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
        auto toOptionalAttr = [](const auto& attr) -> Optional<Attribute>
        {
          if (attr)
            return Optional<Attribute>(*attr);
          return {};
        };

        for (auto eit = mesh.getPolytope(1); !eit.end(); ++eit)
        {
          const auto& ev = eit->getVertices();
          const Optional<Attribute> a = toOptionalAttr(eit->getAttribute());
          if (a)
            inEdgeAttrByVerts[makeEdgeKey(ev(0), ev(1))] = *a;
        }

        for (auto fit = mesh.getPolytope(2); !fit.end(); ++fit)
        {
          if (fit->getGeometry() != Polytope::Type::Triangle)
            continue;
          const auto& fv = fit->getVertices();
          const Optional<Attribute> a = toOptionalAttr(fit->getAttribute());
          if (a)
            inFaceAttrByVerts[makeFaceKey(fv(0), fv(1), fv(2))] = *a;
        }

        // --------------------------------------------------------------------
        // Global sign scale
        // --------------------------------------------------------------------
        static constexpr Real rel0_global = Real(1e-10);
        Real globalPhiScale = Real(0);
        for (Index i = 0; i < nv; ++i)
          if (finite(ls[i]))
            globalPhiScale = std::max(globalPhiScale, std::abs(ls[i]));

        const Real eps_sign = std::max(eps_sign_user, rel0_global * globalPhiScale);
        const Real eps_snap = std::max(eps_snap_user, eps_sign);

        // --------------------------------------------------------------------
        // Global sign classification
        // --------------------------------------------------------------------
        enum class Sign : int8_t
        {
          Negative = -1,
          Zero     =  0,
          Positive =  1,
          Invalid  =  2
        };

        auto classify = [&](Real x) -> Sign
        {
          if (!finite(x))          return Sign::Invalid;
          if (x < -eps_sign)       return Sign::Negative;
          if (x >  eps_sign)       return Sign::Positive;
          return Sign::Zero;
        };

        std::vector<Sign> vsgn(static_cast<size_t>(nv), Sign::Invalid);
        for (Index i = 0; i < nv; ++i)
          vsgn[(size_t)i] = classify(ls[i]);

        auto signToInt = [](Sign s) -> int
        {
          switch (s)
          {
            case Sign::Negative: return -1;
            case Sign::Zero:     return  0;
            case Sign::Positive: return  1;
            default:             return  2;
          }
        };

        // --------------------------------------------------------------------
        // Output vertices: originals first
        // --------------------------------------------------------------------
        PointCloud outVerts;
        outVerts.setDimension(3);
        outVerts.reserve(static_cast<size_t>(nv) + static_cast<size_t>(ne));

        for (Index i = 0; i < nv; ++i)
          outVerts.push_back(mesh.getVertexCoordinates(i));

        auto addVertex = [&](const Point& p) -> Index
        {
          outVerts.push_back(p);
          return static_cast<Index>(outVerts.getCount() - 1);
        };

        // --------------------------------------------------------------------
        // Local geometry helpers
        // --------------------------------------------------------------------
        auto signedVolume = [&](Index a, Index b, Index c, Index d) -> long double
        {
          const auto A = outVerts[a];
          const auto B = outVerts[b];
          const auto C = outVerts[c];
          const auto D = outVerts[d];
          return (long double) (B - A).dot((C - A).cross(D - A));
        };

        auto tetAbsVolume = [&](const Tet& t) -> long double
        {
          return std::abs(signedVolume(t[0], t[1], t[2], t[3])) / 6.0L;
        };

        auto triPush = [](std::vector<Tri>& out, Index a, Index b, Index c)
        {
          if (a == b || a == c || b == c)
            return;
          out.push_back({{ a, b, c }});
        };

        auto sortedTri = [](Tri t) -> Tri
        {
          std::sort(t.begin(), t.end());
          return t;
        };

        auto normalizeTriList = [&](const std::vector<Tri>& tris)
        {
          std::vector<Tri> out;
          out.reserve(tris.size());
          for (const auto& tri : tris)
            out.push_back(sortedTri(tri));
          std::sort(out.begin(), out.end());
          return out;
        };

        auto extractBoundary = [&](const std::vector<Tet>& tets)
        {
          std::map<Tri, int> count;
          auto addFace = [&](Index a, Index b, Index c)
          {
            count[sortedTri(Tri{{a, b, c}})]++;
          };

          for (const auto& t : tets)
          {
            addFace(t[0], t[1], t[2]);
            addFace(t[0], t[1], t[3]);
            addFace(t[0], t[2], t[3]);
            addFace(t[1], t[2], t[3]);
          }

          std::vector<Tri> bdry;
          for (const auto& [tri, n] : count)
            if (n == 1)
              bdry.push_back(tri);

          std::sort(bdry.begin(), bdry.end());
          return bdry;
        };

        struct TetQualityStats
        {
          long double qmin = -1.0L;
          long double qavg = -1.0L;
          long double qmax = -1.0L;
        };

        auto candidateQuality = [&](const std::vector<Tet>& tets) -> TetQualityStats
        {
          TetQualityStats qs;
          qs.qmin = std::numeric_limits<long double>::infinity();
          qs.qmax = 0.0L;
          qs.qavg = 0.0L;

          for (const auto& t : tets)
          {
            const long double q = m_quality_metric(outVerts, t);
            if (!(q > 0.0L) || !std::isfinite((double) q))
              return TetQualityStats{};

            qs.qmin = std::min(qs.qmin, q);
            qs.qmax = std::max(qs.qmax, q);
            qs.qavg += q;
          }

          qs.qavg /= (long double) tets.size();
          return qs;
        };

        struct InterfaceStats
        {
          long double area = 0.0L;
          long double minEdge = 0.0L;
          long double minTriArea = 0.0L;
          Point centroid{3};
          Point normal{3};
          int triCount = 0;
        };

        auto edgeLen = [&](Index a, Index b) -> long double
        {
          return (long double) (outVerts[b] - outVerts[a]).norm();
        };

        auto triangleArea = [&](Index a, Index b, Index c) -> long double
        {
          const auto A = outVerts[a];
          const auto B = outVerts[b];
          const auto C = outVerts[c];
          return (long double) ((B - A).cross(C - A)).norm() / 2.0L;
        };

        auto interfaceStats = [&](const std::vector<Tri>& tris) -> InterfaceStats
        {
          InterfaceStats s;
          s.centroid.setZero();
          s.normal.setZero();
          s.minEdge = std::numeric_limits<long double>::infinity();
          s.minTriArea = std::numeric_limits<long double>::infinity();
          s.triCount = static_cast<int>(tris.size());
          if (tris.empty())
          {
            s.minEdge = 0.0L;
            s.minTriArea = 0.0L;
            return s;
          }

          for (const auto& t : tris)
          {
            const auto A = outVerts[t[0]];
            const auto B = outVerts[t[1]];
            const auto C = outVerts[t[2]];
            const auto n = (B - A).cross(C - A);
            const long double area2 = (long double) n.norm();
            const long double area = area2 / 2.0L;
            s.area += area;
            s.normal += n;
            s.centroid += (A + B + C) * (area / 3.0L);
            s.minTriArea = std::min(s.minTriArea, area);
            s.minEdge = std::min(s.minEdge, edgeLen(t[0], t[1]));
            s.minEdge = std::min(s.minEdge, edgeLen(t[1], t[2]));
            s.minEdge = std::min(s.minEdge, edgeLen(t[2], t[0]));
          }

          if (s.area > 0.0L)
            s.centroid /= s.area;
          else
            s.centroid.setZero();

          if (s.minEdge == std::numeric_limits<long double>::infinity())
            s.minEdge = 0.0L;
          if (s.minTriArea == std::numeric_limits<long double>::infinity())
            s.minTriArea = 0.0L;
          return s;
        };

        struct CandidateScore
        {
          long double qmin = -1.0L;
          long double qavg = -1.0L;
          long double qualitySpreadPenalty = std::numeric_limits<long double>::infinity();
          long double fidelityPenalty = std::numeric_limits<long double>::infinity();
          long double shortEdgePenalty = std::numeric_limits<long double>::infinity();
          long double smallTriPenalty = std::numeric_limits<long double>::infinity();
          long double jaggedPenalty = std::numeric_limits<long double>::infinity();
          long double neighborPenalty = std::numeric_limits<long double>::infinity();
        };

        auto betterScore = [&](const CandidateScore& a, const CandidateScore& b)
        {
          static constexpr long double tol = 1e-18L;
          // Lexicographic policy:
          //   1) enforce minimum quality strongly (higher qmin, then qavg),
          //   2) preserve interface fidelity/smoothness among admissible candidates.
          if (a.qmin > b.qmin + tol) return true;
          if (b.qmin > a.qmin + tol) return false;
          if (a.qavg > b.qavg + tol) return true;
          if (b.qavg > a.qavg + tol) return false;
          if (a.qualitySpreadPenalty + tol < b.qualitySpreadPenalty) return true;
          if (b.qualitySpreadPenalty + tol < a.qualitySpreadPenalty) return false;
          if (a.fidelityPenalty + tol < b.fidelityPenalty) return true;
          if (b.fidelityPenalty + tol < a.fidelityPenalty) return false;
          if (a.shortEdgePenalty + tol < b.shortEdgePenalty) return true;
          if (b.shortEdgePenalty + tol < a.shortEdgePenalty) return false;
          if (a.smallTriPenalty + tol < b.smallTriPenalty) return true;
          if (b.smallTriPenalty + tol < a.smallTriPenalty) return false;
          if (a.jaggedPenalty + tol < b.jaggedPenalty) return true;
          if (b.jaggedPenalty + tol < a.jaggedPenalty) return false;
          if (a.neighborPenalty + tol < b.neighborPenalty) return true;
          return false;
        };

        auto triArea2 = [&](Index a, Index b, Index c) -> long double
        {
          const auto A = outVerts[a];
          const auto B = outVerts[b];
          const auto C = outVerts[c];
          return (long double) ((B - A).cross(C - A)).norm();
        };

        auto sortedPair = [](Index i, Index j) -> std::pair<Index, Index>
        {
          if (i > j) std::swap(i, j);
          return { i, j };
        };

        auto triangulateQuadFace =
          [&](Index q0, Index q1, Index q2, Index q3, std::vector<Tri>& out)
        {
          const long double m02 =
            std::min(triArea2(q0, q1, q2), triArea2(q0, q2, q3));
          const long double m13 =
            std::min(triArea2(q1, q2, q3), triArea2(q1, q3, q0));

          if (m02 > m13)
          {
            triPush(out, q0, q1, q2);
            triPush(out, q0, q2, q3);
            return;
          }

          if (m13 > m02)
          {
            triPush(out, q1, q2, q3);
            triPush(out, q1, q3, q0);
            return;
          }

          const auto d02 = sortedPair(q0, q2);
          const auto d13 = sortedPair(q1, q3);
          if (d02 < d13)
          {
            triPush(out, q0, q1, q2);
            triPush(out, q0, q2, q3);
          }
          else
          {
            triPush(out, q1, q2, q3);
            triPush(out, q1, q3, q0);
          }
        };

        // --------------------------------------------------------------------
        // Edge cuts
        // --------------------------------------------------------------------
        struct EdgeCut
        {
          enum class Kind : uint8_t
          {
            None,
            EndpointA,
            EndpointB,
            Interior,
            Invalid
          };

          Kind kind = Kind::None;
          Index vertex = INVALID;
        };

        std::vector<EdgeCut> edgeCuts(static_cast<size_t>(ne));

        for (Index eid = 0; eid < ne; ++eid)
        {
          const auto& ev = conn.getIncidence({1, 0}, eid);
          const Index a = ev[0];
          const Index b = ev[1];

          const Real fa = ls[a];
          const Real fb = ls[b];
          const Sign sa = vsgn[(size_t)a];
          const Sign sb = vsgn[(size_t)b];

          EdgeCut ec;

          if (sa == Sign::Invalid || sb == Sign::Invalid)
          {
            ec.kind = EdgeCut::Kind::Invalid;
            edgeCuts[(size_t)eid] = ec;
            continue;
          }

          if (sa == Sign::Zero)
          {
            ec.kind   = EdgeCut::Kind::EndpointA;
            ec.vertex = a;
            edgeCuts[(size_t)eid] = ec;
            continue;
          }

          if (sb == Sign::Zero)
          {
            ec.kind   = EdgeCut::Kind::EndpointB;
            ec.vertex = b;
            edgeCuts[(size_t)eid] = ec;
            continue;
          }

          if (sa == sb)
          {
            ec.kind = EdgeCut::Kind::None;
            edgeCuts[(size_t)eid] = ec;
            continue;
          }

          const Real denom = fa - fb;
          if (!finite(denom) || std::abs(denom) <= eps_snap)
          {
            if (std::abs(fa) <= std::abs(fb))
            {
              ec.kind   = EdgeCut::Kind::EndpointA;
              ec.vertex = a;
            }
            else
            {
              ec.kind   = EdgeCut::Kind::EndpointB;
              ec.vertex = b;
            }
            edgeCuts[(size_t)eid] = ec;
            continue;
          }

          Real t = fa / denom;
          t = std::max(Real(0), std::min(Real(1), t));

          if (t <= eps_snap)
          {
            ec.kind   = EdgeCut::Kind::EndpointA;
            ec.vertex = a;
            edgeCuts[(size_t)eid] = ec;
            continue;
          }

          if (t >= Real(1) - eps_snap)
          {
            ec.kind   = EdgeCut::Kind::EndpointB;
            ec.vertex = b;
            edgeCuts[(size_t)eid] = ec;
            continue;
          }

          const auto pa = outVerts[a];
          const auto pb = outVerts[b];

          ec.kind   = EdgeCut::Kind::Interior;
          ec.vertex = addVertex(pa + t * (pb - pa));
          edgeCuts[(size_t)eid] = ec;
        }

        auto getEdgeCutVertex = [&](Index eid) -> Index
        {
          const EdgeCut& ec = edgeCuts[(size_t)eid];
          switch (ec.kind)
          {
            case EdgeCut::Kind::EndpointA:
            case EdgeCut::Kind::EndpointB:
            case EdgeCut::Kind::Interior:
              return ec.vertex;
            default:
              break;
          }

          Alert::MemberFunctionException(*this, __func__)
            << "Requested cut vertex on non-cut edge " << eid << "."
            << Alert::Raise;
          return INVALID;
        };

        auto findFaceEdge = [&](Index fid, Index x, Index y) -> Index
        {
          const auto& feids = conn.getIncidence({2, 1}, fid);
          for (Index eid : feids)
          {
            const auto& ev = conn.getIncidence({1, 0}, eid);
            if ((ev[0] == x && ev[1] == y) || (ev[0] == y && ev[1] == x))
              return eid;
          }

          Alert::MemberFunctionException(*this, __func__)
            << "Could not find edge (" << x << ", " << y << ") in face " << fid << "."
            << Alert::Raise;
          return INVALID;
        };

        // --------------------------------------------------------------------
        // Face restrictions
        // --------------------------------------------------------------------
        struct FaceRestriction
        {
          enum class Kind : uint8_t
          {
            WholeNegative,
            WholePositive,
            Split,
            Invalid
          };

          Kind kind = Kind::Invalid;
          std::vector<Tri> neg;
          std::vector<Tri> pos;
        };

        std::vector<FaceRestriction> faceRestr(static_cast<size_t>(nf));

        for (Index fid = 0; fid < nf; ++fid)
        {
          FaceRestriction fr;

          const auto& fv = conn.getPolytope(2, fid);
          const Index a = fv(0);
          const Index b = fv(1);
          const Index c = fv(2);

          const Sign sa = vsgn[(size_t)a];
          const Sign sb = vsgn[(size_t)b];
          const Sign sc = vsgn[(size_t)c];

          if (sa == Sign::Invalid || sb == Sign::Invalid || sc == Sign::Invalid)
          {
            fr.kind = FaceRestriction::Kind::Invalid;
            faceRestr[(size_t)fid] = std::move(fr);
            continue;
          }

          const int ia = signToInt(sa);
          const int ib = signToInt(sb);
          const int ic = signToInt(sc);

          const int nneg = (ia < 0) + (ib < 0) + (ic < 0);
          const int npos = (ia > 0) + (ib > 0) + (ic > 0);
          const int nzer = (ia == 0) + (ib == 0) + (ic == 0);

          if (nneg == 0 && npos == 0 && nzer == 3)
          {
            fr.kind = FaceRestriction::Kind::Split;
            triPush(fr.neg, a, b, c);
            triPush(fr.pos, a, b, c);
            faceRestr[(size_t)fid] = std::move(fr);
            continue;
          }

          if (npos == 0)
          {
            fr.kind = FaceRestriction::Kind::WholeNegative;
            triPush(fr.neg, a, b, c);
            faceRestr[(size_t)fid] = std::move(fr);
            continue;
          }

          if (nneg == 0)
          {
            fr.kind = FaceRestriction::Kind::WholePositive;
            triPush(fr.pos, a, b, c);
            faceRestr[(size_t)fid] = std::move(fr);
            continue;
          }

          fr.kind = FaceRestriction::Kind::Split;

          auto cutOn = [&](Index x, Index y) -> Index
          {
            const Index eid = findFaceEdge(fid, x, y);
            return getEdgeCutVertex(eid);
          };

          auto handleOneTwo = [&](Index lone, Index s0, Index s1, bool loneIsNeg)
          {
            const Index p0 = cutOn(lone, s0);
            const Index p1 = cutOn(lone, s1);

            if (loneIsNeg)
            {
              triPush(fr.neg, lone, p0, p1);
              triangulateQuadFace(s0, s1, p1, p0, fr.pos);
            }
            else
            {
              triPush(fr.pos, lone, p0, p1);
              triangulateQuadFace(s0, s1, p1, p0, fr.neg);
            }
          };

          if (nzer == 1 && nneg == 1 && npos == 1)
          {
            Index z = INVALID, n = INVALID, p = INVALID;
            if (ia == 0) z = a; else if (ia < 0) n = a; else p = a;
            if (ib == 0) z = b; else if (ib < 0) n = b; else p = b;
            if (ic == 0) z = c; else if (ic < 0) n = c; else p = c;

            const Index q = cutOn(n, p);
            triPush(fr.neg, z, n, q);
            triPush(fr.pos, z, q, p);

            faceRestr[(size_t)fid] = std::move(fr);
            continue;
          }

          if (nneg == 1 && npos == 2)
          {
            if (ia < 0)      handleOneTwo(a, b, c, true);
            else if (ib < 0) handleOneTwo(b, a, c, true);
            else             handleOneTwo(c, a, b, true);

            faceRestr[(size_t)fid] = std::move(fr);
            continue;
          }

          if (nneg == 2 && npos == 1)
          {
            if (ia > 0)      handleOneTwo(a, b, c, false);
            else if (ib > 0) handleOneTwo(b, a, c, false);
            else             handleOneTwo(c, a, b, false);

            faceRestr[(size_t)fid] = std::move(fr);
            continue;
          }

          if (npos == 0)
          {
            fr.kind = FaceRestriction::Kind::WholeNegative;
            fr.neg.clear();
            fr.pos.clear();
            triPush(fr.neg, a, b, c);
            faceRestr[(size_t)fid] = std::move(fr);
            continue;
          }

          if (nneg == 0)
          {
            fr.kind = FaceRestriction::Kind::WholePositive;
            fr.neg.clear();
            fr.pos.clear();
            triPush(fr.pos, a, b, c);
            faceRestr[(size_t)fid] = std::move(fr);
            continue;
          }

          fr.kind = FaceRestriction::Kind::Invalid;
          faceRestr[(size_t)fid] = std::move(fr);
        }

        // --------------------------------------------------------------------
        // Cell patterns
        // --------------------------------------------------------------------
        struct CellPattern
        {
          enum class Kind : uint8_t
          {
            WholeNegative,
            WholePositive,
            Split13,
            Split31,
            Split22,
            Invalid
          };

          Kind kind = Kind::Invalid;
          Optional<Attribute> negAttr;
          Optional<Attribute> posAttr;
        };

        std::vector<CellPattern> cellPatterns(static_cast<size_t>(mesh.getCellCount()));

        for (auto cit = mesh.getCell(); !cit.end(); ++cit)
        {
          const Index cid = cit->getIndex();

          if (cit->getGeometry() != Polytope::Type::Tetrahedron)
            Alert::MemberFunctionException(*this, __func__)
              << "Only tetrahedral meshes are supported."
              << Alert::Raise;

          const Optional<Attribute> cellAttr = toOptionalAttr(cit->getAttribute());
          const auto& cv = cit->getVertices();

          std::array<Sign, 4> ss{
            vsgn[(size_t)cv(0)],
            vsgn[(size_t)cv(1)],
            vsgn[(size_t)cv(2)],
            vsgn[(size_t)cv(3)]
          };

          CellPattern cp;

          if (ss[0] == Sign::Invalid || ss[1] == Sign::Invalid ||
              ss[2] == Sign::Invalid || ss[3] == Sign::Invalid)
          {
            cp.kind = CellPattern::Kind::Invalid;
            cellPatterns[(size_t)cid] = cp;
            continue;
          }

          const int nneg =
            (ss[0] == Sign::Negative) + (ss[1] == Sign::Negative) +
            (ss[2] == Sign::Negative) + (ss[3] == Sign::Negative);

          const int npos =
            (ss[0] == Sign::Positive) + (ss[1] == Sign::Positive) +
            (ss[2] == Sign::Positive) + (ss[3] == Sign::Positive);

          const bool strictCut = (nneg > 0 && npos > 0);

          if (!strictCut)
          {
            if (npos == 0)
            {
              cp.kind    = CellPattern::Kind::WholeNegative;
              cp.negAttr = inferLabel(3, cellAttr, true);
            }
            else
            {
              cp.kind    = CellPattern::Kind::WholePositive;
              cp.posAttr = inferLabel(3, cellAttr, false);
            }
          }
          else
          {
            cp.negAttr = inferLabel(3, cellAttr, true);
            cp.posAttr = inferLabel(3, cellAttr, false);

            if (nneg == 1)      cp.kind = CellPattern::Kind::Split13;
            else if (nneg == 3) cp.kind = CellPattern::Kind::Split31;
            else if (nneg == 2) cp.kind = CellPattern::Kind::Split22;
            else                cp.kind = CellPattern::Kind::Invalid;
          }

          cellPatterns[(size_t)cid] = cp;
        }

        for (Index fid = 0; fid < nf; ++fid)
        {
          const auto& fr = faceRestr[(size_t)fid];
          if (fr.kind == FaceRestriction::Kind::Invalid)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Invalid face restriction on face " << fid << "."
              << Alert::Raise;
          }
        }

        // --------------------------------------------------------------------
        // Output cells
        // --------------------------------------------------------------------
        static constexpr int8_t SidePositive = 0;
        static constexpr int8_t SideNegative = 1;

        std::vector<Tet> outCells;
        std::vector<Optional<Attribute>> outCellAttr;
        std::vector<int8_t> outCellSide;

        outCells.reserve(static_cast<size_t>(mesh.getCellCount()) * 8);
        outCellAttr.reserve(outCells.capacity());
        outCellSide.reserve(outCells.capacity());

        auto emitTet = [&](Index a, Index b, Index c, Index d,
                           const Optional<Attribute>& attr,
                           int8_t side)
        {
          if (a == b || a == c || a == d || b == c || b == d || c == d)
            return;

          Tet t{{ a, b, c, d }};
          if (signedVolume(t[0], t[1], t[2], t[3]) < (long double) 0)
            std::swap(t[1], t[2]);

          outCells.push_back(t);
          outCellAttr.push_back(attr);
          outCellSide.push_back(side);
        };

        auto emitTetList = [&](const std::vector<Tet>& tets,
                               const Optional<Attribute>& attr,
                               int8_t side)
        {
          for (const auto& t : tets)
            emitTet(t[0], t[1], t[2], t[3], attr, side);
        };

        auto conePolyhedron =
          [&](const std::vector<Tri>& boundaryTris,
              const Optional<Attribute>& attr,
              int8_t side)
        {
          if (boundaryTris.empty())
            return;

          std::unordered_set<Index> uniq;
          uniq.reserve(boundaryTris.size() * 2);

          Point c(3);
          c.setZero();

          size_t cnt = 0;
          for (const auto& tri : boundaryTris)
          {
            for (int k = 0; k < 3; ++k)
            {
              const Index vi = tri[(size_t)k];
              if (uniq.insert(vi).second)
              {
                c += outVerts[vi];
                cnt++;
              }
            }
          }

          if (cnt == 0)
            return;

          c /= Real(cnt);
          const Index center = addVertex(c);

          for (const auto& tri : boundaryTris)
            emitTet(center, tri[0], tri[1], tri[2], attr, side);
        };

        auto buildLocalEdgeId = [&](Index cid, const std::array<Index, 4>& v, Index localEdgeId[4][4])
        {
          for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
              localEdgeId[i][j] = INVALID;

          const auto& ceids = conn.getIncidence({3, 1}, cid);
          for (Index eid : ceids)
          {
            const auto& ev = conn.getIncidence({1, 0}, eid);
            int ia = -1, ib = -1;
            for (int i = 0; i < 4; ++i)
            {
              if (v[(size_t)i] == ev[0]) ia = i;
              if (v[(size_t)i] == ev[1]) ib = i;
            }

            if (ia >= 0 && ib >= 0)
            {
              localEdgeId[ia][ib] = eid;
              localEdgeId[ib][ia] = eid;
            }
          }
        };

        auto cutLocal = [&](Index localEdgeId[4][4], int i, int j) -> Index
        {
          const Index eid = localEdgeId[i][j];
          if (eid == INVALID)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Missing local edge (" << i << ", " << j << ") in tetrahedron."
              << Alert::Raise;
          }
          return getEdgeCutVertex(eid);
        };

        auto prismVolumeFromReferenceFill =
          [&](Index u0, Index u1, Index u2,
              Index v0, Index v1, Index v2) -> long double
        {
          const std::array<Tet, 3> ref{{
            Tet{{u0, u1, u2, v0}},
            Tet{{u1, u2, v0, v1}},
            Tet{{u2, v0, v1, v2}}
          }};

          long double vol = 0.0L;
          for (const auto& t : ref)
            vol += tetAbsVolume(t);
          return vol;
        };

        auto tryBestPrismFillExhaustive =
          [&](Index u0, Index u1, Index u2,
              Index v0, Index v1, Index v2,
              const std::vector<Tri>& boundary,
              std::vector<Tet>& bestTets,
              CandidateScore& bestScore,
              const std::vector<Tri>& candidateIface,
              const std::vector<Tri>& referenceIface,
              long double neighborPenalty,
              long double minQuality) -> bool
        {
          bestTets.clear();
          bestScore = CandidateScore{};
          bestScore.qmin = -1.0L;
          bestScore.qavg = -1.0L;

          const auto required = normalizeTriList(boundary);
          const std::array<Index, 6> P{{u0, u1, u2, v0, v1, v2}};
          const long double prismVol =
            prismVolumeFromReferenceFill(u0, u1, u2, v0, v1, v2);

          static constexpr long double volTol = 1e-18L;

          std::vector<Tet> allTets;
          allTets.reserve(15);

          for (int i = 0; i < 6; ++i)
          {
            for (int j = i + 1; j < 6; ++j)
            {
              for (int k = j + 1; k < 6; ++k)
              {
                for (int l = k + 1; l < 6; ++l)
                {
                  Tet t{{P[i], P[j], P[k], P[l]}};
                  if (tetAbsVolume(t) > volTol)
                    allTets.push_back(t);
                }
              }
            }
          }

          auto tryCandidate = [&](const std::vector<Tet>& cand)
          {
            long double sumVol = 0.0L;
            for (const auto& t : cand)
            {
              const long double v = tetAbsVolume(t);
              if (!(v > volTol))
                return;
              sumVol += v;
            }

            if (std::abs(sumVol - prismVol) > 1e-12L * std::max<long double>(1.0L, prismVol))
              return;

            const auto bdry = extractBoundary(cand);
            if (bdry != required)
              return;

            std::array<bool, 6> used{{false, false, false, false, false, false}};
            auto mark = [&](Index x)
            {
              for (int q = 0; q < 6; ++q)
                if (P[q] == x)
                  used[(size_t)q] = true;
            };

            for (const auto& t : cand)
            {
              mark(t[0]); mark(t[1]); mark(t[2]); mark(t[3]);
            }

            if (!std::all_of(used.begin(), used.end(), [](bool x) { return x; }))
              return;

            const auto q = candidateQuality(cand);
            if (q.qmin < minQuality)
              return;

            const auto refStats  = interfaceStats(referenceIface);
            const auto candStats = interfaceStats(candidateIface);

            CandidateScore score;
            score.qmin = q.qmin;
            score.qavg = q.qavg;
            score.qualitySpreadPenalty = (q.qmin > 0.0L) ? (q.qmax / q.qmin) : 1e30L;
            score.fidelityPenalty =
              std::abs(candStats.area - refStats.area) +
              (long double) (candStats.centroid - refStats.centroid).norm() +
              std::abs((long double) candStats.normal.norm() - (long double) refStats.normal.norm());
            const long double ifaceScale = std::max<long double>(1e-18L, std::sqrt(std::max<long double>(candStats.area, refStats.area)));
            score.shortEdgePenalty = (candStats.minEdge > 0.0L) ? (ifaceScale / candStats.minEdge) : 1e30L;
            score.smallTriPenalty = (candStats.minTriArea > 0.0L) ? ((ifaceScale * ifaceScale) / candStats.minTriArea) : 1e30L;
            score.jaggedPenalty = (candStats.triCount > 0 && candStats.area > 0.0L)
              ? std::abs((long double) candStats.normal.norm() / (2.0L * candStats.area) - 1.0L)
              : 1e30L;
            score.neighborPenalty = neighborPenalty;

            if (bestScore.qmin < 0.0L || betterScore(score, bestScore))
            {
              bestScore = score;
              bestTets  = cand;
            }
          };

          const int m = static_cast<int>(allTets.size());
          for (int i = 0; i < m; ++i)
          {
            for (int j = i + 1; j < m; ++j)
            {
              for (int k = j + 1; k < m; ++k)
              {
                const std::vector<Tet> cand{
                  allTets[(size_t)i],
                  allTets[(size_t)j],
                  allTets[(size_t)k]
                };
                tryCandidate(cand);
              }
            }
          }

          return !bestTets.empty();
        };

        for (auto cit = mesh.getCell(); !cit.end(); ++cit)
        {
          const Index cid = cit->getIndex();
          const auto& cp = cellPatterns[(size_t)cid];
          const auto& cv = cit->getVertices();
          const std::array<Index, 4> v{{ cv(0), cv(1), cv(2), cv(3) }};

          switch (cp.kind)
          {
            case CellPattern::Kind::WholeNegative:
              emitTet(v[0], v[1], v[2], v[3], cp.negAttr, SideNegative);
              break;

            case CellPattern::Kind::WholePositive:
              emitTet(v[0], v[1], v[2], v[3], cp.posAttr, SidePositive);
              break;

            case CellPattern::Kind::Invalid:
              Alert::MemberFunctionException(*this, __func__)
                << "Invalid cell pattern on tetrahedron " << cid << "."
                << Alert::Raise;
              break;

            case CellPattern::Kind::Split13:
            case CellPattern::Kind::Split31:
            {
              std::vector<Tri> negFaceBdry, posFaceBdry;
              negFaceBdry.reserve(8);
              posFaceBdry.reserve(8);

              const auto& cfids = conn.getIncidence({3, 2}, cid);
              for (Index fid : cfids)
              {
                const auto& fr = faceRestr[(size_t)fid];
                negFaceBdry.insert(negFaceBdry.end(), fr.neg.begin(), fr.neg.end());
                posFaceBdry.insert(posFaceBdry.end(), fr.pos.begin(), fr.pos.end());
              }

              Index localEdgeId[4][4];
              buildLocalEdgeId(cid, v, localEdgeId);

              std::array<int, 4> ss{
                signToInt(vsgn[(size_t)v[0]]),
                signToInt(vsgn[(size_t)v[1]]),
                signToInt(vsgn[(size_t)v[2]]),
                signToInt(vsgn[(size_t)v[3]])
              };

              const int nneg = (ss[0] < 0) + (ss[1] < 0) + (ss[2] < 0) + (ss[3] < 0);
              const bool loneIsNeg = (nneg == 1);

              int lone = -1;
              int others[3];
              int k = 0;
              for (int i = 0; i < 4; ++i)
              {
                const bool isNeg = (ss[i] < 0);
                if (isNeg == loneIsNeg) lone = i;
                else others[k++] = i;
              }

              const Index p0 = cutLocal(localEdgeId, lone, others[0]);
              const Index p1 = cutLocal(localEdgeId, lone, others[1]);
              const Index p2 = cutLocal(localEdgeId, lone, others[2]);

              Tri iface{{ p0, p1, p2 }};

              const long double loneQ =
                m_quality_metric(outVerts, Tet{{v[(size_t)lone], p0, p1, p2}});
              if (!(loneQ >= this->getMinimumQuality()))
              {
                std::vector<Tri> negBoundary = negFaceBdry;
                std::vector<Tri> posBoundary = posFaceBdry;
                triPush(negBoundary, iface[0], iface[1], iface[2]);
                triPush(posBoundary, iface[0], iface[1], iface[2]);
                conePolyhedron(negBoundary, cp.negAttr, SideNegative);
                conePolyhedron(posBoundary, cp.posAttr, SidePositive);
                break;
              }

              if (loneIsNeg)
              {
                std::vector<Tri> posBoundary = posFaceBdry;
                triPush(posBoundary, iface[0], iface[1], iface[2]);

                const Index u0 = v[(size_t)others[0]];
                const Index u1 = v[(size_t)others[1]];
                const Index u2 = v[(size_t)others[2]];

                std::vector<Tet> best;
                CandidateScore score;
                std::vector<Tri> exactIface;
                triPush(exactIface, iface[0], iface[1], iface[2]);
                bool ok = tryBestPrismFillExhaustive(
                  u0, u1, u2, p0, p1, p2, posBoundary, best, score, exactIface,
                  exactIface,
                  0.0L, this->getMinimumQuality());
                if (!ok)
                {
                  // Tier-2 fallback before coning: keep exact topology but relax
                  // the strict quality gate to avoid hanging entities.
                  ok = tryBestPrismFillExhaustive(
                    u0, u1, u2, p0, p1, p2, posBoundary, best, score, exactIface,
                    exactIface,
                    1e-6L, 0.0L);
                }
                if (ok)
                {
                  emitTet(v[(size_t)lone], p0, p1, p2, cp.negAttr, SideNegative);
                  emitTetList(best, cp.posAttr, SidePositive);
                }
                else
                  conePolyhedron(posBoundary, cp.posAttr, SidePositive);
              }
              else
              {
                std::vector<Tri> negBoundary = negFaceBdry;
                triPush(negBoundary, iface[0], iface[1], iface[2]);

                const Index u0 = v[(size_t)others[0]];
                const Index u1 = v[(size_t)others[1]];
                const Index u2 = v[(size_t)others[2]];

                std::vector<Tet> best;
                CandidateScore score;
                std::vector<Tri> exactIface;
                triPush(exactIface, iface[0], iface[1], iface[2]);
                bool ok = tryBestPrismFillExhaustive(
                  u0, u1, u2, p0, p1, p2, negBoundary, best, score, exactIface,
                  exactIface,
                  0.0L, this->getMinimumQuality());
                if (!ok)
                {
                  ok = tryBestPrismFillExhaustive(
                    u0, u1, u2, p0, p1, p2, negBoundary, best, score, exactIface,
                    exactIface,
                    1e-6L, 0.0L);
                }
                if (ok)
                {
                  emitTet(v[(size_t)lone], p0, p1, p2, cp.posAttr, SidePositive);
                  emitTetList(best, cp.negAttr, SideNegative);
                }
                else
                  conePolyhedron(negBoundary, cp.negAttr, SideNegative);
              }
              break;
            }

            case CellPattern::Kind::Split22:
            {
              std::vector<Tri> negFaceBdry, posFaceBdry;
              negFaceBdry.reserve(10);
              posFaceBdry.reserve(10);

              const auto& cfids = conn.getIncidence({3, 2}, cid);
              for (Index fid : cfids)
              {
                const auto& fr = faceRestr[(size_t)fid];
                negFaceBdry.insert(negFaceBdry.end(), fr.neg.begin(), fr.neg.end());
                posFaceBdry.insert(posFaceBdry.end(), fr.pos.begin(), fr.pos.end());
              }

              Index localEdgeId[4][4];
              buildLocalEdgeId(cid, v, localEdgeId);

              std::array<int, 4> ss{
                signToInt(vsgn[(size_t)v[0]]),
                signToInt(vsgn[(size_t)v[1]]),
                signToInt(vsgn[(size_t)v[2]]),
                signToInt(vsgn[(size_t)v[3]])
              };

              int ineg[2], ipos[2], kn = 0, kp = 0;
              for (int i = 0; i < 4; ++i)
                (ss[i] < 0) ? (ineg[kn++] = i) : (ipos[kp++] = i);

              const Index n0 = v[(size_t)ineg[0]];
              const Index n1 = v[(size_t)ineg[1]];
              const Index p0 = v[(size_t)ipos[0]];
              const Index p1 = v[(size_t)ipos[1]];

              const Index a = cutLocal(localEdgeId, ineg[0], ipos[0]);
              const Index b = cutLocal(localEdgeId, ineg[0], ipos[1]);
              const Index c = cutLocal(localEdgeId, ineg[1], ipos[0]);
              const Index d = cutLocal(localEdgeId, ineg[1], ipos[1]);

              std::array<std::vector<Tri>, 3> ifaceOptions;
              triPush(ifaceOptions[0], a, b, d);
              triPush(ifaceOptions[0], a, d, c);

              triPush(ifaceOptions[1], a, b, c);
              triPush(ifaceOptions[1], b, d, c);

              // Simplified tier: collapse fragile quads to a single triangle.
              // This sacrifices local fidelity but removes sliver-producing facets.
              const long double dAC = edgeLen(a, c);
              const long double dBD = edgeLen(b, d);
              if (dAC <= dBD)
                triPush(ifaceOptions[2], a, b, d);
              else
                triPush(ifaceOptions[2], a, c, d);

              bool anyValid = false;
              std::vector<Tet> bestNeg, bestPos;
              CandidateScore bestScore{};
              bestScore.qmin = -1.0L;
              bestScore.qavg = -1.0L;

              for (int opt = 0; opt < 3; ++opt)
              {
                std::vector<Tri> negBoundary = negFaceBdry;
                std::vector<Tri> posBoundary = posFaceBdry;
                negBoundary.insert(negBoundary.end(),
                                   ifaceOptions[opt].begin(), ifaceOptions[opt].end());
                posBoundary.insert(posBoundary.end(),
                                   ifaceOptions[opt].begin(), ifaceOptions[opt].end());

                std::vector<Tet> negCand;
                CandidateScore negScore;
                const long double coherencePenalty =
                  (opt == 0 ? 0.0L : (opt == 1 ? 1e-8L : 1e-5L));
                const bool negOK = tryBestPrismFillExhaustive(
                  n0, a, b,
                  n1, c, d,
                  negBoundary,
                  negCand,
                  negScore,
                  ifaceOptions[opt],
                  ifaceOptions[0],
                  coherencePenalty,
                  this->getMinimumQuality());

                std::vector<Tet> posCand;
                CandidateScore posScore;
                const bool posOK = tryBestPrismFillExhaustive(
                  p0, a, c,
                  p1, b, d,
                  posBoundary,
                  posCand,
                  posScore,
                  ifaceOptions[opt],
                  ifaceOptions[0],
                  coherencePenalty,
                  this->getMinimumQuality());

                if (negOK && posOK)
                {
                  CandidateScore cellScore;
                  cellScore.qmin = std::min(negScore.qmin, posScore.qmin);
                  cellScore.qavg = 0.5L * (negScore.qavg + posScore.qavg);
                  cellScore.qualitySpreadPenalty =
                    std::max(negScore.qualitySpreadPenalty, posScore.qualitySpreadPenalty);
                  cellScore.fidelityPenalty = negScore.fidelityPenalty + posScore.fidelityPenalty;
                  cellScore.shortEdgePenalty = negScore.shortEdgePenalty + posScore.shortEdgePenalty;
                  cellScore.smallTriPenalty = negScore.smallTriPenalty + posScore.smallTriPenalty;
                  cellScore.jaggedPenalty = negScore.jaggedPenalty + posScore.jaggedPenalty;
                  cellScore.neighborPenalty = negScore.neighborPenalty + posScore.neighborPenalty;

                  if (!anyValid || betterScore(cellScore, bestScore))
                  {
                    anyValid = true;
                    bestNeg = std::move(negCand);
                    bestPos = std::move(posCand);
                    bestScore = cellScore;
                  }
                }
              }

              if (!anyValid)
              {
                for (int opt = 0; opt < 3; ++opt)
                {
                  std::vector<Tri> negBoundary = negFaceBdry;
                  std::vector<Tri> posBoundary = posFaceBdry;
                  negBoundary.insert(negBoundary.end(),
                                     ifaceOptions[opt].begin(), ifaceOptions[opt].end());
                  posBoundary.insert(posBoundary.end(),
                                     ifaceOptions[opt].begin(), ifaceOptions[opt].end());

                  std::vector<Tet> negCand, posCand;
                  CandidateScore negScore, posScore;
                  const long double coherencePenalty =
                    (opt == 0 ? 1e-6L : (opt == 1 ? 1e-5L : 1e-4L));

                  const bool negOK = tryBestPrismFillExhaustive(
                    n0, a, b, n1, c, d,
                    negBoundary, negCand, negScore, ifaceOptions[opt],
                    ifaceOptions[0],
                    coherencePenalty, 0.0L);
                  const bool posOK = tryBestPrismFillExhaustive(
                    p0, a, c, p1, b, d,
                    posBoundary, posCand, posScore, ifaceOptions[opt],
                    ifaceOptions[0],
                    coherencePenalty, 0.0L);

                  if (negOK && posOK)
                  {
                    CandidateScore cellScore;
                    cellScore.qmin = std::min(negScore.qmin, posScore.qmin);
                    cellScore.qavg = 0.5L * (negScore.qavg + posScore.qavg);
                    cellScore.qualitySpreadPenalty =
                      std::max(negScore.qualitySpreadPenalty, posScore.qualitySpreadPenalty);
                    cellScore.fidelityPenalty = negScore.fidelityPenalty + posScore.fidelityPenalty;
                    cellScore.shortEdgePenalty = negScore.shortEdgePenalty + posScore.shortEdgePenalty;
                    cellScore.smallTriPenalty = negScore.smallTriPenalty + posScore.smallTriPenalty;
                    cellScore.jaggedPenalty = negScore.jaggedPenalty + posScore.jaggedPenalty;
                    cellScore.neighborPenalty = negScore.neighborPenalty + posScore.neighborPenalty;
                    if (!anyValid || betterScore(cellScore, bestScore))
                    {
                      anyValid = true;
                      bestNeg = std::move(negCand);
                      bestPos = std::move(posCand);
                      bestScore = cellScore;
                    }
                  }
                }
              }

              if (anyValid)
              {
                emitTetList(bestNeg, cp.negAttr, SideNegative);
                emitTetList(bestPos, cp.posAttr, SidePositive);
              }
              else
              {
                std::vector<Tri> negBoundary = negFaceBdry;
                std::vector<Tri> posBoundary = posFaceBdry;
                negBoundary.insert(negBoundary.end(),
                                   ifaceOptions[0].begin(), ifaceOptions[0].end());
                posBoundary.insert(posBoundary.end(),
                                   ifaceOptions[0].begin(), ifaceOptions[0].end());

                conePolyhedron(negBoundary, cp.negAttr, SideNegative);
                conePolyhedron(posBoundary, cp.posAttr, SidePositive);
              }
              break;
            }
          }
        }

        // --------------------------------------------------------------------
        // Build output mesh
        // --------------------------------------------------------------------
        typename MeshType::Builder builder;
        builder.initialize(3).nodes(outVerts.getCount());
        builder.setVertices(std::move(outVerts));
        builder.reserve(3, outCells.size());

        for (size_t i = 0; i < outCells.size(); ++i)
        {
          const auto& t = outCells[i];
          builder.tetrahedron(IndexArray{{ t[0], t[1], t[2], t[3] }});
          if (outCellAttr[i])
            builder.attribute({ 3, static_cast<Index>(i) }, *outCellAttr[i]);
        }

        auto out = builder.finalize();

        // --------------------------------------------------------------------
        // Output postprocessing: interface marking and attribute transfer
        // --------------------------------------------------------------------
        out.getConnectivity().compute(2, 3);
        out.getConnectivity().compute(1, 3);
        const auto& oconn = out.getConnectivity();

        const auto& ifaceOpt2 = this->getInterface(2);
        const auto& ifaceOpt1 = this->getInterface(1);

        auto isFittedZeroV = [&](Index vi) -> bool
        {
          return (vi < nv) && finite(ls[vi]) && (std::abs(ls[vi]) <= eps_sign_user);
        };

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
            if (ifaceOpt2)
              out.setAttribute({2, fidx}, *ifaceOpt2);
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
            if (ifaceOpt1)
              out.setAttribute({1, eidx}, *ifaceOpt1);
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

    private:
      static QualityMetric defaultQualityMetric()
      {
        return [](const PointCloud& verts, const Tet& t) -> long double
        {
          const auto& A = verts[t[0]];
          const auto& B = verts[t[1]];
          const auto& C = verts[t[2]];
          const auto& D = verts[t[3]];

          const auto AB = B - A;
          const auto AC = C - A;
          const auto AD = D - A;

          const long double vol6 =
            std::abs((long double) AB.dot(AC.cross(AD)));

          long double hmax2 = 0.0L;
          auto upd = [&](const auto& X, const auto& Y)
          {
            const auto d = Y - X;
            hmax2 = std::max(hmax2, (long double) d.squaredNorm());
          };

          upd(A, B); upd(A, C); upd(A, D);
          upd(B, C); upd(B, D); upd(C, D);

          if (!(hmax2 > 0.0L))
            return 0.0L;

          const long double hmax = std::sqrt(hmax2);
          return vol6 / (hmax * hmax * hmax);
        };
      }

      Real m_sign_tolerance;
      Real m_snap_tolerance;
      long double m_min_quality;
      QualityMetric m_quality_metric;

      std::array<Optional<Attribute>, 4> m_old;
      std::array<Optional<Attribute>, 4> m_fallback;
  };

  template <class Mesh, class Data>
  LevelSetDiscretizerTetrahedra(
    const Variational::GridFunction<Variational::P1<Real, Mesh>, Data>& ls)
    -> LevelSetDiscretizerTetrahedra<
         Variational::GridFunction<Variational::P1<Real, Mesh>, Data>>;
}

#endif
