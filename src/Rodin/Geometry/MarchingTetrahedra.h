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
#include <map>
#include <optional>
#include <tuple>
#include <type_traits>
#include <variant>
#include <vector>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Variational/GridFunction.h"

#include "Mesh.h"

namespace Rodin::Geometry
{
  class MarchingTetrahedra
  {
    public:
      struct NoSplitT {};
      static constexpr NoSplitT NoSplit{};

      struct Split
      {
        Attribute negative;
        Attribute positive;
      };

      using SplitMap = std::map<Attribute, std::variant<Split, NoSplitT>>;

      MarchingTetrahedra() = default;

      MarchingTetrahedra& setNegativeAttribute(Attribute attr)
      {
        m_negativeAttribute = attr;
        return *this;
      }

      MarchingTetrahedra& setPositiveAttribute(Attribute attr)
      {
        m_positiveAttribute = attr;
        return *this;
      }

      MarchingTetrahedra& setInterfaceAttribute(Attribute attr)
      {
        m_interfaceAttribute = attr;
        return *this;
      }

      MarchingTetrahedra& setCellSplit(const SplitMap& split)
      {
        m_cellSplit = split;
        return *this;
      }

      MarchingTetrahedra& setFaceSplit(const SplitMap& split)
      {
        m_faceSplit = split;
        return *this;
      }

      MarchingTetrahedra& setEdgeSplit(const SplitMap& split)
      {
        m_edgeSplit = split;
        return *this;
      }

      MarchingTetrahedra& splitCell(Attribute attr, const Split& split)
      {
        m_cellSplit[attr] = split;
        return *this;
      }

      MarchingTetrahedra& splitFace(Attribute attr, const Split& split)
      {
        m_faceSplit[attr] = split;
        return *this;
      }

      MarchingTetrahedra& splitEdge(Attribute attr, const Split& split)
      {
        m_edgeSplit[attr] = split;
        return *this;
      }

      MarchingTetrahedra& noSplitCell(Attribute attr)
      {
        m_cellSplit[attr] = NoSplit;
        return *this;
      }

      MarchingTetrahedra& noSplitFace(Attribute attr)
      {
        m_faceSplit[attr] = NoSplit;
        return *this;
      }

      MarchingTetrahedra& noSplitEdge(Attribute attr)
      {
        m_edgeSplit[attr] = NoSplit;
        return *this;
      }

      template <class FES, class Data>
      Mesh<Context::Local> discretize(const Variational::GridFunction<FES, Data>& ls) const
      {
        const auto& mesh = ls.getFiniteElementSpace().getMesh();

        if (ls.getDimension() != 1)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "MarchingTetrahedra expects a scalar level-set GridFunction."
            << Alert::Raise;
        }
        if (mesh.getDimension() != 3 || mesh.getSpaceDimension() != 3)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "MarchingTetrahedra expects a 3D volume mesh."
            << Alert::Raise;
        }

        using Point = Math::SpatialPoint;

        struct EdgeKey
        {
          Index a;
          Index b;
          bool operator<(const EdgeKey& other) const
          {
            return std::tie(a, b) < std::tie(other.a, other.b);
          }
        };

        struct FaceKey
        {
          Index a, b, c; // sorted
          bool operator<(const FaceKey& other) const
          {
            return std::tie(a, b, c) < std::tie(other.a, other.b, other.c);
          }
        };

        struct Tet
        {
          std::array<Index, 4> v;
          Attribute attr;
          Attribute sourceAttr; // input CELL attribute (kept for fallback)
          enum class Sign { Negative, Positive, Unknown } sign;
        };

        // ---------------------------------------------------------------------
        // Helpers
        // ---------------------------------------------------------------------
        auto isNegative = [](Real x) { return x < 0.0; };

        auto makeFaceKey = [](Index i, Index j, Index k)
        {
          std::array<Index, 3> s{ i, j, k };
          std::sort(s.begin(), s.end());
          return FaceKey{ s[0], s[1], s[2] };
        };

        auto isNoSplit = [](Attribute attr, const SplitMap& splitMap)
        {
          if (splitMap.empty())
            return false;
          const auto it = splitMap.find(attr);
          if (it == splitMap.end())
            return false;
          return std::holds_alternative<NoSplitT>(it->second);
        };

        auto getSplitLabel = [&](Attribute baseAttr, typename Tet::Sign sign, const SplitMap& splitMap)
        {
          const auto getDefault = [&]()
          {
            if (sign == Tet::Sign::Negative)
              return m_negativeAttribute;
            if (sign == Tet::Sign::Positive)
              return m_positiveAttribute;
            return baseAttr;
          };

          if (splitMap.empty())
            return getDefault();

          const auto it = splitMap.find(baseAttr);
          if (it == splitMap.end())
            return baseAttr;

          return std::visit(
            [&](const auto& value) -> Attribute
            {
              using T = std::decay_t<decltype(value)>;
              if constexpr (std::is_same_v<T, NoSplitT>)
              {
                return baseAttr;
              }
              else
              {
                if (sign == Tet::Sign::Negative)
                  return value.negative;
                if (sign == Tet::Sign::Positive)
                  return value.positive;
                return baseAttr;
              }
            },
            it->second);
        };

        auto shouldSplitCell = [&](Attribute cellAttr)
        {
          // Your original policy:
          // - if m_cellSplit empty: split everything
          // - if no rule: do not split
          // - split iff rule is Split
          if (m_cellSplit.empty())
            return true;
          const auto it = m_cellSplit.find(cellAttr);
          if (it == m_cellSplit.end())
            return false;
          return std::holds_alternative<Split>(it->second);
        };

        auto edgePair = [](Index i, Index j)
        {
          return std::make_pair(std::min(i, j), std::max(i, j));
        };

        auto pushUnique = [](std::vector<Index>& v, Index x)
        {
          if (std::find(v.begin(), v.end(), x) == v.end())
            v.push_back(x);
        };

        auto triangulatePolygon = [&](const std::vector<Index>& polygon, std::vector<std::array<Index, 3>>& triangles)
        {
          if (polygon.size() < 3)
            return;
          if (polygon.size() == 3)
          {
            triangles.push_back({ polygon[0], polygon[1], polygon[2] });
            return;
          }
          if (polygon.size() == 4)
          {
            const auto d02 = edgePair(polygon[0], polygon[2]);
            const auto d13 = edgePair(polygon[1], polygon[3]);
            if (d02 <= d13)
            {
              triangles.push_back({ polygon[0], polygon[1], polygon[2] });
              triangles.push_back({ polygon[0], polygon[2], polygon[3] });
            }
            else
            {
              triangles.push_back({ polygon[0], polygon[1], polygon[3] });
              triangles.push_back({ polygon[1], polygon[2], polygon[3] });
            }
            return;
          }
          for (size_t i = 1; i + 1 < polygon.size(); ++i)
            triangles.push_back({ polygon[0], polygon[i], polygon[i + 1] });
        };

        auto orderCoplanarPolygon = [&](std::vector<Index>& polygon, const std::vector<Point>& verts)
        {
          if (polygon.size() < 4)
            return;

          Point center(3);
          center[0] = center[1] = center[2] = 0.0;
          for (const auto i : polygon)
            center += verts[i];
          center /= static_cast<Real>(polygon.size());

          const auto n =
            (verts[polygon[1]] - verts[polygon[0]]).cross(verts[polygon[2]] - verts[polygon[0]]);

          auto u = verts[polygon[0]] - center;
          if (u.norm() == 0.0)
            u = verts[polygon[1]] - center;
          u.normalize();

          auto v = n.cross(u);
          v.normalize();

          std::sort(
            polygon.begin(), polygon.end(),
            [&](Index i, Index j)
            {
              const auto pi = verts[i] - center;
              const auto pj = verts[j] - center;
              const Real ai = std::atan2(pi.dot(v), pi.dot(u));
              const Real aj = std::atan2(pj.dot(v), pj.dot(u));
              return ai < aj;
            });
        };

        // ---------------------------------------------------------------------
        // Copy original vertices (output will reuse all original vertices)
        // ---------------------------------------------------------------------
        std::vector<Point> vertices;
        vertices.reserve(mesh.getVertexCount());
        for (Index i = 0; i < static_cast<Index>(mesh.getVertexCount()); ++i)
          vertices.push_back(mesh.getVertexCoordinates(i));

        const Index originalVertexCount = static_cast<Index>(vertices.size());

        auto addVertex = [&](const Point& p)
        {
          vertices.push_back(p);
          return static_cast<Index>(vertices.size() - 1);
        };

        // ---------------------------------------------------------------------
        // Input lookup maps to RETAIN ORIGINAL FACE/EDGE ATTRIBUTES
        // keyed by original vertex sets.
        // ---------------------------------------------------------------------
        std::map<EdgeKey, Attribute> inEdgeAttrByVerts;
        std::map<FaceKey, Attribute> inFaceAttrByVerts;

        // Input edges
        for (auto eit = mesh.getPolytope(1); !eit.end(); ++eit)
        {
          const auto& ev = eit->getVertices(); // 2 vertices
          EdgeKey ek{ std::min(ev(0), ev(1)), std::max(ev(0), ev(1)) };
          inEdgeAttrByVerts[ek] = eit->getAttribute();
        }

        // Input triangular faces
        for (auto fit = mesh.getFace(); !fit.end(); ++fit)
        {
          const auto& fv = fit->getVertices(); // 3 vertices
          inFaceAttrByVerts[makeFaceKey(fv(0), fv(1), fv(2))] = fit->getAttribute();
        }

        // ---------------------------------------------------------------------
        // Geometry helpers
        // ---------------------------------------------------------------------
        auto signedVolume = [&](Index a, Index b, Index c, Index d)
        {
          const auto u = vertices[b] - vertices[a];
          const auto v = vertices[c] - vertices[a];
          const auto w = vertices[d] - vertices[a];
          return u.dot(v.cross(w));
        };

        auto makeTetrahedron =
          [&](Index a, Index b, Index c, Index d,
              Attribute attr, Attribute sourceAttr,
              typename Tet::Sign sign,
              std::vector<Tet>& tets)
        {
          Tet t{ { a, b, c, d }, attr, sourceAttr, sign };
          if (signedVolume(a, b, c, d) < 0)
            std::swap(t.v[1], t.v[2]);
          tets.push_back(t);
        };

        // ---------------------------------------------------------------------
        // Marching tetrahedra split
        // ---------------------------------------------------------------------
        std::map<EdgeKey, Index> edgeIntersections;
        std::vector<Tet> outTetrahedra;

        const std::array<std::array<int, 2>, 6> tetEdges{{
          {{0, 1}}, {{0, 2}}, {{0, 3}}, {{1, 2}}, {{1, 3}}, {{2, 3}}
        }};
        const std::array<std::array<int, 3>, 4> tetFaces{{
          {{0, 1, 2}}, {{0, 3, 1}}, {{0, 2, 3}}, {{1, 3, 2}}
        }};

        for (auto it = mesh.getCell(); !it.end(); ++it)
        {
          if (it->getGeometry() != Polytope::Type::Tetrahedron)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "MarchingTetrahedra only supports tetrahedral meshes."
              << Alert::Raise;
          }

          const auto& tv = it->getVertices();
          std::array<Index, 4> v{ tv(0), tv(1), tv(2), tv(3) };
          std::array<Real, 4>  phi{ ls[v[0]], ls[v[1]], ls[v[2]], ls[v[3]] };

          const auto cellAttr = it->getAttribute();

          size_t nneg = 0;
          for (size_t i = 0; i < 4; ++i)
            nneg += isNegative(phi[i]) ? 1 : 0;

          if (nneg == 0)
          {
            makeTetrahedron(
              v[0], v[1], v[2], v[3],
              getSplitLabel(cellAttr, Tet::Sign::Positive, m_cellSplit),
              cellAttr,
              Tet::Sign::Positive,
              outTetrahedra);
            continue;
          }

          if (nneg == 4)
          {
            makeTetrahedron(
              v[0], v[1], v[2], v[3],
              getSplitLabel(cellAttr, Tet::Sign::Negative, m_cellSplit),
              cellAttr,
              Tet::Sign::Negative,
              outTetrahedra);
            continue;
          }

          // -------------------------------------------------------------------
          // NEW (conforming "noSplit face/edge"):
          // If the interface would cut an INPUT face/edge whose attribute is
          // configured as NoSplit, then we must NOT split this cell, otherwise
          // we cannot retain original triangles/edges conformingly.
          //
          // This is implemented purely from the tetra's own vertices using
          // inEdgeAttrByVerts/inFaceAttrByVerts (no connectivity required).
          // -------------------------------------------------------------------
          bool blockedByNoSplit = false;

          // Check noSplit edges
          if (!blockedByNoSplit && !m_edgeSplit.empty())
          {
            for (const auto& e : tetEdges)
            {
              const int a = e[0];
              const int b = e[1];

              EdgeKey ek{ std::min(v[a], v[b]), std::max(v[a], v[b]) };
              auto itEdge = inEdgeAttrByVerts.find(ek);
              if (itEdge == inEdgeAttrByVerts.end())
                continue;

              if (!isNoSplit(itEdge->second, m_edgeSplit))
                continue;

              if (isNegative(phi[a]) != isNegative(phi[b]))
              {
                blockedByNoSplit = true;
                break;
              }
            }
          }

          // Check noSplit faces
          if (!blockedByNoSplit && !m_faceSplit.empty())
          {
            for (const auto& f : tetFaces)
            {
              const int a = f[0], b = f[1], c = f[2];
              auto itFace = inFaceAttrByVerts.find(makeFaceKey(v[a], v[b], v[c]));
              if (itFace == inFaceAttrByVerts.end())
                continue;

              if (!isNoSplit(itFace->second, m_faceSplit))
                continue;

              const bool sa = isNegative(phi[a]);
              const bool sb = isNegative(phi[b]);
              const bool sc = isNegative(phi[c]);

              if (!(sa == sb && sb == sc))
              {
                blockedByNoSplit = true;
                break;
              }
            }
          }

          // Your original cell-splitting gate (by cell attribute), plus the new blockers.
          if (blockedByNoSplit || !shouldSplitCell(cellAttr))
          {
            makeTetrahedron(v[0], v[1], v[2], v[3], cellAttr, cellAttr, Tet::Sign::Unknown, outTetrahedra);
            continue;
          }

          // -------------------------------------------------------------------
          // Normal split path (unchanged)
          // -------------------------------------------------------------------
          auto getIntersection = [&](int ia, int ib)
          {
            const auto va = v[ia];
            const auto vb = v[ib];
            const auto fa = phi[ia];
            const auto fb = phi[ib];

            if (fa == 0.0)
              return va;
            if (fb == 0.0)
              return vb;

            const EdgeKey ek{ std::min(va, vb), std::max(va, vb) };
            auto found = edgeIntersections.find(ek);
            if (found != edgeIntersections.end())
              return found->second;

            const Real t = fa / (fa - fb);
            const auto p = vertices[va] + t * (vertices[vb] - vertices[va]);
            const auto idx = addVertex(p);
            edgeIntersections.emplace(ek, idx);
            return idx;
          };

          std::vector<Index> interfacePolygon;
          for (const auto& e : tetEdges)
          {
            const int a = e[0];
            const int b = e[1];
            if (isNegative(phi[a]) != isNegative(phi[b]))
              pushUnique(interfacePolygon, getIntersection(a, b));
          }
          orderCoplanarPolygon(interfacePolygon, vertices);

          auto buildRegion = [&](bool negative)
          {
            const auto isInside = [&](Real x)
            {
              return negative ? isNegative(x) : !isNegative(x);
            };

            std::vector<Index> regionVertices;
            std::vector<std::array<Index, 3>> boundaryTriangles;

            for (size_t i = 0; i < 4; ++i)
            {
              if (isInside(phi[i]))
                pushUnique(regionVertices, v[i]);
            }
            for (const auto i : interfacePolygon)
              pushUnique(regionVertices, i);

            for (const auto& f : tetFaces)
            {
              std::vector<Index> polygon;
              for (int k = 0; k < 3; ++k)
              {
                const int ia = f[k];
                const int ib = f[(k + 1) % 3];
                if (isInside(phi[ia]))
                  polygon.push_back(v[ia]);
                if (isInside(phi[ia]) != isInside(phi[ib]))
                  polygon.push_back(getIntersection(ia, ib));
              }
              if (polygon.size() > 1 && polygon.front() == polygon.back())
                polygon.pop_back();
              triangulatePolygon(polygon, boundaryTriangles);
            }

            triangulatePolygon(interfacePolygon, boundaryTriangles);
            if (regionVertices.size() < 4 || boundaryTriangles.empty())
              return;

            Point centroid(3);
            centroid[0] = centroid[1] = centroid[2] = 0.0;
            for (const auto i : regionVertices)
              centroid += vertices[i];
            centroid /= static_cast<Real>(regionVertices.size());
            const auto c = addVertex(centroid);

            const auto sign = negative ? Tet::Sign::Negative : Tet::Sign::Positive;
            const auto outAttr = getSplitLabel(cellAttr, sign, m_cellSplit);

            for (const auto& t : boundaryTriangles)
              makeTetrahedron(c, t[0], t[1], t[2], outAttr, cellAttr, sign, outTetrahedra);
          };

          buildRegion(true);
          buildRegion(false);
        }

        // ---------------------------------------------------------------------
        // Build output mesh
        // ---------------------------------------------------------------------
        Mesh<Context::Local>::Builder builder;
        builder.initialize(3).nodes(vertices.size());
        for (const auto& p : vertices)
          builder.vertex(p);

        Index cell = 0;
        for (const auto& t : outTetrahedra)
        {
          builder.polytope(Polytope::Type::Tetrahedron, { t.v[0], t.v[1], t.v[2], t.v[3] });
          builder.attribute({ 3, cell++ }, t.attr);
        }

        auto out = builder.finalize();

        out.getConnectivity().compute(2, 3);
        out.getConnectivity().compute(1, 3);
        const auto& conn = out.getConnectivity();

        // ---------------------------------------------------------------------
        // Face attributes on output:
        // - If interface: interface attribute
        // - Else if the face is an ORIGINAL face (all vertices are original and the
        //   triangle exists in input): retain/transform using input face attribute
        // - Else fallback: previous behavior (based on incident cell sourceAttr)
        // ---------------------------------------------------------------------
        for (auto fit = out.getFace(); !fit.end(); ++fit)
        {
          const auto fidx = fit->getIndex();
          const auto& inc = conn.getIncidence({2, 3}, fidx); // face -> incident cells

          bool hasNeg = false;
          bool hasPos = false;
          std::optional<Attribute> fallbackCellAttr;

          for (const auto cidx : inc)
          {
            const auto& t = outTetrahedra[cidx];
            hasNeg = hasNeg || (t.sign == Tet::Sign::Negative);
            hasPos = hasPos || (t.sign == Tet::Sign::Positive);

            if (!fallbackCellAttr)
              fallbackCellAttr = t.sourceAttr;
            else if (*fallbackCellAttr != t.sourceAttr)
              fallbackCellAttr = std::nullopt;
          }

          if (hasNeg && hasPos)
          {
            out.setAttribute({2, fidx}, m_interfaceAttribute);
            continue;
          }

          const auto sign =
            hasNeg ? Tet::Sign::Negative :
            hasPos ? Tet::Sign::Positive : Tet::Sign::Unknown;

          const auto& fv = fit->getVertices(); // 3 vertices
          const bool allOriginal =
            (fv(0) < originalVertexCount) &&
            (fv(1) < originalVertexCount) &&
            (fv(2) < originalVertexCount);

          if (allOriginal)
          {
            auto itIn = inFaceAttrByVerts.find(makeFaceKey(fv(0), fv(1), fv(2)));
            if (itIn != inFaceAttrByVerts.end())
            {
              const Attribute baseFaceAttr = itIn->second; // input FACE attribute
              out.setAttribute({2, fidx}, getSplitLabel(baseFaceAttr, sign, m_faceSplit));
              continue;
            }
          }

          // Fallback (non-original face): keep old logic if you want a label anyway.
if (fallbackCellAttr)
  out.setAttribute({2, fidx}, m_fallbackFaceAttribute);
        }

        // ---------------------------------------------------------------------
        // Edge attributes on output:
        // - If interface: interface attribute
        // - Else if ORIGINAL edge (both vertices original and exists in input):
        //   retain/transform using input edge attribute
        // - Else fallback: previous behavior (based on incident cell sourceAttr)
        // ---------------------------------------------------------------------
        for (auto eit = out.getPolytope(1); !eit.end(); ++eit)
        {
          const auto eidx = eit->getIndex();
          const auto& inc = conn.getIncidence({1, 3}, eidx); // edge -> incident cells

          bool hasNeg = false;
          bool hasPos = false;
          std::optional<Attribute> fallbackCellAttr;

          for (const auto cidx : inc)
          {
            const auto& t = outTetrahedra[cidx];
            hasNeg = hasNeg || (t.sign == Tet::Sign::Negative);
            hasPos = hasPos || (t.sign == Tet::Sign::Positive);

            if (!fallbackCellAttr)
              fallbackCellAttr = t.sourceAttr;
            else if (*fallbackCellAttr != t.sourceAttr)
              fallbackCellAttr = std::nullopt;
          }

          if (hasNeg && hasPos)
          {
            out.setAttribute({1, eidx}, m_interfaceAttribute);
            continue;
          }

          const auto sign =
            hasNeg ? Tet::Sign::Negative :
            hasPos ? Tet::Sign::Positive : Tet::Sign::Unknown;

          const auto& ev = eit->getVertices(); // 2 vertices
          const bool allOriginal =
            (ev(0) < originalVertexCount) &&
            (ev(1) < originalVertexCount);

          if (allOriginal)
          {
            EdgeKey ek{ std::min(ev(0), ev(1)), std::max(ev(0), ev(1)) };
            auto itIn = inEdgeAttrByVerts.find(ek);
            if (itIn != inEdgeAttrByVerts.end())
            {
              const Attribute baseEdgeAttr = itIn->second; // input EDGE attribute
              out.setAttribute({1, eidx}, getSplitLabel(baseEdgeAttr, sign, m_edgeSplit));
              continue;
            }
          }

          // Fallback (non-original edge)
if (fallbackCellAttr)
  out.setAttribute({1, eidx}, m_fallbackEdgeAttribute);
        }

        return out;
      }

MarchingTetrahedra& setFallbackFaceAttribute(Attribute a)
{
  m_fallbackFaceAttribute = a;
  return *this;
}

MarchingTetrahedra& setFallbackEdgeAttribute(Attribute a)
{
  m_fallbackEdgeAttribute = a;
  return *this;
}

    private:
      Attribute m_negativeAttribute  = 1;
      Attribute m_positiveAttribute  = 2;
      Attribute m_interfaceAttribute = 10;

      SplitMap m_cellSplit;
      SplitMap m_faceSplit;
      SplitMap m_edgeSplit;

Attribute m_fallbackFaceAttribute = 0;
Attribute m_fallbackEdgeAttribute = 0;
  };
}

#endif
