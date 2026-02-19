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
#include <tuple>
#include <set>
#include <vector>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Variational/GridFunction.h"

#include "Mesh.h"

namespace Rodin::Geometry
{
  class MarchingTetrahedra
  {
    public:
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

        struct Tet
        {
          std::array<Index, 4> v;
          Attribute attr;
        };

        std::vector<Point> vertices;
        vertices.reserve(mesh.getVertexCount());
        for (Index i = 0; i < static_cast<Index>(mesh.getVertexCount()); ++i)
          vertices.push_back(mesh.getVertexCoordinates(i));

        auto addVertex = [&](const Point& p)
        {
          vertices.push_back(p);
          return static_cast<Index>(vertices.size() - 1);
        };

        auto signedVolume = [&](Index a, Index b, Index c, Index d)
        {
          const auto u = vertices[b] - vertices[a];
          const auto v = vertices[c] - vertices[a];
          const auto w = vertices[d] - vertices[a];
          return u.dot(v.cross(w));
        };

        auto makeTetrahedron = [&](Index a, Index b, Index c, Index d, Attribute attr, std::vector<Tet>& tets)
        {
          Tet t{ { a, b, c, d }, attr };
          if (signedVolume(a, b, c, d) < 0)
            std::swap(t.v[1], t.v[2]);
          tets.push_back(t);
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

        auto orderCoplanarPolygon = [&](std::vector<Index>& polygon)
        {
          if (polygon.size() < 4)
            return;

          Point center(3);
          center[0] = center[1] = center[2] = 0.0;
          for (const auto i : polygon)
            center += vertices[i];
          center /= static_cast<Real>(polygon.size());

          const auto n = (vertices[polygon[1]] - vertices[polygon[0]])
                       .cross(vertices[polygon[2]] - vertices[polygon[0]]);
          auto u = vertices[polygon[0]] - center;
          if (u.norm() == 0.0)
            u = vertices[polygon[1]] - center;
          u.normalize();
          auto v = n.cross(u);
          v.normalize();

          std::sort(
              polygon.begin(), polygon.end(),
              [&](Index i, Index j)
              {
                const auto pi = vertices[i] - center;
                const auto pj = vertices[j] - center;
                const Real ai = std::atan2(pi.dot(v), pi.dot(u));
                const Real aj = std::atan2(pj.dot(v), pj.dot(u));
                return ai < aj;
              });
        };

        std::map<EdgeKey, Index> edgeIntersections;
        std::vector<Tet> outTetrahedra;
        const std::array<std::array<int, 2>, 6> tetEdges{{
          {{0, 1}}, {{0, 2}}, {{0, 3}}, {{1, 2}}, {{1, 3}}, {{2, 3}}
        }};
        const std::array<std::array<int, 3>, 4> tetFaces{{
          {{0, 1, 2}}, {{0, 3, 1}}, {{0, 2, 3}}, {{1, 3, 2}}
        }};

        auto isNegative = [](Real x) { return x < 0.0; };

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
          std::array<Real, 4> phi{ ls[v[0]], ls[v[1]], ls[v[2]], ls[v[3]] };
          size_t nneg = 0;
          for (size_t i = 0; i < 4; ++i)
            nneg += isNegative(phi[i]) ? 1 : 0;

          if (nneg == 0)
          {
            makeTetrahedron(v[0], v[1], v[2], v[3], m_positiveAttribute, outTetrahedra);
            continue;
          }
          if (nneg == 4)
          {
            makeTetrahedron(v[0], v[1], v[2], v[3], m_negativeAttribute, outTetrahedra);
            continue;
          }

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
          orderCoplanarPolygon(interfacePolygon);

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
            if (regionVertices.size() < 4 || boundaryTriangles.size() == 0)
              return;

            Point centroid(3);
            centroid[0] = centroid[1] = centroid[2] = 0.0;
            for (const auto i : regionVertices)
              centroid += vertices[i];
            centroid /= static_cast<Real>(regionVertices.size());
            const auto c = addVertex(centroid);

            const auto attr = negative ? m_negativeAttribute : m_positiveAttribute;
            for (const auto& t : boundaryTriangles)
              makeTetrahedron(c, t[0], t[1], t[2], attr, outTetrahedra);
          };

          buildRegion(true);
          buildRegion(false);
        }

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
        out.trace(Map<std::pair<Attribute, Attribute>, Attribute>{
          { {m_negativeAttribute, m_positiveAttribute}, m_interfaceAttribute },
          { {m_positiveAttribute, m_negativeAttribute}, m_interfaceAttribute }
        });
        return out;
      }

    private:
      Attribute m_negativeAttribute = 1;
      Attribute m_positiveAttribute = 2;
      Attribute m_interfaceAttribute = 10;
  };
}

#endif
