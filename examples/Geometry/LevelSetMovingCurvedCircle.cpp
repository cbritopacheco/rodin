/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Adaptation/TargetMatrixOptimization/Metrics.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include <Eigen/SparseLU>

#include <Rodin/Adaptation.h>
#include <Rodin/Assembly/Default.h>
#include <Rodin/Geometry.h>
#include <Rodin/Geometry/LevelSetDiscretizerTriangles.h>
#include <Rodin/Geometry/TriangleMeshOptimizer.h>
#include <Rodin/IO.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Adaptation::TargetMatrixOptimization;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  constexpr Real Pi = 3.14159265358979323846;
  constexpr Attribute Interface = 30;
  constexpr Attribute Boundary = 40;
  constexpr Attribute Negative = 1;
  constexpr Attribute Positive = 2;

  Math::SpatialPoint center(Real t)
  {
    return Math::SpatialPoint{
      Real(0.5) + Real(0.18) * std::sin(Real(2) * Pi * t),
      Real(0.5) + Real(0.13) * std::sin(Real(4) * Pi * t + Real(0.35))
    };
  }

  Real wavyRadius(Real theta, Real t)
  {
    return Real(0.17)
      + Real(0.035) * std::sin(
          Real(5) * theta + Real(0.65) * std::sin(Real(2) * Pi * t));
  }

  Real wavyRadiusDerivative(Real theta, Real t)
  {
    return Real(0.175) * std::cos(
        Real(5) * theta + Real(0.65) * std::sin(Real(2) * Pi * t));
  }

  Real phiAt(const Math::SpatialPoint& x, Real t)
  {
    const auto c = center(t);
    const Real dx = x[0] - c[0];
    const Real dy = x[1] - c[1];
    const Real theta = std::atan2(dy, dx);
    return std::hypot(dx, dy) - wavyRadius(theta, t);
  }

  Math::SpatialPoint gradPhiAt(const Math::SpatialPoint& x, Real t)
  {
    const auto c = center(t);
    const Real dx = x[0] - c[0];
    const Real dy = x[1] - c[1];
    const Real rho2 = dx * dx + dy * dy;
    const Real rho = std::sqrt(rho2);
    if (rho <= Real(1e-14))
      return Math::SpatialPoint{1, 0};
    const Real theta = std::atan2(dy, dx);
    const Real dr = wavyRadiusDerivative(theta, t);
    return Math::SpatialPoint{
      dx / rho + dr * dy / rho2,
      dy / rho - dr * dx / rho2
    };
  }

  Real boxBoundaryValue(const Math::SpatialPoint& x)
  {
    const std::array<Real, 4> d = {{
      std::abs(x[0]), std::abs(x[0] - Real(1)),
      std::abs(x[1]), std::abs(x[1] - Real(1))
    }};
    const auto side =
      static_cast<size_t>(std::min_element(d.begin(), d.end()) - d.begin());
    if (side == 0)
      return x[0];
    if (side == 1)
      return x[0] - Real(1);
    if (side == 2)
      return x[1];
    return x[1] - Real(1);
  }

  Math::SpatialPoint boxBoundaryGradient(const Math::SpatialPoint& x)
  {
    const std::array<Real, 4> d = {{
      std::abs(x[0]), std::abs(x[0] - Real(1)),
      std::abs(x[1]), std::abs(x[1] - Real(1))
    }};
    const auto side =
      static_cast<size_t>(std::min_element(d.begin(), d.end()) - d.begin());
    if (side < 2)
      return Math::SpatialPoint{1, 0};
    return Math::SpatialPoint{0, 1};
  }

  Math::SpatialPoint projectToInterface(
      const Math::SpatialPoint& x,
      Real t)
  {
    Math::SpatialPoint p = x;
    for (int i = 0; i < 2; ++i)
    {
      const Real f = phiAt(p, t);
      const auto g = gradPhiAt(p, t);
      const Real gg = g[0] * g[0] + g[1] * g[1];
      if (gg < Real(1e-30))
        break;
      p = Math::SpatialPoint{ p[0] - f * g[0] / gg,
                              p[1] - f * g[1] / gg };
    }
    return p;
  }

  void annotateBoundary(LocalMesh& mesh)
  {
    const auto& conn = mesh.getConnectivity();
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto& edge = conn.getPolytope(1, e);
      const auto x0 = mesh.getVertexCoordinates(edge(0));
      const auto x1 = mesh.getVertexCoordinates(edge(1));
      const Real eps = Real(1e-12);
      const bool on0 =
        std::abs(x0[0]) <= eps || std::abs(x0[0] - Real(1)) <= eps
     || std::abs(x0[1]) <= eps || std::abs(x0[1] - Real(1)) <= eps;
      const bool on1 =
        std::abs(x1[0]) <= eps || std::abs(x1[0] - Real(1)) <= eps
     || std::abs(x1[1]) <= eps || std::abs(x1[1] - Real(1)) <= eps;
      if (on0 && on1)
        mesh.setAttribute({ 1, e }, Boundary);
    }
  }

  std::vector<char> interfaceVertexMask(const LocalMesh& mesh)
  {
    std::vector<char> mask(mesh.getVertexCount(), 0);
    const auto& conn = mesh.getConnectivity();
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto attr = mesh.getAttribute(1, e);
      if (!attr || *attr != Interface)
        continue;
      const auto& edge = conn.getPolytope(1, e);
      mask[edge(0)] = 1;
      mask[edge(1)] = 1;
    }
    return mask;
  }

  Math::SpatialPoint triangleReferenceVertex(std::uint8_t local)
  {
    if (local == 1)
      return Math::SpatialPoint{1, 0};
    if (local == 2)
      return Math::SpatialPoint{0, 1};
    return Math::SpatialPoint{0, 0};
  }

  Math::SpatialPoint triangleReferenceEdgePoint(
      std::uint8_t localA,
      std::uint8_t localB,
      Real s)
  {
    return (Real(1) - s) * triangleReferenceVertex(localA)
         + s * triangleReferenceVertex(localB);
  }

  bool referencePointOnTriangleEdge(
      const Math::SpatialPoint& rc,
      std::uint8_t localA,
      std::uint8_t localB)
  {
    const Real eps = Real(1e-12);
    if ((localA == 0 && localB == 1) || (localA == 1 && localB == 0))
      return std::abs(rc[1]) <= eps;
    if ((localA == 1 && localB == 2) || (localA == 2 && localB == 1))
      return std::abs(rc[0] + rc[1] - Real(1)) <= eps;
    return std::abs(rc[0]) <= eps;
  }

  bool referencePointIsTriangleVertex(const Math::SpatialPoint& rc)
  {
    const Real eps = Real(1e-12);
    return (std::abs(rc[0]) <= eps && std::abs(rc[1]) <= eps)
        || (std::abs(rc[0] - Real(1)) <= eps && std::abs(rc[1]) <= eps)
        || (std::abs(rc[0]) <= eps && std::abs(rc[1] - Real(1)) <= eps);
  }

  Optional<std::array<std::uint8_t, 2>> localEdgeForMeshEdge(
      const Geometry::Polytope::Key& cell,
      const Geometry::Polytope::Key& edge)
  {
    for (std::uint8_t a = 0; a < 3; ++a)
      for (std::uint8_t b = a + 1; b < 3; ++b)
        if ((cell(a) == edge(0) && cell(b) == edge(1))
         || (cell(a) == edge(1) && cell(b) == edge(0)))
          return std::array<std::uint8_t, 2>{{a, b}};
    return {};
  }

  std::vector<std::array<std::uint8_t, 2>> interfaceLocalEdges(
      const LocalMesh& mesh,
      Index cellIndex)
  {
    std::vector<std::array<std::uint8_t, 2>> edges;
    const auto& conn = mesh.getConnectivity();
    if (conn.getGeometry(2, cellIndex) != Polytope::Type::Triangle)
      return edges;
    const auto& cell = conn.getPolytope(2, cellIndex);
    for (Index edgeIndex : conn.getIncidence({2, 1}, cellIndex))
    {
      const auto attr = mesh.getAttribute(1, edgeIndex);
      if (!attr || *attr != Interface)
        continue;
      if (const auto edge =
            localEdgeForMeshEdge(cell, conn.getPolytope(1, edgeIndex)))
        edges.push_back(*edge);
    }
    return edges;
  }

  void installQuadraticGeometry(LocalMesh& mesh, Real t)
  {
    if (mesh.getConnectivity().getIncidence(2, 1).size() == 0)
      return;

    Variational::RealH1Element<2> fe(Polytope::Type::Triangle);
    Variational::RealH1Element<2> edgeFe(Polytope::Type::Segment);
    const auto& conn = mesh.getConnectivity();
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto attr = mesh.getAttribute(1, e);
      const bool isInterface = attr && *attr == Interface;
      const auto& edge = conn.getPolytope(1, e);
      const auto x0 = mesh.getVertexCoordinates(edge(0));
      const auto x1 = mesh.getVertexCoordinates(edge(1));
      PointCloud pm(mesh.getSpaceDimension(), edgeFe.getCount());
      for (size_t local = 0; local < edgeFe.getCount(); ++local)
      {
        const Real s = edgeFe.getNode(local)[0];
        Math::SpatialPoint x = (Real(1) - s) * x0 + s * x1;
        if (isInterface
            && s > Real(1e-12)
            && s < Real(1) - Real(1e-12))
          x = projectToInterface(x, t);
        for (size_t d = 0; d < static_cast<size_t>(x.size()); ++d)
          pm(d, local) = x[d];
      }
      mesh.setPolytopeTransformation(
          {1, e},
          new ParametricTransformation<Variational::RealH1Element<2>>(
              pm, Variational::RealH1Element<2>(edgeFe)));
    }

    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      if (conn.getGeometry(2, c) != Polytope::Type::Triangle)
        continue;
      const auto featureEdges = interfaceLocalEdges(mesh, c);
      const auto& cell = conn.getPolytope(2, c);
      const auto x0 = mesh.getVertexCoordinates(cell(0));
      const auto x1 = mesh.getVertexCoordinates(cell(1));
      const auto x2 = mesh.getVertexCoordinates(cell(2));
      PointCloud pm(mesh.getSpaceDimension(), fe.getCount());
      for (size_t local = 0; local < fe.getCount(); ++local)
      {
        const auto& rc = fe.getNode(local);
        Math::SpatialPoint x =
          (Real(1) - rc[0] - rc[1]) * x0 + rc[0] * x1 + rc[1] * x2;
        bool onInterfaceEdge = false;
        for (const auto& edge : featureEdges)
          onInterfaceEdge = onInterfaceEdge
            || referencePointOnTriangleEdge(rc, edge[0], edge[1]);
        if (!referencePointIsTriangleVertex(rc) && onInterfaceEdge)
          x = projectToInterface(x, t);
        for (size_t d = 0; d < static_cast<size_t>(x.size()); ++d)
          pm(d, local) = x[d];
      }
      auto* trans =
        new ParametricTransformation<Variational::RealH1Element<2>>(
            pm, Variational::RealH1Element<2>(fe));
      bool valid = true;
      for (const Math::SpatialPoint& rc : {
            Math::SpatialPoint{Real(1) / Real(3), Real(1) / Real(3)},
            Math::SpatialPoint{Real(0.5), Real(0.25)},
            Math::SpatialPoint{Real(0.25), Real(0.5)},
            Math::SpatialPoint{Real(0.25), Real(0.25)}})
      {
        Math::SpatialMatrix<Real> J;
        trans->jacobian(J, rc);
        valid = valid && J.rows() == 2 && J.cols() == 2
          && std::isfinite(J.determinant()) && J.determinant() > Real(1e-14);
      }
      if (valid)
        mesh.setPolytopeTransformation({2, c}, trans);
      else
        delete trans;
    }
  }

  std::pair<Real, Real> curvedInterfaceFit(const LocalMesh& mesh, Real t)
  {
    Real sq = 0;
    Real max = 0;
    Index count = 0;
    const auto& conn = mesh.getConnectivity();
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      if (conn.getGeometry(2, c) != Polytope::Type::Triangle)
        continue;
      const auto cellIt = mesh.getPolytope(2, c);
      const auto& trans = cellIt->getTransformation();
      for (const auto& edge : interfaceLocalEdges(mesh, c))
      {
        for (Real s : {Real(0.25), Real(0.5), Real(0.75)})
        {
          Math::SpatialPoint x;
          trans.transform(x, triangleReferenceEdgePoint(edge[0], edge[1], s));
          const Real phi = std::abs(phiAt(x, t));
          max = std::max(max, phi);
          sq += phi * phi;
          ++count;
        }
      }
    }
    return {
      count > 0 ? std::sqrt(sq / static_cast<Real>(count)) : Real(0),
      max
    };
  }

  std::pair<Real, Index> curvedJacobianStats(const LocalMesh& mesh)
  {
    Real minJacobian = std::numeric_limits<Real>::infinity();
    Index invalid = 0;
    const auto& conn = mesh.getConnectivity();
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      if (conn.getGeometry(2, c) != Polytope::Type::Triangle)
        continue;
      const auto cellIt = mesh.getPolytope(2, c);
      const auto& trans = cellIt->getTransformation();
      for (const Math::SpatialPoint& rc : {
            Math::SpatialPoint{Real(1) / Real(3), Real(1) / Real(3)},
            Math::SpatialPoint{Real(0.5), Real(0.25)},
            Math::SpatialPoint{Real(0.25), Real(0.5)},
            Math::SpatialPoint{Real(0.25), Real(0.25)}})
      {
        Math::SpatialMatrix<Real> J;
        trans.jacobian(J, rc);
        const bool shaped = J.rows() == 2 && J.cols() == 2;
        const Real det = shaped ? J.determinant() : Real(0);
        if (!shaped || !(det > Real(0)) || !std::isfinite(det))
          ++invalid;
        else
          minJacobian = std::min(minJacobian, det);
      }
    }
    if (!std::isfinite(minJacobian))
      minJacobian = 0;
    return {minJacobian, invalid};
  }

  LocalMesh curvedInterfacePolyline(
      const LocalMesh& mesh,
      Index samplesPerEdge = 4)
  {
    std::vector<Math::SpatialPoint> points;
    std::vector<std::array<Index, 2>> segments;
    const auto& conn = mesh.getConnectivity();
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      if (conn.getGeometry(2, c) != Polytope::Type::Triangle)
        continue;
      const auto cellIt = mesh.getPolytope(2, c);
      const auto& trans = cellIt->getTransformation();
      for (const auto& edge : interfaceLocalEdges(mesh, c))
      {
        const Index first = static_cast<Index>(points.size());
        for (Index i = 0; i <= samplesPerEdge; ++i)
        {
          Math::SpatialPoint x;
          trans.transform(
              x,
              triangleReferenceEdgePoint(
                edge[0], edge[1],
                static_cast<Real>(i) / static_cast<Real>(samplesPerEdge)));
          points.push_back(x);
          if (i > 0)
            segments.push_back({{first + i - 1, first + i}});
        }
      }
    }

    auto builder = LocalMesh::Builder();
    builder
      .initialize(2)
      .nodes(points.size())
      .reserve(1, segments.size());
    for (const auto& x : points)
      builder.vertex(x);
    for (const auto& segment : segments)
      builder.polytope(Polytope::Type::Segment, {segment[0], segment[1]});
    auto out = builder.finalize();
    out.getConnectivity().compute(1, 0);
    for (Index e = 0; e < out.getPolytopeCount(1); ++e)
      out.setAttribute({1, e}, Interface);
    return out;
  }

  LocalMesh curvedVisualizationMesh(
      const LocalMesh& mesh,
      Index subdivisions = 3)
  {
    std::vector<Math::SpatialPoint> points;
    std::vector<std::array<Index, 3>> cells;
    std::vector<Attribute> cellAttributes;
    const auto& conn = mesh.getConnectivity();
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      if (conn.getGeometry(2, c) != Polytope::Type::Triangle)
        continue;
      const auto cellIt = mesh.getPolytope(2, c);
      const auto& trans = cellIt->getTransformation();
      const Index base = static_cast<Index>(points.size());
      for (Index i = 0; i <= subdivisions; ++i)
      {
        for (Index j = 0; j <= subdivisions - i; ++j)
        {
          Math::SpatialPoint x;
          trans.transform(
              x,
              Math::SpatialPoint{
                static_cast<Real>(i) / static_cast<Real>(subdivisions),
                static_cast<Real>(j) / static_cast<Real>(subdivisions)});
          points.push_back(x);
        }
      }

      auto idx = [base, subdivisions](Index i, Index j)
      {
        Index offset = 0;
        for (Index r = 0; r < i; ++r)
          offset += subdivisions - r + 1;
        return base + offset + j;
      };

      for (Index i = 0; i < subdivisions; ++i)
      {
        for (Index j = 0; j < subdivisions - i; ++j)
        {
          cells.push_back({{idx(i, j), idx(i + 1, j), idx(i, j + 1)}});
          cellAttributes.push_back(
              mesh.getAttribute(2, c).value_or(Attribute(0)));
          if (j + i + 1 < subdivisions)
          {
            cells.push_back({{idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)}});
            cellAttributes.push_back(
                mesh.getAttribute(2, c).value_or(Attribute(0)));
          }
        }
      }
    }

    auto builder = LocalMesh::Builder();
    builder
      .initialize(2)
      .nodes(points.size())
      .reserve(2, cells.size());
    for (const auto& x : points)
      builder.vertex(x);
    for (const auto& cell : cells)
      builder.polytope(Polytope::Type::Triangle, {cell[0], cell[1], cell[2]});
    auto out = builder.finalize();
    out.getConnectivity().compute(2, 1);
    out.getConnectivity().compute(1, 0);
    for (Index c = 0;
         c < static_cast<Index>(std::min(cells.size(), cellAttributes.size()));
         ++c)
    {
      if (cellAttributes[static_cast<size_t>(c)] != Attribute(0))
        out.setAttribute({2, c}, cellAttributes[static_cast<size_t>(c)]);
    }
    return out;
  }

  template <class GridFunctionType>
  void applyDisplacement(LocalMesh& mesh, const GridFunctionType& u)
  {
    const auto& data = u.getData();
    const Index nv = mesh.getVertexCount();
    for (Index v = 0; v < nv; ++v)
    {
      auto x = mesh.getVertexCoordinates(v);
      x[0] += data(v);
      x[1] += data(v + nv);
      mesh.setVertexCoordinates(v, x);
    }
  }

  template <class FES, class Data>
  void applyParametricDisplacement(
      LocalMesh& mesh,
      const GridFunction<FES, Data>& u)
  {
    const auto& fes = u.getFiniteElementSpace();
    const auto& data = u.getData();
    const size_t sdim = mesh.getSpaceDimension();
    const auto& conn = mesh.getConnectivity();
    std::vector<PointCloud> edgePointClouds(conn.getCount(1));
    std::vector<PointCloud> cellPointClouds(mesh.getCellCount());
    std::vector<char> hasCellPointCloud(mesh.getCellCount(), 0);

    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto edgeIt = mesh.getPolytope(1, e);
      const auto& edge = *edgeIt;
      const auto& fe = fes.getFiniteElement(edge.getDimension(), edge.getIndex());
      PointCloud pm(sdim, fe.getCount() / sdim);
      for (size_t node = 0; node < fe.getCount() / sdim; ++node)
      {
        const auto& rc = fe.getNode(node * sdim);
        const Geometry::Point point(edge, rc);
        auto x = point.getPhysicalCoordinates();
        for (size_t c = 0; c < sdim; ++c)
        {
          const size_t local = node * sdim + c;
          x[c] += data(fes.getGlobalIndex(
              {edge.getDimension(), edge.getIndex()},
              static_cast<Index>(local)));
        }
        for (size_t c = 0; c < sdim; ++c)
          pm(c, node) = x[c];
      }
      edgePointClouds[static_cast<size_t>(e)] = std::move(pm);
    }

    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      if (conn.getGeometry(2, c) != Polytope::Type::Triangle)
        continue;
      const auto cellIt = mesh.getPolytope(2, c);
      const auto& cell = *cellIt;
      const auto& fe = fes.getFiniteElement(cell.getDimension(), cell.getIndex());
      PointCloud pm(sdim, fe.getCount() / sdim);
      for (size_t node = 0; node < fe.getCount() / sdim; ++node)
      {
        const auto& rc = fe.getNode(node * sdim);
        const Geometry::Point point(cell, rc);
        auto x = point.getPhysicalCoordinates();
        for (size_t d = 0; d < sdim; ++d)
        {
          const size_t local = node * sdim + d;
          x[d] += data(fes.getGlobalIndex(
              {cell.getDimension(), cell.getIndex()},
              static_cast<Index>(local)));
        }
        for (size_t d = 0; d < sdim; ++d)
          pm(d, node) = x[d];
      }
      cellPointClouds[static_cast<size_t>(c)] = std::move(pm);
      hasCellPointCloud[static_cast<size_t>(c)] = 1;
    }

    for (Index v = 0; v < static_cast<Index>(mesh.getVertexCount()); ++v)
    {
      auto x = mesh.getVertexCoordinates(v);
      for (size_t d = 0; d < sdim; ++d)
        x[d] += data(fes.getGlobalIndex({0, v}, static_cast<Index>(d)));
      mesh.setVertexCoordinates(v, x);
    }

    // setVertexCoordinates() flushes cached/custom transformations. Reinstall
    // the displaced P2 geometry after the vertex update so midside/interface
    // nodes remain part of the actual mesh map used for diagnostics/export.
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      mesh.setPolytopeTransformation(
          {1, e},
          new ParametricTransformation<Variational::RealH1Element<2>>(
              edgePointClouds[static_cast<size_t>(e)],
              Variational::RealH1Element<2>(Polytope::Type::Segment)));
    }
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      if (!hasCellPointCloud[static_cast<size_t>(c)])
        continue;
      mesh.setPolytopeTransformation(
          {2, c},
          new ParametricTransformation<Variational::RealH1Element<2>>(
              cellPointClouds[static_cast<size_t>(c)],
              Variational::RealH1Element<2>(Polytope::Type::Triangle)));
    }
  }

  std::vector<Math::SpatialPoint> vertexCoordinates(const LocalMesh& mesh)
  {
    std::vector<Math::SpatialPoint> coordinates;
    coordinates.reserve(static_cast<size_t>(mesh.getVertexCount()));
    for (Index v = 0; v < mesh.getVertexCount(); ++v)
      coordinates.push_back(mesh.getVertexCoordinates(v));
    return coordinates;
  }

  void restoreVertexCoordinates(
      LocalMesh& mesh,
      const std::vector<Math::SpatialPoint>& coordinates)
  {
    const Index count = std::min(
        mesh.getVertexCount(),
        static_cast<Index>(coordinates.size()));
    for (Index v = 0; v < count; ++v)
      mesh.setVertexCoordinates(v, coordinates[static_cast<size_t>(v)]);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
  }

  void demoteInterface(LocalMesh& mesh)
  {
    auto& conn = mesh.getConnectivity();
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto attr = mesh.getAttribute(1, e);
      if (attr && *attr == Interface)
        mesh.setAttribute({ 1, e }, Attribute(0));
    }
  }

  std::pair<Real, Index> meshQuality(const LocalMesh& mesh)
  {
    const auto& conn = mesh.getConnectivity();
    Real qmin = std::numeric_limits<Real>::infinity();
    Index inverted = 0;
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      const auto& cell = conn.getPolytope(2, c);
      const auto x0 = mesh.getVertexCoordinates(cell(0));
      const auto x1 = mesh.getVertexCoordinates(cell(1));
      const auto x2 = mesh.getVertexCoordinates(cell(2));
      const Real orient =
        (x1[0] - x0[0]) * (x2[1] - x0[1])
      - (x1[1] - x0[1]) * (x2[0] - x0[0]);
      if (orient <= 0)
        ++inverted;
      const Real area = Real(0.5) * std::abs(orient);
      const Real den = (x1 - x0).squaredNorm()
                     + (x2 - x1).squaredNorm()
                     + (x0 - x2).squaredNorm();
      if (den > 0)
        qmin = std::min(qmin, Real(4) * std::sqrt(Real(3)) * area / den);
    }
    if (!std::isfinite(qmin))
      qmin = 0;
    return {qmin, inverted};
  }
}

int main(int, char**)
{
  constexpr size_t resolution = 10;
  constexpr Index steps = 40;
  constexpr Real h = Real(1) / static_cast<Real>(resolution - 1);

  LocalMesh background =
    LocalMesh::UniformGrid(Polytope::Type::Triangle, { resolution, resolution });
  background.scale(h);

  IO::XDMF xdmf("LevelSetMovingCurvedCircle");

  // Diagnostics are accumulated and reported as an average over the whole
  // process, not just the final step.
  Real sumQmin = 0, sumFit = 0;
  Index sumInverted = 0, rejectedTMOP = 0, diagSteps = 0;

  for (Index step = 0; step < steps; ++step)
  {
    const Real t = static_cast<Real>(step)
      / static_cast<Real>(std::max<Index>(1, steps - 1));

    background.getConnectivity().compute(2, 1);
    background.getConnectivity().compute(1, 0);

    P1<Real, LocalMesh> phiSpace(background);
    GridFunction phi(phiSpace);
    for (Index i = 0; i < background.getVertexCount(); ++i)
      phi[i] = phiAt(background.getVertexCoordinates(i), t);

    auto cut = LevelSetDiscretizerTriangles(phi)
      .setSignTolerance(1e-12)
      .setInterfaceAttribute(Interface)
      .setNegativeCellAttribute(Negative)
      .setPositiveCellAttribute(Positive)
      .setCrossingSnapTolerance(0.05)
      .setMinCutQuality(0.05)
      .discretize();
    annotateBoundary(cut.mesh);

    LocalMesh optimized = cut.mesh;
    {
      const auto cutInterfaceMask = interfaceVertexMask(optimized);
      const auto parameters =
        TriangleMeshOptimizerParameters::levelSetCarryForward(h);
      TriangleMeshOptimizer()
        .setParameters(parameters)
        .setFeatureProjection([t](const Math::SpatialPoint& x)
          { return projectToInterface(x, t); })
        .setProtectedVertices(cutInterfaceMask)
        .optimize(optimized);
    }
    optimized.getConnectivity().compute(2, 1);
    optimized.getConnectivity().compute(1, 0);
    installQuadraticGeometry(optimized, t);

    const auto featureMask = interfaceVertexMask(optimized);
    CurvedQualityTargetJacobian target(optimized, 0.75);

    VectorH1<2, LocalMesh> space(std::integral_constant<size_t, 2>{}, optimized, 2);
    optimized.getConnectivity().compute(1, 2);
    GridFunction u(space);
    u.getData().setZero();
    TrialFunction du(space);
    TestFunction v(space);

    auto phiValue = [t](const Math::SpatialPoint& x) { return phiAt(x, t); };
    auto phiGradient =
      [t](const Math::SpatialPoint& x) { return gradPhiAt(x, t); };

    QualityTerm quality(SquaredDistanceMetric{}, target, 0.04);
    quality.setQuadratureOrder(4);
    DeviationTerm deviation(0.25);
    AnalyticLevelSetFitTerm fit(
        phiValue, phiGradient, Optional<Attribute>(Interface), 2.0);
    AnalyticLevelSetFitTerm boundaryFit(
        boxBoundaryValue, boxBoundaryGradient,
        Optional<Attribute>(Boundary), 1.0);

    auto makeTangent = [&]()
    {
      return quality.tangent(u, du, v)
           + deviation.tangent(du, v)
           + fit.tangent(u, du, v)
           + boundaryFit.tangent(u, du, v);
    };
    auto makeResidual = [&]()
    {
      return quality.residual(u, v)
           + deviation.residual(u, v)
           + fit.residual(u, v)
           + boundaryFit.residual(u, v);
    };
    auto meritEnergy = [&]()
    {
      return quality.energy(u)
           + deviation.energy(u)
           + fit.energy(u)
           + boundaryFit.energy(u);
    };

    bool solvedTMOP = true;
    Real lastAlpha = 1;
    Index tmopIterations = 0;
    try
    {
      for (Index it = 0; it < 5; ++it)
      {
        ++tmopIterations;
        LinearForm R(v);
        R = makeResidual();
        R.assemble();
        const auto r = R.getVector();
        if (!r.allFinite())
        {
          solvedTMOP = false;
          break;
        }
        if (r.norm() <= Real(1e-10))
          break;

        BilinearForm J(du, v);
        J = makeTangent();
        J.assemble();

        Eigen::SparseLU<std::decay_t<decltype(J.getOperator())>> lu;
        lu.compute(J.getOperator());
        if (lu.info() != Eigen::Success)
        {
          solvedTMOP = false;
          break;
        }
        const Math::Vector<Real> dx = lu.solve(-r);
        if (lu.info() != Eigen::Success || !dx.allFinite())
        {
          solvedTMOP = false;
          break;
        }

        const Real e0 = meritEnergy();
        const Math::Vector<Real> u0 = u.getData();
        Real alpha = 1;
        bool acceptedStep = false;
        for (int ls = 0; ls < 30; ++ls)
        {
          u.getData() = u0 + alpha * dx;
          const Real e = meritEnergy();
          if (std::isfinite(e) && e <= e0 * (Real(1) + Real(1e-12)))
          {
            acceptedStep = true;
            break;
          }
          alpha *= Real(0.5);
        }
        lastAlpha = alpha;
        if (!acceptedStep)
        {
          u.getData() = u0;
          break;
        }
        if (alpha * dx.norm() <= Real(1e-10))
          break;
      }
    }
    catch (...)
    {
      solvedTMOP = false;
    }

    LocalMesh beforeTMOP(optimized);
    if (solvedTMOP && u.getData().allFinite())
      applyParametricDisplacement(optimized, u);
    else
      solvedTMOP = false;

    optimized.getConnectivity().compute(2, 1);
    optimized.getConnectivity().compute(1, 0);
    const auto [tmopQmin, tmopInverted] = meshQuality(optimized);
    auto [tmopMinJacobian, tmopInvalidJacobianSamples] =
      curvedJacobianStats(optimized);
    const bool acceptedTMOP =
      solvedTMOP
      && tmopInverted == 0
      && std::isfinite(tmopQmin)
      && tmopQmin > Real(0)
      && tmopInvalidJacobianSamples == 0
      && tmopMinJacobian > Real(0)
      && std::isfinite(tmopMinJacobian);
    if (!acceptedTMOP)
    {
      optimized = std::move(beforeTMOP);
      ++rejectedTMOP;
      const auto fallbackJacobianStats = curvedJacobianStats(optimized);
      tmopMinJacobian = fallbackJacobianStats.first;
      tmopInvalidJacobianSamples = fallbackJacobianStats.second;
    }
    optimized.getConnectivity().compute(2, 1);
    optimized.getConnectivity().compute(1, 0);

    P1<Real, LocalMesh> outputSpace(optimized);
    GridFunction outputPhi(outputSpace);
    for (Index i = 0; i < optimized.getVertexCount(); ++i)
      outputPhi[i] = phiAt(optimized.getVertexCoordinates(i), t);

    auto grid = xdmf.grid("curved-p2-xdmf");
    grid.setMesh(optimized, IO::XDMF::MeshPolicy::Transient);
    grid.clear();
    grid.add("phi", outputPhi, IO::XDMF::Center::Node);

    auto curvedMesh = curvedVisualizationMesh(optimized, 3);
    auto curvedGrid = xdmf.grid("curved-p2-sampled");
    curvedGrid.setMesh(curvedMesh, IO::XDMF::MeshPolicy::Transient);
    curvedGrid.clear();

    auto curvedInterface = curvedInterfacePolyline(optimized, 4);
    auto curveGrid = xdmf.grid("curved-interface");
    curveGrid.setMesh(curvedInterface, IO::XDMF::MeshPolicy::Transient);
    curveGrid.clear();
    xdmf.write(t).flush();

    const auto [qmin, inverted] = meshQuality(optimized);
    const auto [curvedFitRms, curvedFitMax] = curvedInterfaceFit(optimized, t);
    Real fitErr = 0;
    for (Index i = 0; i < optimized.getVertexCount(); ++i)
      if (i < featureMask.size() && featureMask[i])
        fitErr = std::max(fitErr,
            std::abs(phiAt(optimized.getVertexCoordinates(i), t)));
    sumQmin += qmin;
    sumInverted += inverted;
    sumFit += fitErr;
    ++diagSteps;
    std::cout << "step " << step
              << " t=" << t
              << " cells=" << optimized.getCellCount()
              << " qmin=" << qmin
              << " inverted=" << inverted
              << " fit=" << fitErr
              << " curved_fit_rms=" << curvedFitRms
              << " curved_fit_max=" << curvedFitMax
              << " curved_min_jac=" << tmopMinJacobian
              << " tmop_iterations=" << tmopIterations
              << " tmop_last_alpha=" << lastAlpha
              << " tmop_accepted=" << (acceptedTMOP ? 1 : 0)
              << std::endl;

    demoteInterface(optimized);
    optimized.flush();
    background = std::move(optimized);
  }

  xdmf.close();

  if (diagSteps > 0)
    std::cout << "AVERAGE over " << diagSteps << " steps:"
              << " qmin=" << sumQmin / static_cast<Real>(diagSteps)
              << " inverted="
              << static_cast<Real>(sumInverted) / static_cast<Real>(diagSteps)
              << " fit=" << sumFit / static_cast<Real>(diagSteps)
              << " rejected_tmop=" << rejectedTMOP
              << std::endl;
  return 0;
}
