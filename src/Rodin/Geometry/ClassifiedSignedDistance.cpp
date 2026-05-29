/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "ClassifiedSignedDistance.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "Mesh.h"
#include "Point.h"
#include "Polytope.h"

namespace Rodin::Geometry
{
  ClassifiedSignedDistance::ClassifiedSignedDistance(
      const MeshBase& mesh,
      std::vector<int> labels,
      std::vector<Index> interfaceFacets)
    : m_dimension(mesh.getDimension()),
      m_labels(std::move(labels)),
      m_interfaceFacets(std::move(interfaceFacets))
  {
    if (mesh.getDimension() != 2 || mesh.getSpaceDimension() != 2)
      throw std::invalid_argument(
          "ClassifiedSignedDistance currently supports 2D meshes embedded in 2D.");
    if (m_labels.size() != mesh.getCellCount())
      throw std::invalid_argument(
          "ClassifiedSignedDistance labels must match the number of cells.");
    if (m_interfaceFacets.empty())
      throw std::invalid_argument(
          "ClassifiedSignedDistance requires at least one interface facet.");

    m_segments.reserve(m_interfaceFacets.size());
    for (const Index facet : m_interfaceFacets)
    {
      const auto face = mesh.getFace(facet);
      const auto& vertices = face->getVertices();
      if (vertices.size() != 2)
        throw std::invalid_argument(
            "ClassifiedSignedDistance expected segment facets in a 2D mesh.");

      Segment segment;
      segment.facet = facet;
      segment.a = mesh.getVertexCoordinates(vertices[0]);
      segment.b = mesh.getVertexCoordinates(vertices[1]);
      m_segments.push_back(std::move(segment));
    }

    m_cellCentroids.resize(mesh.getCellCount());
    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const auto& vertices = cellIt->getVertices();
      Math::SpatialPoint centroid(2);
      centroid.setZero();
      for (const Index vertex : vertices)
        centroid += mesh.getVertexCoordinates(vertex);
      centroid /= static_cast<Real>(vertices.size());
      m_cellCentroids[cellIt->getIndex()] = std::move(centroid);
    }
  }

  Real ClassifiedSignedDistance::operator()(const Math::SpatialPoint& p) const
  {
    const Real sign = getNearestCellLabel(p) < 0 ? -1 : 1;
    return sign * distance(p);
  }

  Real ClassifiedSignedDistance::operator()(const Point& p) const
  {
    int label = getNearestCellLabel(p.getPhysicalCoordinates());
    const Polytope& polytope = p.getPolytope();
    if (polytope.getDimension() == m_dimension)
      label = getCellLabel(polytope.getIndex());
    const Real sign = label < 0 ? -1 : 1;
    return sign * distance(p.getPhysicalCoordinates());
  }

  Real ClassifiedSignedDistance::distance(const Math::SpatialPoint& p) const
  {
    Real best = std::numeric_limits<Real>::infinity();
    for (const Segment& segment : m_segments)
      best = std::min(best, pointSegmentDistance(p, segment.a, segment.b));
    return best;
  }

  int ClassifiedSignedDistance::getCellLabel(Index cell) const
  {
    if (cell >= m_labels.size())
      throw std::out_of_range("ClassifiedSignedDistance cell index is out of range.");
    return m_labels[cell];
  }

  int ClassifiedSignedDistance::getNearestCellLabel(const Math::SpatialPoint& p) const
  {
    Real best = std::numeric_limits<Real>::infinity();
    Index closest = 0;
    for (Index cell = 0; cell < m_cellCentroids.size(); ++cell)
    {
      const Real dx = p(0) - m_cellCentroids[cell](0);
      const Real dy = p(1) - m_cellCentroids[cell](1);
      const Real squaredDistance = dx * dx + dy * dy;
      if (squaredDistance < best)
      {
        best = squaredDistance;
        closest = cell;
      }
    }
    return getCellLabel(closest);
  }

  Real ClassifiedSignedDistance::pointSegmentDistance(
      const Math::SpatialPoint& p,
      const Math::SpatialPoint& a,
      const Math::SpatialPoint& b)
  {
    const Real vx = b(0) - a(0);
    const Real vy = b(1) - a(1);
    const Real wx = p(0) - a(0);
    const Real wy = p(1) - a(1);
    const Real lengthSquared = vx * vx + vy * vy;
    Real t = 0;
    if (lengthSquared > 0)
      t = std::clamp((wx * vx + wy * vy) / lengthSquared, Real(0), Real(1));
    const Real px = a(0) + t * vx;
    const Real py = a(1) + t * vy;
    const Real dx = p(0) - px;
    const Real dy = p(1) - py;
    return std::sqrt(dx * dx + dy * dy);
  }
}
