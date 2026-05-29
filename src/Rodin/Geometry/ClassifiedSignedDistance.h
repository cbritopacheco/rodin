/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_CLASSIFIEDSIGNEDDISTANCE_H
#define RODIN_GEOMETRY_CLASSIFIEDSIGNEDDISTANCE_H

#include <vector>

#include "Rodin/Types.h"
#include "Rodin/Math/SpatialVector.h"

#include "ForwardDecls.h"

namespace Rodin::Geometry
{
  /**
   * @brief Brute-force signed distance to the interface of a binary cell cut.
   *
   * This first implementation supports 2D meshes and stores the interface
   * skeleton as line segments.  Distances are negative on cells labelled -1 and
   * positive on cells labelled +1.
   */
  class ClassifiedSignedDistance
  {
    public:
      ClassifiedSignedDistance(
          const MeshBase& mesh,
          std::vector<int> labels,
          std::vector<Index> interfaceFacets);

      Real operator()(const Math::SpatialPoint& p) const;

      Real operator()(const Point& p) const;

      Real distance(const Math::SpatialPoint& p) const;

      int getCellLabel(Index cell) const;

      const std::vector<Index>& getInterfaceFacets() const noexcept
      {
        return m_interfaceFacets;
      }

    private:
      struct Segment
      {
        Index facet;
        Math::SpatialPoint a;
        Math::SpatialPoint b;
      };

      int getNearestCellLabel(const Math::SpatialPoint& p) const;

      static Real pointSegmentDistance(
          const Math::SpatialPoint& p,
          const Math::SpatialPoint& a,
          const Math::SpatialPoint& b);

      size_t m_dimension;
      std::vector<int> m_labels;
      std::vector<Index> m_interfaceFacets;
      std::vector<Segment> m_segments;
      std::vector<Math::SpatialPoint> m_cellCentroids;
  };
}

#endif
