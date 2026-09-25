/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Mesh-hierarchy utilities specific to h-convergence studies.
 */

#ifndef RODIN_TESTS_CONVERGENCE_H_CONVERGENCE_H
#define RODIN_TESTS_CONVERGENCE_H_CONVERGENCE_H

#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <vector>

#include "../Convergence.h"

namespace Rodin::Tests::Convergence::H
{
  class UniformGridHierarchy
  {
    public:
      UniformGridHierarchy(
        Geometry::Polytope::Type geometry, std::initializer_list<size_t> pointsPerAxis)
        : m_grid(geometry),
          m_pointsPerAxis(pointsPerAxis)
      {
        assert(!m_pointsPerAxis.empty());
        assert(std::all_of(m_pointsPerAxis.begin(), m_pointsPerAxis.end(),
          [](size_t points) { return points >= 2; }));
      }

      size_t getDimension() const
      {
        return m_grid.getDimension();
      }

      Geometry::Polytope::Type getGeometry() const
      {
        return m_grid.getGeometry();
      }

      const std::vector<size_t>& getLevels() const
      {
        return m_pointsPerAxis;
      }

      Real getMeshSize(size_t pointsPerAxis) const
      {
        assert(pointsPerAxis >= 2);
        return Real(1) / Real(pointsPerAxis - 1);
      }

      Geometry::LocalMesh makeMesh(size_t pointsPerAxis) const
      {
        return m_grid.makeMesh(pointsPerAxis);
      }

      static std::string_view getGeometryName(Geometry::Polytope::Type geometry)
      {
        return UniformGrid::getGeometryName(geometry);
      }

    private:
      UniformGrid m_grid;
      std::vector<size_t> m_pointsPerAxis;
  };
}

#endif
