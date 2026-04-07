/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2024.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TESTS_UNIT_VARIATIONAL_INTERFACETESTUTILS_H
#define RODIN_TESTS_UNIT_VARIATIONAL_INTERFACETESTUTILS_H

#include <utility>
#include "Rodin/Variational.h"

namespace Rodin::Tests::Unit
{
  /**
   * @brief Finds an interior face in the mesh and returns a point at its midpoint.
   * @param[in] mesh The mesh (face-to-cell connectivity will be computed)
   * @returns A pair: (true, Point on interior face) or (false, dummy Point)
   */
  inline std::pair<bool, Geometry::Point>
  findInteriorFacePoint(Geometry::Mesh<Context::Local>& mesh)
  {
    const size_t d = mesh.getDimension();
    mesh.getConnectivity().compute(d - 1, d);
    for (auto it = mesh.getFace(); it; ++it)
    {
      if (!it->isBoundary())
      {
        Math::SpatialPoint rc(1);
        rc(0) = 0.5;
        return { true, Geometry::Point(*it, std::move(rc)) };
      }
    }
    return { false, Geometry::Point(*(mesh.getFace()), Math::SpatialPoint(1)) };
  }
}

#endif
