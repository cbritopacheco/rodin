/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_ORDERING_MORTON_H
#define RODIN_GEOMETRY_ORDERING_MORTON_H

#include "Rodin/Geometry/SpaceFillingCurveOrdering.h"

namespace Rodin::Geometry::Ordering
{
  /**
   * @brief Morton (Z-curve) ordering of local mesh vertices.
   *
   * Integer coordinate bits are interleaved, most significant first.
   * Bounding-box normalization, configurable bit depth and stable tie handling
   * are provided by Internal::SpaceFillingCurveOrdering. Construction uses
   * automatic precision. Apply before constructing finite element spaces or
   * transferring a mesh to a device.
   */
  class Morton : public Internal::SpaceFillingCurveOrdering<Morton>
  {
    public:
      /**
       * @brief Computes the spatial vertex ordering.
       * @param[in] mesh Local mesh with finite vertex coordinates.
       * @returns Permutation with newToOld[newIndex] equal to the old vertex index.
       */
      IndexVector operator()(const Mesh<Context::Local>& mesh) const;
  };
}

#endif
