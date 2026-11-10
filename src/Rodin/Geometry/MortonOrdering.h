/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_MORTONORDERING_H
#define RODIN_GEOMETRY_MORTONORDERING_H

#include "SpaceFillingCurveOrdering.h"

namespace Rodin::Geometry
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
  class MortonOrdering : public Internal::SpaceFillingCurveOrdering<MortonOrdering>
  {
    friend class Internal::SpaceFillingCurveOrdering<MortonOrdering>;

    private:
      static void transform(std::vector<std::uint64_t>&, size_t) {}
  };
}

#endif
