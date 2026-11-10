/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_HILBERTORDERING_H
#define RODIN_GEOMETRY_HILBERTORDERING_H

#include "SpaceFillingCurveOrdering.h"

namespace Rodin::Geometry
{
  /**
   * @brief Hilbert ordering of local mesh vertices.
   *
   * Integer coordinates are transformed using Skilling's axes-to-transpose
   * algorithm before their bits are interleaved (doi:10.1063/1.1751381).
   * Bounding-box normalization, configurable bit depth and stable tie handling
   * are provided by Internal::SpaceFillingCurveOrdering. Construction uses
   * automatic precision. Apply before constructing finite element spaces or
   * transferring a mesh to a device.
   */
  class HilbertOrdering : public Internal::SpaceFillingCurveOrdering<HilbertOrdering>
  {
    friend class Internal::SpaceFillingCurveOrdering<HilbertOrdering>;

    private:
      static void transform(std::vector<std::uint64_t>& coordinates, size_t bits);
  };
}

#endif
