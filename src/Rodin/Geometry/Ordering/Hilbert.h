/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_ORDERING_HILBERT_H
#define RODIN_GEOMETRY_ORDERING_HILBERT_H

#include "Rodin/Geometry/SpaceFillingCurveOrdering.h"

namespace Rodin::Geometry::Ordering
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
  class Hilbert : public Internal::SpaceFillingCurveOrdering<Hilbert>
  {
    public:
      /**
       * @brief Computes the spatial vertex ordering.
       * @param[in] mesh Local mesh with finite vertex coordinates.
       * @returns Permutation with newToOld[newIndex] equal to the old vertex index.
       */
      IndexVector operator()(const Mesh<Context::Local>& mesh) const;

    private:
      static void transform(std::vector<std::uint64_t>& coordinates, size_t bits);
  };
}

#endif
