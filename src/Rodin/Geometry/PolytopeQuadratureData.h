/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_POLYTOPEQUADRATUREDATA_H
#define RODIN_GEOMETRY_POLYTOPEQUADRATUREDATA_H

/**
 * @file
 * @brief Cached quadrature data for a single polytope and quadrature formula.
 */

#include <vector>

#include "Rodin/QF/ForwardDecls.h"

#include "Point.h"

namespace Rodin::Geometry
{
  /**
   * @brief Cached quadrature data for a single polytope and quadrature formula.
   *
   * Stores the quadrature formula pointer together with the pre-mapped
   * quadrature points.  Each Point carries pre-populated geometric caches
   * (Jacobian, distortion, etc.) so that concurrent const reads are safe
   * without further mutation.
   */
  struct PolytopeQuadratureData
  {
    const QF::QuadratureFormulaBase* qf; ///< Quadrature formula used
    std::vector<Point> points;           ///< Mapped quadrature points
  };
}

#endif // RODIN_GEOMETRY_POLYTOPEQUADRATUREDATA_H
