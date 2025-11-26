/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Implementation of the Centroid quadrature formula.
 */

#include "Centroid.h"

namespace Rodin::QF
{
  /**
   * @brief Static table of centroid coordinates for each geometry type.
   *
   * The centroid is the barycenter of the reference polytope:
   * - Point: single vertex at origin
   * - Segment: midpoint @f$ (0.5) @f$
   * - Triangle: barycenter @f$ (1/3, 1/3) @f$
   * - Quadrilateral: center @f$ (0.5, 0.5) @f$
   * - Tetrahedron: barycenter @f$ (0.25, 0.25, 0.25) @f$
   * - Wedge (prism): @f$ (1/3, 1/3, 0.5) @f$
   */
  const Geometry::GeometryIndexed<Math::SpatialVector<Real>> Centroid::s_points =
  {
    { Geometry::Polytope::Type::Point,
      Math::SpatialVector<Real>{{ 0 }} },
    { Geometry::Polytope::Type::Segment,
      Math::SpatialVector<Real>{{ 0.5 }} },
    { Geometry::Polytope::Type::Triangle,
      Math::SpatialVector<Real>{{ Real(1) / Real(3), Real(1) / Real(3) }} },
    { Geometry::Polytope::Type::Quadrilateral,
      Math::SpatialVector<Real>{{ 0.5, 0.5 }} },
    { Geometry::Polytope::Type::Tetrahedron,
      Math::SpatialVector<Real>{{ 0.25, 0.25, 0.25 }} },
    { Geometry::Polytope::Type::Wedge,
      Math::SpatialVector<Real>{{ Real(1) / Real(3), Real(1) / Real(3), 0.5 }} }
  };

  /**
   * @brief Static table of quadrature weights for each geometry type.
   *
   * The weight equals the measure (volume) of the reference polytope:
   * - Point: @f$ 1 @f$
   * - Segment: @f$ 1 @f$ (length)
   * - Triangle: @f$ 1/2 @f$ (area)
   * - Quadrilateral: @f$ 1 @f$ (area)
   * - Tetrahedron: @f$ 1/6 @f$ (volume)
   * - Wedge (prism): @f$ 1/2 @f$ (volume)
   */
  const Geometry::GeometryIndexed<Real> Centroid::s_weights =
  {
    { Geometry::Polytope::Type::Point, 1 },
    { Geometry::Polytope::Type::Segment, 1 },
    { Geometry::Polytope::Type::Triangle, 0.5 },
    { Geometry::Polytope::Type::Quadrilateral, 1 },
    { Geometry::Polytope::Type::Tetrahedron, Real(1.0) / Real(6) },
    { Geometry::Polytope::Type::Wedge, 0.5 }
  };
}
