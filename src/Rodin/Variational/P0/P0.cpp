/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "P0.h"

namespace Rodin::Variational
{
  const Geometry::GeometryIndexed<ScalarP0Element>
  P0<Scalar, Context::Serial, Geometry::Mesh<Context::Serial>>::s_elements =
  {
    { Geometry::Polytope::Type::Point, ScalarP0Element(Geometry::Polytope::Type::Point) },
    { Geometry::Polytope::Type::Segment, ScalarP0Element(Geometry::Polytope::Type::Segment) },
    { Geometry::Polytope::Type::Triangle, ScalarP0Element(Geometry::Polytope::Type::Triangle) },
    { Geometry::Polytope::Type::Quadrilateral, ScalarP0Element(Geometry::Polytope::Type::Quadrilateral) },
    { Geometry::Polytope::Type::Tetrahedron, ScalarP0Element(Geometry::Polytope::Type::Tetrahedron) }
  };
}
