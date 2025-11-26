/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file P0.cpp
 * @brief Implementation of P0 finite element space static members.
 *
 * This file initializes the static P0Element instances for each supported
 * polytope geometry type (Point, Segment, Triangle, Quadrilateral, 
 * Tetrahedron, Wedge).
 */
#include "P0.h"

namespace Rodin::Variational
{
  const Geometry::GeometryIndexed<P0Element<Real>>
  P0<Real, Geometry::Mesh<Context::Local>>::s_elements =
  {
    { Geometry::Polytope::Type::Point, P0Element<Real>(Geometry::Polytope::Type::Point) },
    { Geometry::Polytope::Type::Segment, P0Element<Real>(Geometry::Polytope::Type::Segment) },
    { Geometry::Polytope::Type::Triangle, P0Element<Real>(Geometry::Polytope::Type::Triangle) },
    { Geometry::Polytope::Type::Quadrilateral, P0Element<Real>(Geometry::Polytope::Type::Quadrilateral) },
    { Geometry::Polytope::Type::Tetrahedron, P0Element<Real>(Geometry::Polytope::Type::Tetrahedron) },
    { Geometry::Polytope::Type::Wedge, P0Element<Real>(Geometry::Polytope::Type::Wedge) }
  };
}
