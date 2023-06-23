/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "P1.h"

namespace Rodin::Variational
{
  const Geometry::GeometryIndexed<P1Element<Scalar>>
  P1<Scalar, Context::Serial, Geometry::Mesh<Context::Serial>>::s_elements =
  {
    { Geometry::Polytope::Geometry::Point, ScalarP1Element(Geometry::Polytope::Geometry::Point) },
    { Geometry::Polytope::Geometry::Segment, ScalarP1Element(Geometry::Polytope::Geometry::Segment) },
    { Geometry::Polytope::Geometry::Triangle, ScalarP1Element(Geometry::Polytope::Geometry::Triangle) },
    { Geometry::Polytope::Geometry::Quadrilateral, ScalarP1Element(Geometry::Polytope::Geometry::Quadrilateral) },
    { Geometry::Polytope::Geometry::Tetrahedron, ScalarP1Element(Geometry::Polytope::Geometry::Tetrahedron) }
  };

  P1<Scalar, Context::Serial, Geometry::Mesh<Context::Serial>>
  ::P1(const Geometry::Mesh<Context>& mesh)
    : m_mesh(mesh)
  {}

  P1<Math::Vector, Context::Serial, Geometry::Mesh<Context::Serial>>
  ::P1(const Geometry::Mesh<Context>& mesh, size_t vdim)
    : m_mesh(mesh), m_vdim(vdim)
  {}
}
