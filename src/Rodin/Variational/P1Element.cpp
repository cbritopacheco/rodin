#include "P1Element.h"

namespace Rodin::Variational
{
  const Geometry::GeometryIndexed<Math::Matrix> P1Element<Scalar>::s_dofs =
  {
    { Geometry::Polytope::Geometry::Point,
      Math::Matrix{{0}} },
    { Geometry::Polytope::Geometry::Segment,
      Math::Matrix{{0},
                   {1}} },
    { Geometry::Polytope::Geometry::Triangle,
      Math::Matrix{{0, 1, 0},
                   {0, 0, 1}} },
    { Geometry::Polytope::Geometry::Quadrilateral,
      Math::Matrix{{0, 1, 0, 1},
                   {0, 0, 1, 1}} },
    { Geometry::Polytope::Geometry::Tetrahedron,
      Math::Matrix{{0, 1, 0, 0},
                   {0, 0, 1, 0},
                   {0, 0, 0, 1}} },
  };

  const Geometry::GeometryIndexed<Math::Matrix> P1Element<Math::Vector>::s_dofs =
  {
    { Geometry::Polytope::Geometry::Point,
      Math::Matrix{{0}, {0}} },
    { Geometry::Polytope::Geometry::Segment,
      Math::Matrix{{0, 0},
                   {1, 1}} },
    { Geometry::Polytope::Geometry::Triangle,
      Math::Matrix{{0, 0, 1, 1, 0, 0},
                   {0, 0, 0, 0, 1, 1}} },
    { Geometry::Polytope::Geometry::Quadrilateral,
      Math::Matrix{{0, 0, 1, 1, 0, 0, 1, 1},
                   {0, 0, 0, 0, 1, 1, 1, 1}} },
    { Geometry::Polytope::Geometry::Tetrahedron,
      Math::Matrix{{0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                   {0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0},
                   {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0}} },
  };
}

