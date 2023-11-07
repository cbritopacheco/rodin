#include <cmath>

#include "QF1P1.h"

namespace Rodin::QF
{
  const Geometry::GeometryIndexed<Math::SpatialVector> QF1P1::s_points =
  {
    { Geometry::Polytope::Type::Point,
      Math::SpatialVector{{0}} },
    { Geometry::Polytope::Type::Segment,
      Math::SpatialVector{{0.5}} },
    { Geometry::Polytope::Type::Triangle,
      Math::SpatialVector{{ Scalar(1) / Scalar(3), Scalar(1) / Scalar(3) }} },
    { Geometry::Polytope::Type::Quadrilateral,
      Math::SpatialVector{{0.5, 0.5}} },
    { Geometry::Polytope::Type::Tetrahedron,
      Math::SpatialVector{{0.25, 0.25, 0.25}} }
  };

  const Geometry::GeometryIndexed<Scalar> QF1P1::s_weights =
  {
    { Geometry::Polytope::Type::Point, 1 },
    { Geometry::Polytope::Type::Segment, 1 },
    { Geometry::Polytope::Type::Triangle, 0.5 },
    { Geometry::Polytope::Type::Quadrilateral, 1 },
    { Geometry::Polytope::Type::Tetrahedron, Scalar(1.0 / (6 * std::sqrt(2))) }
  };
}
