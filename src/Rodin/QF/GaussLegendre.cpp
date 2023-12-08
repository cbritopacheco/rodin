#include <cmath>

#include "GaussLegendre.h"

namespace Rodin::QF
{
  const Geometry::GeometryIndexed<Math::SpatialVector> GaussLegendre::s_points =
  {
    { Geometry::Polytope::Type::Point,
      Math::SpatialVector{{0}} },
    { Geometry::Polytope::Type::Segment,
      Math::SpatialVector{{ (1 - (1.0 / sqrt(3))) / 2, (1 + (1.0 / sqrt(3))) / 2 }} },
    { Geometry::Polytope::Type::Triangle,
      Math::SpatialVector{{ Scalar(1) / Scalar(3), Scalar(1) / Scalar(3) }} },
    { Geometry::Polytope::Type::Quadrilateral,
      Math::SpatialVector{{0.5, 0.5}} },
    { Geometry::Polytope::Type::Tetrahedron,
      Math::SpatialVector{{0.25, 0.25, 0.25}} }
  };

  const Geometry::GeometryIndexed<Math::Vector> GaussLegendre::s_weights =
  {
    { Geometry::Polytope::Type::Point, Math::Vector{{ 1 }} },
    { Geometry::Polytope::Type::Segment, Math::Vector{{0.5, 0.5}} },
    { Geometry::Polytope::Type::Triangle, Math::Vector{{ 0.5 }} },
    { Geometry::Polytope::Type::Quadrilateral, Math::Vector{{1}} },
    { Geometry::Polytope::Type::Tetrahedron, Math::Vector{{1.0 / (6 * std::sqrt(2))}} }
  };
}

