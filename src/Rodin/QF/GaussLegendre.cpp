#include <cmath>

#include "GaussLegendre.h"

namespace Rodin::QF
{
  const Geometry::GeometryIndexed<Math::SpatialVector<Scalar>> GaussLegendre::s_points =
  {
    { Geometry::Polytope::Type::Point,
      Math::SpatialVector<Scalar>{{0}} },
    { Geometry::Polytope::Type::Segment,
      Math::SpatialVector<Scalar>{{ (1 - (1.0 / sqrt(3))) / 2, (1 + (1.0 / sqrt(3))) / 2 }} },
    { Geometry::Polytope::Type::Triangle,
      Math::SpatialVector<Scalar>{{ Scalar(1) / Scalar(3), Scalar(1) / Scalar(3) }} },
    { Geometry::Polytope::Type::Quadrilateral,
      Math::SpatialVector<Scalar>{{0.5, 0.5}} },
    { Geometry::Polytope::Type::Tetrahedron,
      Math::SpatialVector<Scalar>{{0.25, 0.25, 0.25}} }
  };

  const Geometry::GeometryIndexed<Math::Vector<Scalar>> GaussLegendre::s_weights =
  {
    { Geometry::Polytope::Type::Point, Math::Vector<Scalar>{{ 1 }} },
    { Geometry::Polytope::Type::Segment, Math::Vector<Scalar>{{0.5, 0.5}} },
    { Geometry::Polytope::Type::Triangle, Math::Vector<Scalar>{{ 0.5 }} },
    { Geometry::Polytope::Type::Quadrilateral, Math::Vector<Scalar>{{1}} },
    { Geometry::Polytope::Type::Tetrahedron, Math::Vector<Scalar>{{1.0 / (6 * std::sqrt(2))}} }
  };
}

