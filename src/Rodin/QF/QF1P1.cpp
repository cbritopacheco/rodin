#include <cmath>

#include "QF1P1.h"

namespace Rodin::QF
{
  const Geometry::GeometryIndexed<Math::SpatialVector<Real>> QF1P1::s_points =
  {
    { Geometry::Polytope::Type::Point,
      Math::SpatialVector<Real>{{0}} },
    { Geometry::Polytope::Type::Segment,
      Math::SpatialVector<Real>{{0.5}} },
    { Geometry::Polytope::Type::Triangle,
      Math::SpatialVector<Real>{{ Real(1) / Real(3), Real(1) / Real(3) }} },
    { Geometry::Polytope::Type::Quadrilateral,
      Math::SpatialVector<Real>{{0.5, 0.5}} },
    { Geometry::Polytope::Type::Tetrahedron,
      Math::SpatialVector<Real>{{0.25, 0.25, 0.25}} }
  };

  const Geometry::GeometryIndexed<Real> QF1P1::s_weights =
  {
    { Geometry::Polytope::Type::Point, 1 },
    { Geometry::Polytope::Type::Segment, 1 },
    { Geometry::Polytope::Type::Triangle, 0.5 },
    { Geometry::Polytope::Type::Quadrilateral, 1 },
    { Geometry::Polytope::Type::Tetrahedron, Real(1.0 / (6 * std::sqrt(2))) }
  };
}
