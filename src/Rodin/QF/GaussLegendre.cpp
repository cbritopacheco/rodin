#include <cmath>

#include "GaussLegendre.h"

namespace Rodin::QF
{
  const Geometry::GeometryIndexed<Math::SpatialVector<Real>> GaussLegendre::s_points =
  {
    { Geometry::Polytope::Type::Point,
      Math::SpatialVector<Real>{{0}} },
    { Geometry::Polytope::Type::Segment,
      Math::SpatialVector<Real>{{ Real(1 - (1.0 / sqrt(3))) / Real(2), Real(1 + (1.0 / sqrt(3))) / Real(2) }} },
    { Geometry::Polytope::Type::Triangle,
      Math::SpatialVector<Real>{{ Real(1) / Real(3), Real(1) / Real(3) }} },
    { Geometry::Polytope::Type::Quadrilateral,
      Math::SpatialVector<Real>{{0.5, 0.5}} },
    { Geometry::Polytope::Type::Tetrahedron,
      Math::SpatialVector<Real>{{0.25, 0.25, 0.25}} }
  };

  const Geometry::GeometryIndexed<Math::Vector<Real>> GaussLegendre::s_weights =
  {
    { Geometry::Polytope::Type::Point, Math::Vector<Real>{{ 1 }} },
    { Geometry::Polytope::Type::Segment, Math::Vector<Real>{{0.5, 0.5}} },
    { Geometry::Polytope::Type::Triangle, Math::Vector<Real>{{ 0.5 }} },
    { Geometry::Polytope::Type::Quadrilateral, Math::Vector<Real>{{1}} },
    { Geometry::Polytope::Type::Tetrahedron, Math::Vector<Real>{{ Real(1) / Real(6 * std::sqrt(2))}} }
  };
}

