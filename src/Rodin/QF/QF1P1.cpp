#include <cmath>

#include "QF1P1.h"

namespace Rodin::QF
{
  const Geometry::GeometryIndexed<std::vector<Math::SpatialVector>> QF1P1::s_points =
  {
    { Geometry::Polytope::Geometry::Point, std::vector{ Math::SpatialVector{{0}} } },
    { Geometry::Polytope::Geometry::Segment, std::vector{ Math::SpatialVector{{0.5}} } },
    { Geometry::Polytope::Geometry::Triangle, std::vector{ Math::SpatialVector{{1.0 / 3.0, 1.0 / 3.0}} } },
    { Geometry::Polytope::Geometry::Quadrilateral, std::vector{ Math::SpatialVector{{0.5, 0.5}} } },
    { Geometry::Polytope::Geometry::Tetrahedron, std::vector{ Math::SpatialVector{{0.25, 0.25, 0.25}} } }
  };

  const Geometry::GeometryIndexed<Math::Vector> QF1P1::s_weights =
  {
    { Geometry::Polytope::Geometry::Point, Math::Vector{{1}} },
    { Geometry::Polytope::Geometry::Segment, Math::Vector{{1}} },
    { Geometry::Polytope::Geometry::Triangle, Math::Vector{{0.5}} },
    { Geometry::Polytope::Geometry::Quadrilateral, Math::Vector{{1}} },
    { Geometry::Polytope::Geometry::Tetrahedron, Math::Vector{{1.0 / (6 * std::sqrt(2))}} }
  };
}
