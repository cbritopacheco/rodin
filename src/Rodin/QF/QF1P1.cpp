#include <cmath>

#include "QF1P1.h"

namespace Rodin::QF
{
  const Geometry::GeometryIndexed<std::vector<Math::SpatialVector>> QF1P1::s_points =
  {
    { Geometry::Polytope::Type::Point, std::vector{ Math::SpatialVector{{0}} } },
    { Geometry::Polytope::Type::Segment, std::vector{ Math::SpatialVector{{0.5}} } },
    { Geometry::Polytope::Type::Triangle, std::vector{ Math::SpatialVector{{1.0 / 3.0, 1.0 / 3.0}} } },
    { Geometry::Polytope::Type::Quadrilateral, std::vector{ Math::SpatialVector{{0.5, 0.5}} } },
    { Geometry::Polytope::Type::Tetrahedron, std::vector{ Math::SpatialVector{{0.25, 0.25, 0.25}} } }
  };

  const Geometry::GeometryIndexed<Math::Vector> QF1P1::s_weights =
  {
    { Geometry::Polytope::Type::Point, Math::Vector{{1}} },
    { Geometry::Polytope::Type::Segment, Math::Vector{{1}} },
    { Geometry::Polytope::Type::Triangle, Math::Vector{{0.5}} },
    { Geometry::Polytope::Type::Quadrilateral, Math::Vector{{1}} },
    { Geometry::Polytope::Type::Tetrahedron, Math::Vector{{Scalar(1.0 / (6 * std::sqrt(2)))}} }
  };
}
