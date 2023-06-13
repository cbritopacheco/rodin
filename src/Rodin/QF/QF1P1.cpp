#include <cmath>

#include "QF1P1.h"

namespace Rodin::QF
{
  const Geometry::GeometryIndexed<Math::Matrix> QF1P1::s_points =
  {
    { Geometry::Polytope::Geometry::Point, Math::Matrix{{0}} },
    { Geometry::Polytope::Geometry::Segment, Math::Matrix{{0.5}} },
    { Geometry::Polytope::Geometry::Triangle, Math::Matrix{{1.0 / 3.0}, {1.0 / 3.0}} },
    { Geometry::Polytope::Geometry::Quadrilateral, Math::Matrix{{0.5}, {0.5}} },
    { Geometry::Polytope::Geometry::Tetrahedron, Math::Matrix{{0.25}, {0.25}, {0.25}} },
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
