#include "Rodin/Math/Common.h"

#include "QF2P1.h"

namespace Rodin::QF
{
  const Geometry::GeometryIndexed<size_t> QF2P1::s_size =
  {
    { Geometry::Polytope::Type::Point,         1 },
    { Geometry::Polytope::Type::Segment,       1 },
    { Geometry::Polytope::Type::Triangle,      1 },
    { Geometry::Polytope::Type::Quadrilateral, 4 },
    { Geometry::Polytope::Type::Tetrahedron,   1 },
    { Geometry::Polytope::Type::Wedge,         1 }
  };

  const Geometry::GeometryIndexed<std::vector<Math::SpatialVector<Real>>> QF2P1::s_points =
  {
    { Geometry::Polytope::Type::Point,
      { Math::SpatialVector<Real>{{ 0.0 }} } },
    { Geometry::Polytope::Type::Segment,
      { Math::SpatialVector<Real>{{ 0.5 }} } },
    { Geometry::Polytope::Type::Triangle,
      { Math::SpatialVector<Real>{{ Real(1)/Real(3), Real(1)/Real(3) }} } },
    { Geometry::Polytope::Type::Quadrilateral,
      {
        Math::SpatialVector<Real>{{ Real(0.5) * (Real(1) - Real(1) / std::sqrt(3)), Real(0.5) * (Real(1) + Real(1) / std::sqrt(3)) }},
        Math::SpatialVector<Real>{{ Real(0.5) * (Real(1) - Real(1) / std::sqrt(3)), Real(0.5) * (Real(1) - Real(1) / std::sqrt(3)) }},
        Math::SpatialVector<Real>{{ Real(0.5) * (Real(1) + Real(1) / std::sqrt(3)), Real(0.5) * (Real(1) + Real(1) / std::sqrt(3)) }},
        Math::SpatialVector<Real>{{ Real(0.5) * (Real(1) + Real(1) / std::sqrt(3)), Real(0.5) * (Real(1) - Real(1) / std::sqrt(3)) }}
      } },
    { Geometry::Polytope::Type::Tetrahedron,
      { Math::SpatialVector<Real>{{ 0.25, 0.25, 0.25 }} } },
    { Geometry::Polytope::Type::Wedge,
      { Math::SpatialVector<Real>{{ Real(1)/Real(3), Real(1)/Real(3), 0.5 }} } }
  };

  const Geometry::GeometryIndexed<Math::Vector<Real>> QF2P1::s_weights =
  {
    { Geometry::Polytope::Type::Point,         Math::Vector<Real>{{ 1.0 }}         },
    { Geometry::Polytope::Type::Segment,       Math::Vector<Real>{{ 1.0 }}         },
    { Geometry::Polytope::Type::Triangle,      Math::Vector<Real>{{ 0.5 }}         },
    { Geometry::Polytope::Type::Quadrilateral, Math::Vector<Real>{{ 0.25, 0.25, 0.25, 0.25 }} },
    { Geometry::Polytope::Type::Tetrahedron,   Math::Vector<Real>{{ Real(1) / Real(6) }}},
    { Geometry::Polytope::Type::Wedge,         Math::Vector<Real>{{ 0.5 }}         }
  };
}
