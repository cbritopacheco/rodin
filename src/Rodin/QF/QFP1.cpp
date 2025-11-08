#include "Rodin/Math/Vector.h"
#include "QFP1.h"

namespace Rodin::QF
{
  const Geometry::GeometryIndexed<size_t> QFVertexP1::s_size =
  {
    { Geometry::Polytope::Type::Point,         1 },
    { Geometry::Polytope::Type::Segment,       2 },
    { Geometry::Polytope::Type::Triangle,      3 },
    { Geometry::Polytope::Type::Quadrilateral, 4 },
    { Geometry::Polytope::Type::Tetrahedron,   4 },
    { Geometry::Polytope::Type::Wedge,         6 }
  };

  const Geometry::GeometryIndexed<std::vector<Math::SpatialVector<Real>>> QFVertexP1::s_points =
  {
    { Geometry::Polytope::Type::Point,
      { Math::SpatialVector<Real>{{ 0.0 }} } },

    { Geometry::Polytope::Type::Segment,
      {
        Math::SpatialVector<Real>{{ 0.0 }},
        Math::SpatialVector<Real>{{ 1.0 }}
      } },

    { Geometry::Polytope::Type::Triangle,
      {
        Math::SpatialVector<Real>{{ 0.0, 0.0 }},
        Math::SpatialVector<Real>{{ 1.0, 0.0 }},
        Math::SpatialVector<Real>{{ 0.0, 1.0 }}
      } },

    { Geometry::Polytope::Type::Quadrilateral,
      {
        Math::SpatialVector<Real>{{ 0.0, 0.0 }},
        Math::SpatialVector<Real>{{ 1.0, 0.0 }},
        Math::SpatialVector<Real>{{ 1.0, 1.0 }},
        Math::SpatialVector<Real>{{ 0.0, 1.0 }}
      } },

    { Geometry::Polytope::Type::Tetrahedron,
      {
        Math::SpatialVector<Real>{{ 0.0, 0.0, 0.0 }},
        Math::SpatialVector<Real>{{ 1.0, 0.0, 0.0 }},
        Math::SpatialVector<Real>{{ 0.0, 1.0, 0.0 }},
        Math::SpatialVector<Real>{{ 0.0, 0.0, 1.0 }}
      } },

    { Geometry::Polytope::Type::Wedge, // triangle x segment
      {
        Math::SpatialVector<Real>{{ 0.0, 0.0, 0.0 }},
        Math::SpatialVector<Real>{{ 1.0, 0.0, 0.0 }},
        Math::SpatialVector<Real>{{ 0.0, 1.0, 0.0 }},
        Math::SpatialVector<Real>{{ 0.0, 0.0, 1.0 }},
        Math::SpatialVector<Real>{{ 1.0, 0.0, 1.0 }},
        Math::SpatialVector<Real>{{ 0.0, 1.0, 1.0 }}
      } }
  };

  const Geometry::GeometryIndexed<Math::Vector<Real>> QFVertexP1::s_weights =
  {
    { Geometry::Polytope::Type::Point,
      Math::Vector<Real>{{ 1.0 }} },

    { Geometry::Polytope::Type::Segment,
      Math::Vector<Real>{{ 0.5, 0.5 }} },

    { Geometry::Polytope::Type::Triangle,
      Math::Vector<Real>{{ Real(1)/Real(6),
                           Real(1)/Real(6),
                           Real(1)/Real(6) }} },

    { Geometry::Polytope::Type::Quadrilateral,
      Math::Vector<Real>{{ 0.25, 0.25, 0.25, 0.25 }} },

    { Geometry::Polytope::Type::Tetrahedron,
      Math::Vector<Real>{{ Real(1)/Real(24),
                           Real(1)/Real(24),
                           Real(1)/Real(24),
                           Real(1)/Real(24) }} },

    { Geometry::Polytope::Type::Wedge,
      Math::Vector<Real>{{ Real(1)/Real(12),
                           Real(1)/Real(12),
                           Real(1)/Real(12),
                           Real(1)/Real(12),
                           Real(1)/Real(12),
                           Real(1)/Real(12) }} }
  };
}
