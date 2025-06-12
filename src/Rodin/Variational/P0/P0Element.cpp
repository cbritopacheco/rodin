#include "P0Element.h"

namespace Rodin::Variational
{
  const Geometry::GeometryIndexed<Math::PointMatrix> RealP0Element::s_nodes =
  {
    { Geometry::Polytope::Type::Point,
      Math::PointMatrix{{ 0 }} },
    { Geometry::Polytope::Type::Segment,
      Math::PointMatrix{{ 0.5 }} },
    { Geometry::Polytope::Type::Triangle,
      Math::PointMatrix{{ Real(1) / Real(3) },
                        { Real(1) / Real(3) }}},
    { Geometry::Polytope::Type::Quadrilateral,
      Math::PointMatrix{{ 0.5 },
                        { 0.5 }} },
    { Geometry::Polytope::Type::Tetrahedron,
      Math::PointMatrix{{ 0.25 },
                        { 0.25 },
                        { 0.25 }} },
    { Geometry::Polytope::Type::Wedge,
      Math::PointMatrix{{ Real(1) / Real(3) },
                        { Real(1) / Real(3) },
                        { 0.5 }} }
  };

  const Geometry::GeometryIndexed<RealP0Element::BasisFunction>
  RealP0Element::s_basis =
  {
    { Geometry::Polytope::Type::Point, Geometry::Polytope::Type::Point },
    { Geometry::Polytope::Type::Segment, Geometry::Polytope::Type::Segment },
    { Geometry::Polytope::Type::Triangle, Geometry::Polytope::Type::Triangle },
    { Geometry::Polytope::Type::Quadrilateral, Geometry::Polytope::Type::Quadrilateral },
    { Geometry::Polytope::Type::Tetrahedron, Geometry::Polytope::Type::Tetrahedron },
    { Geometry::Polytope::Type::Wedge, Geometry::Polytope::Type::Wedge }
  };

  const Geometry::GeometryIndexed<RealP0Element::LinearForm>
  RealP0Element::s_ls =
  {
    { Geometry::Polytope::Type::Point, Geometry::Polytope::Type::Point },
    { Geometry::Polytope::Type::Segment, Geometry::Polytope::Type::Segment },
    { Geometry::Polytope::Type::Triangle, Geometry::Polytope::Type::Triangle },
    { Geometry::Polytope::Type::Quadrilateral, Geometry::Polytope::Type::Quadrilateral },
    { Geometry::Polytope::Type::Tetrahedron, Geometry::Polytope::Type::Tetrahedron },
    { Geometry::Polytope::Type::Wedge, Geometry::Polytope::Type::Wedge }
  };

  const Geometry::GeometryIndexed<RealP0Element::GradientFunction>
  RealP0Element::s_gradient =
  {
    { Geometry::Polytope::Type::Point, Geometry::Polytope::Type::Point },
    { Geometry::Polytope::Type::Segment, Geometry::Polytope::Type::Segment },
    { Geometry::Polytope::Type::Triangle, Geometry::Polytope::Type::Triangle },
    { Geometry::Polytope::Type::Quadrilateral, Geometry::Polytope::Type::Quadrilateral },
    { Geometry::Polytope::Type::Tetrahedron, Geometry::Polytope::Type::Tetrahedron },
    { Geometry::Polytope::Type::Wedge, Geometry::Polytope::Type::Wedge }
  };
}
