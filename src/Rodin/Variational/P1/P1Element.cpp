#include "P1Element.h"

namespace Rodin::Variational
{
  const Geometry::GeometryIndexed<Math::PointMatrix> RealP1Element::s_nodes =
  {
    { Geometry::Polytope::Type::Point,
      Math::PointMatrix{{0}} },
    { Geometry::Polytope::Type::Segment,
      Math::PointMatrix{{0, 1}} },
    { Geometry::Polytope::Type::Triangle,
      Math::PointMatrix{{0, 1, 0},
                        {0, 0, 1}} },
    { Geometry::Polytope::Type::Quadrilateral,
      Math::PointMatrix{{0, 1, 0, 1},
                        {0, 0, 1, 1}} },
    { Geometry::Polytope::Type::Tetrahedron,
      Math::PointMatrix{{0, 1, 0, 0},
                        {0, 0, 1, 0},
                        {0, 0, 0, 1}} },
    { Geometry::Polytope::Type::Wedge,
      Math::PointMatrix{{0, 1, 0, 0, 1, 0},
                        {0, 0, 1, 0, 0, 1},
                        {0, 0, 0, 1, 1, 1}} },
  };

  // ComplexP1Element --------------------------------------------------------

  const Geometry::GeometryIndexed<Math::PointMatrix> ComplexP1Element::s_nodes =
  {
    { Geometry::Polytope::Type::Point,
      Math::PointMatrix{{0}} },
    { Geometry::Polytope::Type::Segment,
      Math::PointMatrix{{0, 1}} },
    { Geometry::Polytope::Type::Triangle,
      Math::PointMatrix{{0, 1, 0},
                        {0, 0, 1}} },
    { Geometry::Polytope::Type::Quadrilateral,
      Math::PointMatrix{{0, 1, 0, 1},
                        {0, 0, 1, 1}} },
    { Geometry::Polytope::Type::Tetrahedron,
      Math::PointMatrix{{0, 1, 0, 0},
                        {0, 0, 1, 0},
                        {0, 0, 0, 1}} },
    { Geometry::Polytope::Type::Wedge,
      Math::PointMatrix{{0, 1, 0, 0, 1, 0},
                        {0, 0, 1, 0, 0, 1},
                        {0, 0, 0, 1, 1, 1}} },
  };

  namespace Internal
  {
    Geometry::GeometryIndexed<Math::PointMatrix> InitVectorP1Nodes()
    {
      Geometry::GeometryIndexed<Math::PointMatrix> res;
      for (auto g : Geometry::Polytope::Types)
      {
        switch (g)
        {
          case Geometry::Polytope::Type::Point:
          {
            res[g] = Math::PointMatrix::Zero(1, 1);
            break;
          }
          default:
          {
            const auto sfe = RealP1Element(g);
            const size_t n = Geometry::Polytope::getVertexCount(g);
            const size_t d = Geometry::Polytope::getGeometryDimension(g);
            res[g].resize(d, n * d);
            for (size_t i = 0; i < n * d; i++)
              res[g].col(i) = sfe.getNode(i / d);
            break;
          }
        }
      }
      return res;
    }
  }
}
