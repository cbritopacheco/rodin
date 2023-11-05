#include "P0Element.h"

namespace Rodin::Variational
{
  const Geometry::GeometryIndexed<Math::PointMatrix> ScalarP0Element::s_nodes =
  {
    { Geometry::Polytope::Type::Point,
      Math::PointMatrix{{0}} },
    { Geometry::Polytope::Type::Segment,
      Math::PointMatrix{{0.5}} },
    { Geometry::Polytope::Type::Triangle,
      Math::PointMatrix{{Scalar(1) / Scalar(3), Scalar(1) / Scalar(3)}} },
    { Geometry::Polytope::Type::Quadrilateral,
      Math::PointMatrix{{0.5, 0.5}} },
    { Geometry::Polytope::Type::Tetrahedron,
      Math::PointMatrix{{Scalar(1) / Scalar(3), Scalar(1) / Scalar(3), Scalar(1) / Scalar(3)}} }
  };

  const Geometry::GeometryIndexed<ScalarP0Element::BasisFunction>
  ScalarP0Element::s_basis =
  {
    { Geometry::Polytope::Type::Point, Geometry::Polytope::Type::Point },
    { Geometry::Polytope::Type::Segment, Geometry::Polytope::Type::Segment },
    { Geometry::Polytope::Type::Triangle, Geometry::Polytope::Type::Triangle },
    { Geometry::Polytope::Type::Quadrilateral, Geometry::Polytope::Type::Quadrilateral },
    { Geometry::Polytope::Type::Tetrahedron, Geometry::Polytope::Type::Tetrahedron }
  };

  const Geometry::GeometryIndexed<ScalarP0Element::LinearForm>
  ScalarP0Element::s_ls =
  {
    { Geometry::Polytope::Type::Point, Geometry::Polytope::Type::Point },
    { Geometry::Polytope::Type::Segment, Geometry::Polytope::Type::Segment },
    { Geometry::Polytope::Type::Triangle, Geometry::Polytope::Type::Triangle },
    { Geometry::Polytope::Type::Quadrilateral, Geometry::Polytope::Type::Quadrilateral },
    { Geometry::Polytope::Type::Tetrahedron, Geometry::Polytope::Type::Tetrahedron }
  };

  const Geometry::GeometryIndexed<ScalarP0Element::GradientFunction>
  ScalarP0Element::s_gradient =
  {
    { Geometry::Polytope::Type::Point, Geometry::Polytope::Type::Point },
    { Geometry::Polytope::Type::Segment, Geometry::Polytope::Type::Segment },
    { Geometry::Polytope::Type::Triangle, Geometry::Polytope::Type::Triangle },
    { Geometry::Polytope::Type::Quadrilateral, Geometry::Polytope::Type::Quadrilateral },
    { Geometry::Polytope::Type::Tetrahedron, Geometry::Polytope::Type::Tetrahedron }
  };

  // const std::array<Geometry::GeometryIndexed<Math::PointMatrix>, RODIN_P0_MAX_VECTOR_DIMENSION>
  // VectorP0Element::s_nodes = Internal::initVectorP0Nodes();

  // namespace Internal
  // {
  //   std::array<Geometry::GeometryIndexed<Math::PointMatrix>, RODIN_P0_MAX_VECTOR_DIMENSION>
  //   initVectorP0Nodes()
  //   {
  //     std::array<Geometry::GeometryIndexed<Math::PointMatrix>, RODIN_P0_MAX_VECTOR_DIMENSION> res;
  //     for (size_t vdim = 1; vdim < RODIN_P0_MAX_VECTOR_DIMENSION; vdim++)
  //     {
  //       for (auto g : Geometry::Polytope::Types)
  //       {
  //         switch (g)
  //         {
  //           case Geometry::Polytope::Type::Point:
  //           {
  //             res[vdim][g] = Math::SpatialMatrix::Zero(1, vdim);
  //             break;
  //           }
  //           default:
  //           {
  //             auto sfe = ScalarP0Element(g);
  //             const size_t n = Geometry::Polytope::getVertexCount(g);
  //             const size_t d = Geometry::Polytope::getGeometryDimension(g);
  //             res[vdim][g].resize(d, n * vdim);
  //             for (size_t i = 0; i < n * vdim; i++)
  //               res[vdim][g].col(i) = sfe.getNode(i / vdim);
  //             break;
  //           }
  //         }
  //       }
  //     }
  //     return res;
  //   }
  // }
}
