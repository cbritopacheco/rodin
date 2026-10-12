/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Orientation and trace contracts for codimension-one normals. */
#include <gtest/gtest.h>
#include "Rodin/Variational/FaceNormal.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  void checkOrientation(
    const FaceNormal& normal, const Point& point, const Polytope& cell)
  {
    const auto value = normal.getValue(point);
    const Point center(cell, Polytope::Traits(cell.getGeometry()).getCentroid());
    EXPECT_NEAR(value.norm(), 1, 1e-12);
    EXPECT_LT(value.dot(center.getCoordinates() - point.getCoordinates()), 0);
    for (Eigen::Index j = 0; j < point.getJacobian().cols(); ++j)
      EXPECT_NEAR(value.dot(point.getJacobian().col(j)), 0, 1e-12);
  }
}

TEST(Rodin_Variational_FaceNormal, BoundaryOrientationAndInteriorTraceAcrossGeometries)
{
  for (const auto geometry :
    {Polytope::Type::Triangle, Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
      Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    auto mesh = LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 3));
    mesh.getConnectivity().compute(dimension - 1, dimension);
    for (auto cell = mesh.getPolytope(dimension); cell; ++cell)
      mesh.setAttribute({dimension, cell->getIndex()}, cell->getIndex() + 1);
    const auto& adjacent = mesh.getConnectivity().getIncidence(dimension - 1, dimension);
    FaceNormal normal(mesh);
    EXPECT_EQ(normal.getDimension(), dimension);
    size_t boundary = 0, interior = 0;
    for (auto face = mesh.getPolytope(dimension - 1); face; ++face)
    {
      const Point point(*face, Polytope::Traits(face->getGeometry()).getCentroid());
      const auto& cells = adjacent.at(face->getIndex());
      EXPECT_EQ(normal.getOrder(*face), 0);
      if (cells.size() == 1)
      {
        checkOrientation(normal, point, *mesh.getPolytope(dimension, cells[0]));
        std::unique_ptr<FaceNormal> clone(normal.copy());
        FaceNormal moved(std::move(*clone));
        EXPECT_LT((normal(point) - moved(point)).norm(), 1e-12);
        ++boundary;
      }
      else
      {
        ASSERT_EQ(cells.size(), 2);
        if (interior == 0)
        {
          EXPECT_THROW(normal.getValue(point), Alert::Exception);
        }
        auto first = normal;
        first.traceOf(cells[0] + 1);
        auto second = normal;
        second.traceOf(cells[1] + 1);
        checkOrientation(first, point, *mesh.getPolytope(dimension, cells[0]));
        checkOrientation(second, point, *mesh.getPolytope(dimension, cells[1]));
        EXPECT_LT((first(point) + second(point)).norm(), 1e-12);
        if (interior == 0)
        {
          auto invalid = normal;
          invalid.traceOf(999999);
          EXPECT_THROW(invalid.getValue(point), Alert::Exception);
        }
        ++interior;
      }
    }
    EXPECT_GT(boundary, 0);
    EXPECT_GT(interior, 0);
  }
}

TEST(Rodin_Variational_FaceNormal, UnsupportedOneDimensionalJacobianRaises)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {2});
  mesh.getConnectivity().compute(0, 1);
  const auto vertex = mesh.getPolytope(0, 0);
  const Point point(*vertex, Math::SpatialPoint(0));
  FaceNormal normal(mesh);
  EXPECT_THROW(normal.getValue(point), Alert::Exception);
}
