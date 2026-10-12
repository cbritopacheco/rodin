/**
 * @file BoundaryNormalTest.cpp
 * @brief Tests for the BoundaryNormal class.
 */
#include <gtest/gtest.h>
#include <cmath>

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

/// @brief Verifies construction for variational boundary normal by checking exact expected values.
TEST(Rodin_Variational_BoundaryNormal, Construction)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  BoundaryNormal n(mesh);

  EXPECT_EQ(n.getDimension(), mesh.getSpaceDimension());
}

/// @brief Verifies copy for variational boundary normal by checking exact expected values, copy semantics.
TEST(Rodin_Variational_BoundaryNormal, Copy)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  BoundaryNormal n(mesh);
  auto copy = n;

  EXPECT_EQ(copy.getDimension(), n.getDimension());
}

/// @brief Verifies get order for variational boundary normal by checking exact expected values, true predicates.
TEST(Rodin_Variational_BoundaryNormal, GetOrder)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);

  BoundaryNormal n(mesh);
  auto order = n.getOrder(*it);
  EXPECT_TRUE(order.has_value());
  EXPECT_EQ(*order, 0u);
}

/// @brief Verifies unit normal on boundary for variational boundary normal by checking tolerance-based numerical results, true predicates.
TEST(Rodin_Variational_BoundaryNormal, UnitNormalOnBoundary)
{
  // Create a simple mesh and check normal at a boundary face
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  mesh.getConnectivity().compute(1, 2);

  BoundaryNormal n(mesh);

  // Find a boundary edge (codimension-1 entity)
  const size_t dimension = mesh.getDimension();
  bool foundBoundary = false;
  for (auto it = mesh.getPolytope(dimension - 1); !it.end(); ++it)
  {
    const auto& polytope = *it;
    if (mesh.isBoundary(polytope.getIndex()))
    {
      // Create a point on this boundary face
      const Math::Vector<Real> rc{{ 0.5 }};
      Point p(polytope, rc);

      auto val = n.getValue(p);

      // Check that it's a unit vector (norm ≈ 1)
      Real norm = 0.0;
      for (size_t i = 0; i < static_cast<size_t>(val.size()); ++i)
        norm += val(i) * val(i);
      norm = std::sqrt(norm);

      EXPECT_NEAR(norm, 1.0, 1e-10);
      foundBoundary = true;
      break;
    }
  }
  EXPECT_TRUE(foundBoundary);
}

/// @brief On the unit cube, outward normals agree with the coordinate planes.
TEST(Rodin_Variational_BoundaryNormal, AllCubeFacesHaveExactOutwardNormals)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Hexahedron, {2, 2, 2});
  mesh.getConnectivity().compute(2, 3);
  BoundaryNormal normal(mesh);
  size_t count = 0;
  for (auto face = mesh.getBoundary(); face; ++face)
  {
    Point point(*face, Polytope::Traits(face->getGeometry()).getCentroid());
    const auto actual = normal.getValue(point);
    const auto coordinates = point.getCoordinates();
    Math::SpatialVector<Real> expected(3);
    expected.setZero();
    for (size_t axis = 0; axis < 3; ++axis)
    {
      if (std::abs(coordinates[axis]) < 1e-12)
        expected[axis] = -1;
      else if (std::abs(coordinates[axis] - 1) < 1e-12)
        expected[axis] = 1;
    }
    EXPECT_LT((actual - expected).norm(), 1e-12);
    ++count;
  }
  EXPECT_EQ(count, 6);
}

/// @brief Corner normals require a selected incident boundary trace.
TEST(Rodin_Variational_BoundaryNormal, CornerTraceSelectsIncidentEdge)
{
  auto mesh = LocalMesh::Builder()
                .initialize(2)
                .nodes(3)
                .vertex({0, 0})
                .vertex({1, 0})
                .vertex({0, 1})
                .polytope(Polytope::Type::Triangle, {0, 1, 2})
                .finalize();
  mesh.getConnectivity().compute(0, 1);
  mesh.getConnectivity().compute(1, 2);
  for (auto edge = mesh.getPolytope(1); edge; ++edge)
    mesh.setAttribute({1, edge->getIndex()}, edge->getIndex() + 1);
  auto vertex = mesh.getPolytope(0, 0);
  Point corner(*vertex, Math::SpatialPoint(0));
  const auto& incident = mesh.getConnectivity().getIncidence({0, 1}, 0);
  size_t checked = 0;
  for (auto index : incident)
  {
    if (!mesh.isBoundary(index))
      continue;
    auto edge = mesh.getPolytope(1, index);
    Point midpoint(*edge, Polytope::Traits(edge->getGeometry()).getCentroid());
    BoundaryNormal normal(mesh);
    const auto expected = normal.getValue(midpoint);
    normal.traceOf(index + 1);
    EXPECT_LT((normal.getValue(corner) - expected).norm(), 1e-12);
    ++checked;
  }
  EXPECT_EQ(checked, 2);
  BoundaryNormal invalid(mesh);
  invalid.traceOf(999);
  EXPECT_THROW(invalid.getValue(corner), Alert::Exception);
}
