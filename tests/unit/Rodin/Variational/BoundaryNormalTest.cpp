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

TEST(Rodin_Variational_BoundaryNormal, Construction)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  BoundaryNormal n(mesh);

  EXPECT_EQ(n.getDimension(), mesh.getSpaceDimension());
}

TEST(Rodin_Variational_BoundaryNormal, Copy)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  BoundaryNormal n(mesh);
  auto copy = n;

  EXPECT_EQ(copy.getDimension(), n.getDimension());
}

TEST(Rodin_Variational_BoundaryNormal, GetOrder)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);

  BoundaryNormal n(mesh);
  auto order = n.getOrder(*it);
  EXPECT_TRUE(order.has_value());
  EXPECT_EQ(*order, 0u);
}

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
