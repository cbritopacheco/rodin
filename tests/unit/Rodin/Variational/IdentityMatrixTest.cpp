/**
 * @file IdentityMatrixTest.cpp
 * @brief Tests for the IdentityMatrix function.
 */
#include <gtest/gtest.h>

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

TEST(Rodin_Variational_IdentityMatrix, TwoByTwo)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  IdentityMatrix I(2);
  auto val = I.getValue(p);

  EXPECT_EQ(val.rows(), 2);
  EXPECT_EQ(val.cols(), 2);
  EXPECT_NEAR(val(0, 0), 1.0, 1e-15);
  EXPECT_NEAR(val(0, 1), 0.0, 1e-15);
  EXPECT_NEAR(val(1, 0), 0.0, 1e-15);
  EXPECT_NEAR(val(1, 1), 1.0, 1e-15);
}

TEST(Rodin_Variational_IdentityMatrix, ThreeByThree)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  IdentityMatrix I(3);
  auto val = I.getValue(p);

  EXPECT_EQ(val.rows(), 3);
  EXPECT_EQ(val.cols(), 3);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      EXPECT_NEAR(val(i, j), (i == j) ? 1.0 : 0.0, 1e-15);
}

TEST(Rodin_Variational_IdentityMatrix, Rows)
{
  IdentityMatrix I(4);
  EXPECT_EQ(I.getRows(), 4u);
}

TEST(Rodin_Variational_IdentityMatrix, Columns)
{
  IdentityMatrix I(4);
  EXPECT_EQ(I.getColumns(), 4u);
}

TEST(Rodin_Variational_IdentityMatrix, GetOrder)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);

  IdentityMatrix I(2);
  auto order = I.getOrder(*it);

  ASSERT_TRUE(order.has_value());
  EXPECT_EQ(*order, 0u);
}

TEST(Rodin_Variational_IdentityMatrix, Copy)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  IdentityMatrix I(2);
  auto copy = I;
  auto val = copy.getValue(p);

  EXPECT_NEAR(val(0, 0), 1.0, 1e-15);
  EXPECT_NEAR(val(1, 1), 1.0, 1e-15);
}
