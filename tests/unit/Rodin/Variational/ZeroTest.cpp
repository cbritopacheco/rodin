/**
 * @file ZeroTest.cpp
 * @brief Tests for the Zero function (scalar and vector).
 */
#include <gtest/gtest.h>

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

/// @brief Verifies scalar zero for variational zero by checking tolerance-based numerical results.
TEST(Rodin_Variational_Zero, ScalarZero)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  Zero z;

  EXPECT_NEAR(z.getValue(p), 0.0, 1e-15);
}

/// @brief Verifies scalar copy for variational zero by checking tolerance-based numerical results, copy semantics.
TEST(Rodin_Variational_Zero, ScalarCopy)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  Zero z;
  auto copy = z;

  EXPECT_NEAR(copy.getValue(p), 0.0, 1e-15);
}

/// @brief Verifies vector zero for variational zero by checking tolerance-based numerical results, exact expected values.
TEST(Rodin_Variational_Zero, VectorZero)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  Zero z(3);  // 3D zero vector

  auto val = z.getValue(p);
  EXPECT_EQ(val.size(), 3);
  EXPECT_NEAR(val(0), 0.0, 1e-15);
  EXPECT_NEAR(val(1), 0.0, 1e-15);
  EXPECT_NEAR(val(2), 0.0, 1e-15);
}

/// @brief Verifies vector copy for variational zero by checking tolerance-based numerical results, exact expected values, copy semantics.
TEST(Rodin_Variational_Zero, VectorCopy)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  Zero z(2);
  auto copy = z;

  auto val = copy.getValue(p);
  EXPECT_EQ(val.size(), 2);
  EXPECT_NEAR(val.norm(), 0.0, 1e-15);
}

/// @brief Verifies scalar get order for variational zero by checking exact expected values, true predicates.
TEST(Rodin_Variational_Zero, ScalarGetOrder)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);

  Zero z;
  auto order = z.getOrder(*it);

  ASSERT_TRUE(order.has_value());
  EXPECT_EQ(*order, 0u);
}

/// @brief Verifies vector get order for variational zero by checking exact expected values, true predicates.
TEST(Rodin_Variational_Zero, VectorGetOrder)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);

  Zero z(2);
  auto order = z.getOrder(*it);

  ASSERT_TRUE(order.has_value());
  EXPECT_EQ(*order, 0u);
}
