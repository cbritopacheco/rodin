/**
 * @file CoordinateFunctionTest.cpp
 * @brief Tests for the F::X, F::Y, F::Z coordinate projection functions.
 */
#include <gtest/gtest.h>
#include <cmath>

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

// --- F::X tests ---

TEST(Rodin_Variational_F_X, ReturnsXCoordinate)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  P1 Vh(mesh);
  GridFunction gf(Vh);
  gf = F::x;

  // Check that F::x projects the x-coordinate at each vertex
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.0, 0.0 }};
  Point p(*it, rc);

  Real val = F::x.getValue(p);
  Real xCoord = p.x();
  EXPECT_NEAR(val, xCoord, 1e-10);
}

TEST(Rodin_Variational_F_X, DifferentPoints)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);

  // Test at multiple reference coordinates within the element
  {
    const Math::Vector<Real> rc{{ 0.1, 0.1 }};
    Point p(*it, rc);
    EXPECT_NEAR(F::x.getValue(p), p.x(), 1e-10);
  }
  {
    const Math::Vector<Real> rc{{ 0.5, 0.2 }};
    Point p(*it, rc);
    EXPECT_NEAR(F::x.getValue(p), p.x(), 1e-10);
  }
}

TEST(Rodin_Variational_F_X, Copy)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  F::X xFunc;
  auto copy = xFunc;
  EXPECT_NEAR(copy.getValue(p), xFunc.getValue(p), 1e-10);
}

TEST(Rodin_Variational_F_X, GetOrder)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);

  auto order = F::x.getOrder(*it);
  EXPECT_TRUE(order.has_value());
  // For affine elements, coordinate projection is order 1
  EXPECT_GE(*order, 1u);
}

// --- F::Y tests ---

TEST(Rodin_Variational_F_Y, ReturnsYCoordinate)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.3, 0.3 }};
  Point p(*it, rc);

  EXPECT_NEAR(F::y.getValue(p), p.y(), 1e-10);
}

TEST(Rodin_Variational_F_Y, Copy)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  F::Y yFunc;
  auto copy = yFunc;
  EXPECT_NEAR(copy.getValue(p), yFunc.getValue(p), 1e-10);
}

// --- Composition test ---

TEST(Rodin_Variational_F, XPlusY)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  auto sum = F::x + F::y;
  Real expected = p.x() + p.y();
  EXPECT_NEAR(sum.getValue(p), expected, 1e-10);
}

TEST(Rodin_Variational_F, XTimesY)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  auto prod = F::x * F::y;
  Real expected = p.x() * p.y();
  EXPECT_NEAR(prod.getValue(p), expected, 1e-10);
}
