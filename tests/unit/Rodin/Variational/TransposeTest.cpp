/**
 * @file TransposeTest.cpp
 * @brief Tests for the Transpose operator on matrix-valued functions.
 */
#include <gtest/gtest.h>

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

TEST(Rodin_Variational_Transpose, IdentityIsSymmetric)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  IdentityMatrix I(2);
  auto It = Transpose(I);

  auto val = It.getValue(p);
  auto orig = I.getValue(p);

  EXPECT_NEAR((val - orig).norm(), 0.0, 1e-10);
}

TEST(Rodin_Variational_Transpose, Copy)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  IdentityMatrix I(3);
  auto It = Transpose(I);
  auto copy = It;

  auto val = copy.getValue(p);
  EXPECT_NEAR(val(0, 0), 1.0, 1e-10);
  EXPECT_NEAR(val(0, 1), 0.0, 1e-10);
}
