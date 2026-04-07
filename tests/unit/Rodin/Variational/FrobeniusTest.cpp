/**
 * @file FrobeniusTest.cpp
 * @brief Tests for the Frobenius norm operator.
 */
#include <gtest/gtest.h>
#include <cmath>

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

TEST(Rodin_Variational_Frobenius, Identity2x2)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  IdentityMatrix I(2);
  auto frob = Frobenius(I);

  // ||I_2||_F = sqrt(1² + 1²) = sqrt(2)
  EXPECT_NEAR(frob.getValue(p), std::sqrt(2.0), 1e-10);
}

TEST(Rodin_Variational_Frobenius, Identity3x3)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  IdentityMatrix I(3);
  auto frob = Frobenius(I);

  // ||I_3||_F = sqrt(3)
  EXPECT_NEAR(frob.getValue(p), std::sqrt(3.0), 1e-10);
}

TEST(Rodin_Variational_Frobenius, Copy)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  IdentityMatrix I(2);
  auto frob = Frobenius(I);
  auto copy = frob;

  EXPECT_NEAR(copy.getValue(p), std::sqrt(2.0), 1e-10);
}
