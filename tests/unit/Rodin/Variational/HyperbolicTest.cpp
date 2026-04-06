/**
 * @file HyperbolicTest.cpp
 * @brief Tests for hyperbolic functions: Cosh, Sinh.
 */
#include <gtest/gtest.h>
#include <cmath>

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

// ============================================================
//  Cosh
// ============================================================

TEST(Rodin_Variational_Cosh, Zero)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  RealFunction f(0.0);
  auto result = Cosh(f);

  EXPECT_NEAR(result.getValue(p), 1.0, 1e-10);
}

TEST(Rodin_Variational_Cosh, One)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  RealFunction f(1.0);
  auto result = Cosh(f);

  EXPECT_NEAR(result.getValue(p), std::cosh(1.0), 1e-10);
}

TEST(Rodin_Variational_Cosh, Negative)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  RealFunction f(-2.0);
  auto result = Cosh(f);

  // cosh is even: cosh(-x) = cosh(x)
  EXPECT_NEAR(result.getValue(p), std::cosh(2.0), 1e-10);
}

TEST(Rodin_Variational_Cosh, EvenFunction)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  RealFunction pos(1.5);
  RealFunction neg(-1.5);

  EXPECT_NEAR(Cosh(pos).getValue(p), Cosh(neg).getValue(p), 1e-10);
}

TEST(Rodin_Variational_Cosh, Copy)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  RealFunction f(1.0);
  auto result = Cosh(f);
  auto copy = result;

  EXPECT_NEAR(copy.getValue(p), std::cosh(1.0), 1e-10);
}

TEST(Rodin_Variational_Cosh, HelperFunction)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  RealFunction f(1.0);

  EXPECT_NEAR(cosh(f).getValue(p), std::cosh(1.0), 1e-10);
}

// ============================================================
//  Sinh
// ============================================================

TEST(Rodin_Variational_Sinh, Zero)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  RealFunction f(0.0);
  auto result = Sinh(f);

  EXPECT_NEAR(result.getValue(p), 0.0, 1e-10);
}

TEST(Rodin_Variational_Sinh, One)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  RealFunction f(1.0);
  auto result = Sinh(f);

  EXPECT_NEAR(result.getValue(p), std::sinh(1.0), 1e-10);
}

TEST(Rodin_Variational_Sinh, Negative)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  RealFunction f(-2.0);
  auto result = Sinh(f);

  // sinh is odd: sinh(-x) = -sinh(x)
  EXPECT_NEAR(result.getValue(p), -std::sinh(2.0), 1e-10);
}

TEST(Rodin_Variational_Sinh, OddFunction)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  RealFunction pos(1.5);
  RealFunction neg(-1.5);

  EXPECT_NEAR(Sinh(pos).getValue(p), -Sinh(neg).getValue(p), 1e-10);
}

TEST(Rodin_Variational_Sinh, CoshSinhIdentity)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  RealFunction f(1.0);
  auto c = Cosh(f);
  auto s = Sinh(f);

  // cosh²(x) - sinh²(x) = 1
  auto identity = c * c - s * s;
  EXPECT_NEAR(identity.getValue(p), 1.0, 1e-10);
}

TEST(Rodin_Variational_Sinh, Copy)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  RealFunction f(1.0);
  auto result = Sinh(f);
  auto copy = result;

  EXPECT_NEAR(copy.getValue(p), std::sinh(1.0), 1e-10);
}

TEST(Rodin_Variational_Sinh, HelperFunction)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  RealFunction f(1.0);

  EXPECT_NEAR(sinh(f).getValue(p), std::sinh(1.0), 1e-10);
}
