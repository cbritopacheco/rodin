/**
 * @file LogicalTest.cpp
 * @brief Tests for logical operators: AND, OR.
 */
#include <gtest/gtest.h>

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  Point makePoint(Mesh<Context::Local>& mesh)
  {
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const Math::Vector<Real> rc{{ 0.25, 0.25 }};
    return Point(*it, rc);
  }
}

// ============================================================
//  AND
// ============================================================

TEST(Rodin_Variational_AND, TrueAndTrue)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(5.0);
  RealFunction b(3.0);
  auto lhs = (a > b);    // true
  auto rhs = (a > 0.0);  // true

  auto result = (lhs && rhs);
  EXPECT_TRUE(result.getValue(p));
}

TEST(Rodin_Variational_AND, TrueAndFalse)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(5.0);
  RealFunction b(3.0);
  auto lhs = (a > b);   // true
  auto rhs = (a < b);   // false

  auto result = (lhs && rhs);
  EXPECT_FALSE(result.getValue(p));
}

TEST(Rodin_Variational_AND, FalseAndFalse)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(1.0);
  RealFunction b(3.0);
  auto lhs = (a > b);  // false
  auto rhs = (a > b);  // false

  auto result = (lhs && rhs);
  EXPECT_FALSE(result.getValue(p));
}

TEST(Rodin_Variational_AND, WithBooleanConstant)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(5.0);
  RealFunction b(3.0);
  auto cond = (a > b);  // true

  EXPECT_TRUE((true && cond).getValue(p));
  EXPECT_FALSE((false && cond).getValue(p));
  EXPECT_TRUE((cond && true).getValue(p));
  EXPECT_FALSE((cond && false).getValue(p));
}

TEST(Rodin_Variational_AND, Copy)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(5.0);
  RealFunction b(3.0);
  auto lhs = (a > b);
  auto rhs = (a > 0.0);
  auto result = (lhs && rhs);
  auto copy = result;

  EXPECT_TRUE(copy.getValue(p));
}

// ============================================================
//  OR
// ============================================================

TEST(Rodin_Variational_OR, TrueOrFalse)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(5.0);
  RealFunction b(3.0);
  auto lhs = (a > b);  // true
  auto rhs = (a < b);  // false

  auto result = (lhs || rhs);
  EXPECT_TRUE(result.getValue(p));
}

TEST(Rodin_Variational_OR, FalseOrFalse)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(1.0);
  RealFunction b(3.0);
  auto lhs = (a > b);  // false
  auto rhs = (a > b);  // false

  auto result = (lhs || rhs);
  EXPECT_FALSE(result.getValue(p));
}

TEST(Rodin_Variational_OR, TrueOrTrue)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(5.0);
  RealFunction b(3.0);
  auto lhs = (a > b);
  auto rhs = (a > 0.0);

  auto result = (lhs || rhs);
  EXPECT_TRUE(result.getValue(p));
}

TEST(Rodin_Variational_OR, WithBooleanConstant)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(1.0);
  RealFunction b(3.0);
  auto cond = (a < b);  // true

  EXPECT_TRUE((true || cond).getValue(p));
  EXPECT_TRUE((cond || false).getValue(p));
  EXPECT_TRUE((false || cond).getValue(p));
}

TEST(Rodin_Variational_OR, Copy)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(5.0);
  RealFunction b(3.0);
  auto lhs = (a > b);
  auto rhs = (a < b);
  auto result = (lhs || rhs);
  auto copy = result;

  EXPECT_TRUE(copy.getValue(p));
}
