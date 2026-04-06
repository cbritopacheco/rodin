/**
 * @file ComparisonTest.cpp
 * @brief Tests for comparison operators: EQ, NEQ, LT, GT, LEQ, GEQ.
 */
#include <gtest/gtest.h>

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  /**
   * @brief Helper to build a mesh and a point for evaluation.
   */
  Point makePoint(Mesh<Context::Local>& mesh)
  {
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const Math::Vector<Real> rc{{ 0.25, 0.25 }};
    return Point(*it, rc);
  }
}

// ============================================================
//  EQ
// ============================================================

TEST(Rodin_Variational_EQ, EqualConstants)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(3.0);
  RealFunction b(3.0);
  auto eq = (a == b);

  EXPECT_TRUE(eq.getValue(p));
}

TEST(Rodin_Variational_EQ, UnequalConstants)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(3.0);
  RealFunction b(4.0);
  auto eq = (a == b);

  EXPECT_FALSE(eq.getValue(p));
}

TEST(Rodin_Variational_EQ, NumberLHS)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(5.0);
  auto eq = (5.0 == a);

  EXPECT_TRUE(eq.getValue(p));
}

TEST(Rodin_Variational_EQ, NumberRHS)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(5.0);
  auto eq = (a == 6.0);

  EXPECT_FALSE(eq.getValue(p));
}

TEST(Rodin_Variational_EQ, Copy)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(7.0);
  RealFunction b(7.0);
  auto eq = (a == b);
  auto eq_copy = eq;

  EXPECT_TRUE(eq_copy.getValue(p));
}

// ============================================================
//  NEQ
// ============================================================

TEST(Rodin_Variational_NEQ, UnequalConstants)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(3.0);
  RealFunction b(4.0);
  auto neq = (a != b);

  EXPECT_TRUE(neq.getValue(p));
}

TEST(Rodin_Variational_NEQ, EqualConstants)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(3.0);
  RealFunction b(3.0);
  auto neq = (a != b);

  EXPECT_FALSE(neq.getValue(p));
}

TEST(Rodin_Variational_NEQ, NumberLHS)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction b(5.0);
  auto neq = (4.0 != b);

  EXPECT_TRUE(neq.getValue(p));
}

TEST(Rodin_Variational_NEQ, NumberRHS)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(5.0);
  auto neq = (a != 5.0);

  EXPECT_FALSE(neq.getValue(p));
}

TEST(Rodin_Variational_NEQ, Copy)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(1.0);
  RealFunction b(2.0);
  auto neq = (a != b);
  auto neq_copy = neq;

  EXPECT_TRUE(neq_copy.getValue(p));
}

// ============================================================
//  LT
// ============================================================

TEST(Rodin_Variational_LT, LessThan)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(2.0);
  RealFunction b(3.0);

  EXPECT_TRUE((a < b).getValue(p));
  EXPECT_FALSE((b < a).getValue(p));
}

TEST(Rodin_Variational_LT, Equal)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(5.0);
  RealFunction b(5.0);

  EXPECT_FALSE((a < b).getValue(p));
}

TEST(Rodin_Variational_LT, NumberOperands)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(2.0);

  EXPECT_TRUE((a < 3.0).getValue(p));
  EXPECT_TRUE((1.0 < a).getValue(p));
  EXPECT_FALSE((a < 1.0).getValue(p));
}

TEST(Rodin_Variational_LT, Copy)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(1.0);
  RealFunction b(2.0);
  auto lt = (a < b);
  auto lt_copy = lt;

  EXPECT_TRUE(lt_copy.getValue(p));
}

// ============================================================
//  GT
// ============================================================

TEST(Rodin_Variational_GT, GreaterThan)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(5.0);
  RealFunction b(3.0);

  EXPECT_TRUE((a > b).getValue(p));
  EXPECT_FALSE((b > a).getValue(p));
}

TEST(Rodin_Variational_GT, Equal)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(5.0);
  RealFunction b(5.0);

  EXPECT_FALSE((a > b).getValue(p));
}

TEST(Rodin_Variational_GT, NumberOperands)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(2.0);

  EXPECT_TRUE((a > 1.0).getValue(p));
  EXPECT_TRUE((3.0 > a).getValue(p));
  EXPECT_FALSE((a > 3.0).getValue(p));
}

TEST(Rodin_Variational_GT, Copy)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(5.0);
  RealFunction b(3.0);
  auto gt = (a > b);
  auto gt_copy = gt;

  EXPECT_TRUE(gt_copy.getValue(p));
}

// ============================================================
//  LEQ
// ============================================================

TEST(Rodin_Variational_LEQ, LessOrEqual)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(2.0);
  RealFunction b(3.0);
  RealFunction c(3.0);

  EXPECT_TRUE((a <= b).getValue(p));
  EXPECT_TRUE((b <= c).getValue(p));   // equal
  EXPECT_FALSE((b <= a).getValue(p));  // greater
}

TEST(Rodin_Variational_LEQ, NumberOperands)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(2.0);

  EXPECT_TRUE((a <= 2.0).getValue(p));
  EXPECT_TRUE((a <= 3.0).getValue(p));
  EXPECT_FALSE((a <= 1.0).getValue(p));
  EXPECT_TRUE((1.0 <= a).getValue(p));
}

TEST(Rodin_Variational_LEQ, Copy)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(2.0);
  RealFunction b(3.0);
  auto leq = (a <= b);
  auto leq_copy = leq;

  EXPECT_TRUE(leq_copy.getValue(p));
}

// ============================================================
//  GEQ
// ============================================================

TEST(Rodin_Variational_GEQ, GreaterOrEqual)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(5.0);
  RealFunction b(3.0);
  RealFunction c(5.0);

  EXPECT_TRUE((a >= b).getValue(p));
  EXPECT_TRUE((a >= c).getValue(p));   // equal
  EXPECT_FALSE((b >= a).getValue(p));  // less
}

TEST(Rodin_Variational_GEQ, NumberOperands)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(3.0);

  EXPECT_TRUE((a >= 3.0).getValue(p));
  EXPECT_TRUE((a >= 2.0).getValue(p));
  EXPECT_FALSE((a >= 4.0).getValue(p));
  EXPECT_TRUE((4.0 >= a).getValue(p));
}

TEST(Rodin_Variational_GEQ, Copy)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  RealFunction a(5.0);
  RealFunction b(3.0);
  auto geq = (a >= b);
  auto geq_copy = geq;

  EXPECT_TRUE(geq_copy.getValue(p));
}
