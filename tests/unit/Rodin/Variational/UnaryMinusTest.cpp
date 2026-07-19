/**
 * @file UnaryMinusTest.cpp
 * @brief Dedicated tests for the UnaryMinus operator.
 */
#include <gtest/gtest.h>
#include <cmath>

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

/// @brief Verifies scalar constant for variational unary minus by checking tolerance-based numerical results.
TEST(Rodin_Variational_UnaryMinus, ScalarConstant)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  RealFunction f(3.0);
  auto neg = -f;

  EXPECT_NEAR(neg.getValue(p), -3.0, 1e-10);
}

/// @brief Verifies scalar zero for variational unary minus by checking tolerance-based numerical results.
TEST(Rodin_Variational_UnaryMinus, ScalarZero)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  RealFunction f(0.0);
  auto neg = -f;

  EXPECT_NEAR(neg.getValue(p), 0.0, 1e-10);
}

/// @brief Verifies scalar negative for variational unary minus by checking tolerance-based numerical results.
TEST(Rodin_Variational_UnaryMinus, ScalarNegative)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  RealFunction f(-7.5);
  auto neg = -f;

  EXPECT_NEAR(neg.getValue(p), 7.5, 1e-10);
}

/// @brief Verifies double negation for variational unary minus by checking tolerance-based numerical results.
TEST(Rodin_Variational_UnaryMinus, DoubleNegation)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  RealFunction f(5.0);
  auto doubleNeg = -(-f);

  EXPECT_NEAR(doubleNeg.getValue(p), 5.0, 1e-10);
}

/// @brief Verifies copy for variational unary minus by checking tolerance-based numerical results, copy semantics.
TEST(Rodin_Variational_UnaryMinus, Copy)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  RealFunction f(4.0);
  auto neg = -f;
  auto copy = neg;

  EXPECT_NEAR(copy.getValue(p), -4.0, 1e-10);
}

/// @brief Verifies get order for variational unary minus by checking exact expected values, true predicates.
TEST(Rodin_Variational_UnaryMinus, GetOrder)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);

  RealFunction f(2.0);
  auto neg = -f;

  auto order = neg.getOrder(*it);
  EXPECT_TRUE(order.has_value());
  EXPECT_EQ(*order, 0u);
}
