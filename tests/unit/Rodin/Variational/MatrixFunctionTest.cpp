/*
 * @file MatrixFunctionTest.cpp
 * @brief Unit tests for MatrixFunction class.
 */
#include <gtest/gtest.h>
#include <cmath>

#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

static Geometry::Point makePoint(Mesh<Context::Local>& mesh)
{
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  return Geometry::Point(*it, rc);
}

TEST(Rodin_Variational_MatrixFunction, ConstantMatrix)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  Math::Matrix<Real> m(2, 2);
  m << 1.0, 2.0,
       3.0, 4.0;
  MatrixFunction mf(m);

  const auto& val = mf.getValue(p);
  EXPECT_DOUBLE_EQ(val(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(val(0, 1), 2.0);
  EXPECT_DOUBLE_EQ(val(1, 0), 3.0);
  EXPECT_DOUBLE_EQ(val(1, 1), 4.0);
}

TEST(Rodin_Variational_MatrixFunction, RowsAndColumns)
{
  Math::Matrix<Real> m(3, 2);
  m.setZero();
  MatrixFunction mf(m);

  EXPECT_EQ(mf.getRows(), 3u);
  EXPECT_EQ(mf.getColumns(), 2u);
}

TEST(Rodin_Variational_MatrixFunction, Copy)
{
  Math::Matrix<Real> m(2, 2);
  m << 5.0, 6.0,
       7.0, 8.0;
  MatrixFunction mf(m);

  auto* mfCopy = mf.copy();
  ASSERT_NE(mfCopy, nullptr);

  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  const auto& val = mfCopy->getValue(p);
  EXPECT_DOUBLE_EQ(val(0, 0), 5.0);
  EXPECT_DOUBLE_EQ(val(1, 1), 8.0);

  delete mfCopy;
}

TEST(Rodin_Variational_MatrixFunction, GetOrder)
{
  Math::Matrix<Real> m(2, 2);
  m.setIdentity();
  MatrixFunction mf(m);

  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getCell();
  auto polytope = *it;

  auto order = mf.getOrder(polytope);
  ASSERT_TRUE(order.has_value());
  EXPECT_EQ(order.value(), 0u);
}

TEST(Rodin_Variational_MatrixFunction, IdentityMatrix2x2)
{
  Math::Matrix<Real> m(2, 2);
  m.setIdentity();
  MatrixFunction mf(m);

  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto p = makePoint(mesh);

  const auto& val = mf.getValue(p);
  EXPECT_DOUBLE_EQ(val(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(val(0, 1), 0.0);
  EXPECT_DOUBLE_EQ(val(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(val(1, 1), 1.0);
}
