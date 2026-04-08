/**
 * @file RelativeErrorTest.cpp
 * @brief Tests for the RelativeError utility class.
 */
#include <gtest/gtest.h>

#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#include "Rodin/Variational/RelativeError.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

TEST(Rodin_Variational_RelativeError, L2_ExactMatch)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  P1 Vh(mesh);
  GridFunction uh(Vh);
  uh = RealFunction([](const Geometry::Point& p) { return p.x() + p.y(); });

  auto exact = RealFunction([](const Geometry::Point& p) { return p.x() + p.y(); });
  Real err = RelativeError::l2(uh, exact);
  EXPECT_NEAR(err, 0.0, 1e-10);
}

TEST(Rodin_Variational_RelativeError, L1_ExactMatch)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  P1 Vh(mesh);
  GridFunction uh(Vh);
  uh = RealFunction(1.0);

  auto exact = RealFunction(1.0);
  Real err = RelativeError::l1(uh, exact);
  EXPECT_NEAR(err, 0.0, 1e-10);
}

TEST(Rodin_Variational_RelativeError, LInf_ExactMatch)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  P1 Vh(mesh);
  GridFunction uh(Vh);
  uh = RealFunction(2.0);

  auto exact = RealFunction(2.0);
  Real err = RelativeError::lInf(uh, exact);
  EXPECT_NEAR(err, 0.0, 1e-10);
}

TEST(Rodin_Variational_RelativeError, Compute_WithNorm)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  P1 Vh(mesh);
  GridFunction uh(Vh);
  uh = RealFunction(1.0);

  auto exact = RealFunction(1.0);
  Real err = RelativeError::compute(uh, exact, RelativeError::Norm::L2);
  EXPECT_NEAR(err, 0.0, 1e-10);
}

TEST(Rodin_Variational_RelativeError, L2_NonZeroError)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 8, 8 });
  P1 Vh(mesh);
  GridFunction uh(Vh);
  uh = RealFunction(1.0);

  auto exact = RealFunction(2.0);
  Real err = RelativeError::l2(uh, exact);
  EXPECT_NEAR(err, 0.5, 1e-10);
}
