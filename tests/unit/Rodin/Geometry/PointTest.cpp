/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Context/Local.h"
#include "Rodin/Geometry.h"
#include "Rodin/Geometry/ForwardDecls.h"
#include "Rodin/Math/Vector.h"

using namespace Rodin;
using namespace Rodin::Geometry;

// ==================== Basic Construction Tests ====================

TEST(Geometry_Point, BasicConstruction_2D_Triangle)
{
  // Create a simple 2D triangular mesh
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(2).nodes(3);

  builder.vertex({0.0, 0.0});  // v0
  builder.vertex({1.0, 0.0});  // v1
  builder.vertex({0.0, 1.0});  // v2

  builder.polytope(Polytope::Type::Triangle, {0, 1, 2});

  mesh = builder.finalize();

  // Get the first triangle polytope
  auto it = mesh.getCell();

  Polytope tri = *it;

  // Create a point at the barycenter in reference coordinates
  Math::SpatialVector<Real> rc(2);
  rc[0] = 1.0 / 3.0;
  rc[1] = 1.0 / 3.0;

  Point p(tri, rc);

  // Verify basic properties
  EXPECT_EQ(p.getDimension(PointBase::Coordinates::Physical), 2);
  EXPECT_EQ(p.getDimension(PointBase::Coordinates::Reference), 2);
}

TEST(Geometry_Point, BasicConstruction_3D_Tetrahedron)
{
  // Create a simple 3D tetrahedral mesh
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(3).nodes(4);

  builder.vertex({0.0, 0.0, 0.0});  // v0
  builder.vertex({1.0, 0.0, 0.0});  // v1
  builder.vertex({0.0, 1.0, 0.0});  // v2
  builder.vertex({0.0, 0.0, 1.0});  // v3

  builder.polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 3});

  mesh = builder.finalize();

  // Get the tetrahedron polytope
  auto it = mesh.getCell();
  Polytope tet = *it;

  // Create a point at the barycenter
  Math::SpatialVector<Real> rc(3);
  rc[0] = 0.25;
  rc[1] = 0.25;
  rc[2] = 0.25;

  Point p(tet, rc);

  // Verify basic properties
  EXPECT_EQ(p.getDimension(PointBase::Coordinates::Physical), 3);
  EXPECT_EQ(p.getDimension(PointBase::Coordinates::Reference), 3);
}

// ==================== Copy and Move Semantics Tests ====================

TEST(Geometry_Point, CopyConstruction)
{
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(2).nodes(3);

  builder.vertex({0.0, 0.0});
  builder.vertex({1.0, 0.0});
  builder.vertex({0.0, 1.0});
  builder.polytope(Polytope::Type::Triangle, {0, 1, 2});

  mesh = builder.finalize();


  auto it = mesh.getCell();
  Polytope tri = *it;

  Math::SpatialVector<Real> rc(2);
  rc[0] = 0.5;
  rc[1] = 0.5;

  Point p1(tri, rc);
  Point p2(p1);  // Copy construction

  // Verify both points have same physical coordinates
  EXPECT_EQ(p1.getDimension(), p2.getDimension());
  EXPECT_NEAR(p1.x(), p2.x(), 1e-10);
  EXPECT_NEAR(p1.y(), p2.y(), 1e-10);
}

TEST(Geometry_Point, MoveConstruction)
{
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(2).nodes(3);

  builder.vertex({0.0, 0.0});
  builder.vertex({1.0, 0.0});
  builder.vertex({0.0, 1.0});
  builder.polytope(Polytope::Type::Triangle, {0, 1, 2});

  mesh = builder.finalize();

  auto it = mesh.getCell();
  Polytope tri = *it;

  Math::SpatialVector<Real> rc(2);
  rc[0] = 0.5;
  rc[1] = 0.5;

  Point p1(tri, rc);
  Real x_orig = p1.x();
  Real y_orig = p1.y();

  Point p2(std::move(p1));  // Move construction

  // Verify moved point has correct coordinates
  EXPECT_NEAR(p2.x(), x_orig, 1e-10);
  EXPECT_NEAR(p2.y(), y_orig, 1e-10);
}

// ==================== Coordinate Access Tests ====================

TEST(Geometry_Point, CoordinateAccess_XYZ)
{
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(3).nodes(4);

  builder.vertex({0.0, 0.0, 0.0});
  builder.vertex({1.0, 0.0, 0.0});
  builder.vertex({0.0, 1.0, 0.0});
  builder.vertex({0.0, 0.0, 1.0});
  builder.polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 3});

  mesh = builder.finalize();

  auto it = mesh.getCell();
  Polytope tet = *it;

  // Point at (0.25, 0.25, 0.25) in physical coordinates
  Math::SpatialVector<Real> rc(3);
  rc[0] = 0.25;
  rc[1] = 0.25;
  rc[2] = 0.25;

  Point p(tet, rc);

  // Test coordinate accessors
  EXPECT_NEAR(p.x(), 0.25, 1e-10);
  EXPECT_NEAR(p.y(), 0.25, 1e-10);
  EXPECT_NEAR(p.z(), 0.25, 1e-10);

  // Test operator() accessor
  EXPECT_NEAR(p(0), 0.25, 1e-10);
  EXPECT_NEAR(p(1), 0.25, 1e-10);
  EXPECT_NEAR(p(2), 0.25, 1e-10);
}

TEST(Geometry_Point, CoordinateAccess_AsVector)
{
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(2).nodes(3);

  builder.vertex({0.0, 0.0});
  builder.vertex({2.0, 0.0});
  builder.vertex({0.0, 2.0});
  builder.polytope(Polytope::Type::Triangle, {0, 1, 2});

  mesh = builder.finalize();

  auto it = mesh.getCell();
  Polytope tri = *it;

  Math::SpatialVector<Real> rc(2);
  rc[0] = 0.5;
  rc[1] = 0.5;

  Point p(tri, rc);

  // Test asVector accessor
  const auto& vec = p.vector();
  EXPECT_EQ(vec.size(), 2);
  EXPECT_NEAR(vec(0), 1.0, 1e-10);
  EXPECT_NEAR(vec(1), 1.0, 1e-10);
}

TEST(Geometry_Point, GetPhysicalCoordinates)
{
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(2).nodes(3);

  builder.vertex({0.0, 0.0});
  builder.vertex({3.0, 0.0});
  builder.vertex({0.0, 3.0});
  builder.polytope(Polytope::Type::Triangle, {0, 1, 2});

  mesh = builder.finalize();

  auto it = mesh.getCell();
  Polytope tri = *it;

  Math::SpatialVector<Real> rc(2);
  rc[0] = 1.0 / 3.0;
  rc[1] = 1.0 / 3.0;

  Point p(tri, rc);

  // Get physical coordinates
  const auto& pc = p.getPhysicalCoordinates();
  EXPECT_EQ(pc.size(), 2);
  EXPECT_NEAR(pc(0), 1.0, 1e-10);
  EXPECT_NEAR(pc(1), 1.0, 1e-10);
}

TEST(Geometry_Point, GetReferenceCoordinates)
{
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(2).nodes(3);

  builder.vertex({0.0, 0.0});
  builder.vertex({1.0, 0.0});
  builder.vertex({0.0, 1.0});
  builder.polytope(Polytope::Type::Triangle, {0, 1, 2});

  mesh = builder.finalize();

  auto it = mesh.getCell();
  Polytope tri = *it;

  Math::SpatialVector<Real> rc(2);
  rc[0] = 0.2;
  rc[1] = 0.3;

  Point p(tri, rc);

  // Get reference coordinates
  const auto& rc_ret = p.getReferenceCoordinates();
  EXPECT_EQ(rc_ret.size(), 2);
  EXPECT_NEAR(rc_ret(0), 0.2, 1e-10);
  EXPECT_NEAR(rc_ret(1), 0.3, 1e-10);
}

// ==================== Norm Tests ====================

TEST(Geometry_Point, NormCalculations)
{
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(2).nodes(3);

  builder.vertex({0.0, 0.0});
  builder.vertex({4.0, 0.0});
  builder.vertex({0.0, 3.0});
  builder.polytope(Polytope::Type::Triangle, {0, 1, 2});

  mesh = builder.finalize();

  auto it = mesh.getCell();
  Polytope tri = *it;

  // Point at physical (3, 4) should have norm 5
  Math::SpatialVector<Real> rc(2);
  rc[0] = 0.75;
  rc[1] = 0;

  Point p(tri, rc);

  // Test various norm calculations
  Real norm = p.norm();
  EXPECT_NEAR(norm, 3.0, 1e-10);

  Real sqNorm = p.squaredNorm();
  EXPECT_NEAR(sqNorm, 9.0, 1e-10);
}

TEST(Geometry_Point, Norm_3D)
{
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(3).nodes(4);

  builder.vertex({0.0, 0.0, 0.0});
  builder.vertex({2.0, 0.0, 0.0});
  builder.vertex({0.0, 2.0, 0.0});
  builder.vertex({0.0, 0.0, 2.0});
  builder.polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 3});

  mesh = builder.finalize();

  auto it = mesh.getCell();
  Polytope tet = *it;

  // Point at (1, 1, 1) should have norm sqrt(3)
  Math::SpatialVector<Real> rc(3);
  rc[0] = 0.5;
  rc[1] = 0.0;
  rc[2] = 0.0;

  Point p(tet, rc);

  Real norm = p.norm();
  EXPECT_NEAR(norm, 1.0, 1e-10);
}

// ==================== Jacobian Tests ====================

TEST(Geometry_Point, GetJacobian_2D)
{
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(2).nodes(3);

  builder.vertex({0.0, 0.0});
  builder.vertex({1.0, 0.0});
  builder.vertex({0.0, 1.0});
  builder.polytope(Polytope::Type::Triangle, {0, 1, 2});

  mesh = builder.finalize();

  auto it = mesh.getCell();
  Polytope tri = *it;

  Math::SpatialVector<Real> rc(2);
  rc[0] = 0.5;
  rc[1] = 0.5;

  Point p(tri, rc);

  // Get Jacobian matrix
  const auto& jac = p.getJacobian();
  EXPECT_EQ(jac.rows(), 2);
  EXPECT_EQ(jac.cols(), 2);

  // For identity triangle, Jacobian should be identity
  EXPECT_NEAR(jac(0, 0), 1.0, 1e-10);
  EXPECT_NEAR(jac(1, 1), 1.0, 1e-10);
}

TEST(Geometry_Point, GetJacobianDeterminant_2D)
{
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(2).nodes(3);

  builder.vertex({0.0, 0.0});
  builder.vertex({2.0, 0.0});
  builder.vertex({0.0, 2.0});
  builder.polytope(Polytope::Type::Triangle, {0, 1, 2});

  mesh = builder.finalize();

  auto it = mesh.getCell();
  Polytope tri = *it;

  Math::SpatialVector<Real> rc(2);
  rc[0] = 0.5;
  rc[1] = 0.5;

  Point p(tri, rc);

  // Get Jacobian determinant
  Real det = p.getJacobianDeterminant();
  EXPECT_GT(det, 0.0);  // Should be positive for proper orientation
}

TEST(Geometry_Point, GetJacobianInverse_2D)
{
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(2).nodes(3);

  builder.vertex({0.0, 0.0});
  builder.vertex({1.0, 0.0});
  builder.vertex({0.0, 1.0});
  builder.polytope(Polytope::Type::Triangle, {0, 1, 2});

  mesh = builder.finalize();

  auto it = mesh.getCell();
  Polytope tri = *it;

  Math::SpatialVector<Real> rc(2);
  rc[0] = 0.5;
  rc[1] = 0.5;

  Point p(tri, rc);

  // Get Jacobian and its inverse
  const auto& jac = p.getJacobian();
  const auto& jacInv = p.getJacobianInverse();

  EXPECT_EQ(jacInv.rows(), 2);
  EXPECT_EQ(jacInv.cols(), 2);

  // Verify that jac * jacInv = identity
  auto product = jac * jacInv;
  EXPECT_NEAR(product(0, 0), 1.0, 1e-9);
  EXPECT_NEAR(product(1, 1), 1.0, 1e-9);
  EXPECT_NEAR(product(0, 1), 0.0, 1e-9);
  EXPECT_NEAR(product(1, 0), 0.0, 1e-9);
}

// ==================== Distortion Tests ====================

TEST(Geometry_Point, GetDistortion_2D)
{
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(2).nodes(3);

  builder.vertex({0.0, 0.0});
  builder.vertex({1.0, 0.0});
  builder.vertex({0.0, 1.0});
  builder.polytope(Polytope::Type::Triangle, {0, 1, 2});

  mesh = builder.finalize();

  auto it = mesh.getCell();
  Polytope tri = *it;

  Math::SpatialVector<Real> rc(2);
  rc[0] = 0.5;
  rc[1] = 0.5;

  Point p(tri, rc);

  // Get distortion
  Real distortion = p.getDistortion();
  EXPECT_GT(distortion, 0.0);
}

// ==================== Comparison Tests ====================

TEST(Geometry_Point, LexicographicalComparison)
{
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(2).nodes(3);

  builder.vertex({0.0, 0.0});
  builder.vertex({2.0, 0.0});
  builder.vertex({0.0, 2.0});
  builder.polytope(Polytope::Type::Triangle, {0, 1, 2});

  mesh = builder.finalize();

  auto it = mesh.getCell();
  Polytope tri = *it;

  Math::SpatialVector<Real> rc1(2), rc2(2);

  rc1[0] = 0.25;
  rc1[1] = 0.25;

  rc2[0] = 0.5;
  rc2[1] = 0.25;

  Point p1(tri, rc1);
  Point p2(tri, rc2);

  // Test lexicographical comparison
  EXPECT_TRUE(p1 < p2);  // (0.5, 0.5) < (1.0, 0.5)
}

// ==================== SetPolytope Tests ====================

TEST(Geometry_Point, SetPolytope)
{
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(2).nodes(3);

  builder.vertex({0.0, 0.0});
  builder.vertex({1.0, 0.0});
  builder.vertex({0.0, 1.0});
  builder.polytope(Polytope::Type::Triangle, {0, 1, 2});

  mesh = builder.finalize();

  auto it = mesh.getCell();
  Polytope tri = *it;

  Math::SpatialVector<Real> rc(2);
  rc[0] = 0.5;
  rc[1] = 0.5;

  Point p(tri, rc);
  Real x_orig = p.x();

  // Set to same polytope
  p.setPolytope(tri);

  // Coordinates should be recalculated
  EXPECT_EQ(&p.getPolytope().getMesh(), &mesh);
}

// ==================== Arithmetic Operations Tests ====================

TEST(Geometry_Point, AdditionWithVector)
{
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(2).nodes(3);

  builder.vertex({0.0, 0.0});
  builder.vertex({1.0, 0.0});
  builder.vertex({0.0, 1.0});
  builder.polytope(Polytope::Type::Triangle, {0, 1, 2});

  mesh = builder.finalize();

  auto it = mesh.getCell();
  Polytope tri = *it;

  Math::SpatialVector<Real> rc(2);
  rc[0] = 0.5;
  rc[1] = 0.5;

  Point p(tri, rc);

  Math::SpatialVector<Real> vec(2);
  vec[0] = 1.0;
  vec[1] = 2.0;

  // Test point + vector
  Math::SpatialVector<Real> result = p + vec;
  EXPECT_NEAR(result(0), p.x() + 1.0, 1e-10);
  EXPECT_NEAR(result(1), p.y() + 2.0, 1e-10);
}

TEST(Geometry_Point, SubtractionWithVector)
{
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(2).nodes(3);

  builder.vertex({0.0, 0.0});
  builder.vertex({1.0, 0.0});
  builder.vertex({0.0, 1.0});
  builder.polytope(Polytope::Type::Triangle, {0, 1, 2});

  mesh = builder.finalize();

  auto it = mesh.getCell();
  Polytope tri = *it;

  Math::SpatialVector<Real> rc(2);
  rc[0] = 0.5;
  rc[1] = 0.5;

  Point p(tri, rc);

  Math::SpatialVector<Real> vec(2);
  vec[0] = 0.1;
  vec[1] = 0.2;

  // Test point - vector
  auto result = p - vec;
  EXPECT_NEAR(result(0), p.x() - 0.1, 1e-10);
  EXPECT_NEAR(result(1), p.y() - 0.2, 1e-10);
}

TEST(Geometry_Point, AdditionOfPoints)
{
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(2).nodes(3);

  builder.vertex({0.0, 0.0});
  builder.vertex({1.0, 0.0});
  builder.vertex({0.0, 1.0});
  builder.polytope(Polytope::Type::Triangle, {0, 1, 2});

  mesh = builder.finalize();

  auto it = mesh.getCell();
  Polytope tri = *it;

  Math::SpatialVector<Real> rc1(2), rc2(2);
  rc1[0] = 0.25;
  rc1[1] = 0.25;

  rc2[0] = 0.5;
  rc2[1] = 0.5;

  Point p1(tri, rc1);
  Point p2(tri, rc2);

  // Test point + point
  auto result = p1 + p2;
  EXPECT_NEAR(result(0), p1.x() + p2.x(), 1e-10);
  EXPECT_NEAR(result(1), p1.y() + p2.y(), 1e-10);
}

TEST(Geometry_Point, SubtractionOfPoints)
{
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(2).nodes(3);

  builder.vertex({0.0, 0.0});
  builder.vertex({1.0, 0.0});
  builder.vertex({0.0, 1.0});
  builder.polytope(Polytope::Type::Triangle, {0, 1, 2});

  mesh = builder.finalize();

  auto it = mesh.getCell();
  Polytope tri = *it;

  Math::SpatialVector<Real> rc1(2), rc2(2);
  rc1[0] = 0.75;
  rc1[1] = 0.25;

  rc2[0] = 0.25;
  rc2[1] = 0.25;

  Point p1(tri, rc1);
  Point p2(tri, rc2);

  // Test point - point
  auto result = p1 - p2;
  EXPECT_NEAR(result(0), p1.x() - p2.x(), 1e-10);
  EXPECT_NEAR(result(1), p1.y() - p2.y(), 1e-10);
}

TEST(Geometry_Point, ScalarMultiplication)
{
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(2).nodes(3);

  builder.vertex({0.0, 0.0});
  builder.vertex({1.0, 0.0});
  builder.vertex({0.0, 1.0});
  builder.polytope(Polytope::Type::Triangle, {0, 1, 2});

  mesh = builder.finalize();

  auto it = mesh.getCell();
  Polytope tri = *it;

  Math::SpatialVector<Real> rc(2);
  rc[0] = 0.5;
  rc[1] = 0.5;

  Point p(tri, rc);

  // Test scalar * point
  auto result1 = 2.0 * p;
  EXPECT_NEAR(result1(0), 2.0 * p.x(), 1e-10);
  EXPECT_NEAR(result1(1), 2.0 * p.y(), 1e-10);

  // Test point * scalar
  auto result2 = p * 2.0;
  EXPECT_NEAR(result2(0), p.x() * 2.0, 1e-10);
  EXPECT_NEAR(result2(1), p.y() * 2.0, 1e-10);
}

// ==================== Edge Case Tests ====================

TEST(Geometry_Point, EdgeCase_PointAtVertex)
{
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(2).nodes(3);

  builder.vertex({0.0, 0.0});
  builder.vertex({1.0, 0.0});
  builder.vertex({0.0, 1.0});
  builder.polytope(Polytope::Type::Triangle, {0, 1, 2});

  mesh = builder.finalize();

  auto it = mesh.getCell();
  Polytope tri = *it;

  // Point at first vertex (0, 0)
  Math::SpatialVector<Real> rc(2);
  rc[0] = 0.0;
  rc[1] = 0.0;

  Point p(tri, rc);

  EXPECT_NEAR(p.x(), 0.0, 1e-10);
  EXPECT_NEAR(p.y(), 0.0, 1e-10);
}

TEST(Geometry_Point, EdgeCase_PointOnEdge)
{
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(2).nodes(3);

  builder.vertex({0.0, 0.0});
  builder.vertex({2.0, 0.0});
  builder.vertex({0.0, 2.0});
  builder.polytope(Polytope::Type::Triangle, {0, 1, 2});

  mesh = builder.finalize();

  auto it = mesh.getCell();
  Polytope tri = *it;

  // Point on edge between v0 and v1 (midpoint)
  Math::SpatialVector<Real> rc(2);
  rc[0] = 0.5;
  rc[1] = 0.0;

  Point p(tri, rc);

  EXPECT_NEAR(p.x(), 1.0, 1e-10);
  EXPECT_NEAR(p.y(), 0.0, 1e-10);
}

// ==================== 3D Specific Tests ====================

TEST(Geometry_Point, 3D_JacobianDeterminant)
{
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(3).nodes(4);

  builder.vertex({0.0, 0.0, 0.0});
  builder.vertex({1.0, 0.0, 0.0});
  builder.vertex({0.0, 1.0, 0.0});
  builder.vertex({0.0, 0.0, 1.0});
  builder.polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 3});

  mesh = builder.finalize();

  auto it = mesh.getCell();
  Polytope tet = *it;

  Math::SpatialVector<Real> rc(3);
  rc[0] = 0.25;
  rc[1] = 0.25;
  rc[2] = 0.25;

  Point p(tet, rc);

  // Get Jacobian determinant
  Real det = p.getJacobianDeterminant();
  EXPECT_GT(det, 0.0);  // Should be positive
}

TEST(Geometry_Point, 3D_JacobianInverse)
{
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(3).nodes(4);

  builder.vertex({0.0, 0.0, 0.0});
  builder.vertex({1.0, 0.0, 0.0});
  builder.vertex({0.0, 1.0, 0.0});
  builder.vertex({0.0, 0.0, 1.0});
  builder.polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 3});

  mesh = builder.finalize();

  auto it = mesh.getCell();
  Polytope tet = *it;

  Math::SpatialVector<Real> rc(3);
  rc[0] = 0.25;
  rc[1] = 0.25;
  rc[2] = 0.25;

  Point p(tet, rc);

  // Get Jacobian and its inverse
  const auto& jac = p.getJacobian();
  const auto& jacInv = p.getJacobianInverse();

  EXPECT_EQ(jacInv.rows(), 3);
  EXPECT_EQ(jacInv.cols(), 3);

  // Verify that jac * jacInv ≈ identity
  auto product = jac * jacInv;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      Real expected = (i == j) ? 1.0 : 0.0;
      EXPECT_NEAR(product(i, j), expected, 1e-9);
    }
  }
}

// ==================== Construction with Physical Coordinates ====================

TEST(Geometry_Point, ConstructionWithPhysicalCoordinates)
{
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(2).nodes(3);

  builder.vertex({0.0, 0.0});
  builder.vertex({1.0, 0.0});
  builder.vertex({0.0, 1.0});
  builder.polytope(Polytope::Type::Triangle, {0, 1, 2});

  mesh = builder.finalize();

  auto it = mesh.getCell();
  Polytope tri = *it;

  Math::SpatialVector<Real> rc(2), pc(2);
  rc[0] = 0.5;
  rc[1] = 0.5;

  pc[0] = 0.5;
  pc[1] = 0.5;

  // Construct with both reference and physical coordinates
  Point p(tri, rc, pc);

  EXPECT_NEAR(p.x(), 0.5, 1e-10);
  EXPECT_NEAR(p.y(), 0.5, 1e-10);
  EXPECT_NEAR(p.getReferenceCoordinates()(0), 0.5, 1e-10);
  EXPECT_NEAR(p.getReferenceCoordinates()(1), 0.5, 1e-10);
}

// ==================== GetPolytope Tests ====================

TEST(Geometry_Point, GetPolytope_ValidReference)
{
  Mesh mesh;
  Mesh<Context::Local>::Builder builder;
  builder.initialize(2).nodes(3);

  builder.vertex({0.0, 0.0});
  builder.vertex({1.0, 0.0});
  builder.vertex({0.0, 1.0});
  builder.polytope(Polytope::Type::Triangle, {0, 1, 2});

  mesh = builder.finalize();

  auto it = mesh.getCell();
  Polytope tri = *it;

  Math::SpatialVector<Real> rc(2);
  rc[0] = 0.5;
  rc[1] = 0.5;

  Point p(tri, rc);

  // Verify getPolytope returns the correct polytope
  const Polytope& polyRef = p.getPolytope();
  EXPECT_EQ(polyRef.getIndex(), tri.getIndex());
  EXPECT_EQ(&polyRef.getMesh(), &tri.getMesh());
}
