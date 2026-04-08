/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <Eigen/Core>

#include "Rodin/Geometry/Polytope.h"

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  // ---- isTensorProduct ----

  TEST(Geometry_PolytopeTraitsExt, Point_IsTensorProduct)
  {
    Polytope::Traits traits(Polytope::Type::Point);
    // A point is classified as both simplex and tensor product
    EXPECT_TRUE(traits.isSimplex());
    EXPECT_TRUE(traits.isTensorProduct());
  }

  TEST(Geometry_PolytopeTraitsExt, Segment_IsBoth)
  {
    Polytope::Traits traits(Polytope::Type::Segment);
    // A segment is both simplex and tensor product
    EXPECT_TRUE(traits.isSimplex());
    EXPECT_TRUE(traits.isTensorProduct());
  }

  TEST(Geometry_PolytopeTraitsExt, Triangle_IsNotTensorProduct)
  {
    Polytope::Traits traits(Polytope::Type::Triangle);
    EXPECT_TRUE(traits.isSimplex());
    EXPECT_FALSE(traits.isTensorProduct());
  }

  TEST(Geometry_PolytopeTraitsExt, Quadrilateral_IsTensorProduct)
  {
    Polytope::Traits traits(Polytope::Type::Quadrilateral);
    EXPECT_FALSE(traits.isSimplex());
    EXPECT_TRUE(traits.isTensorProduct());
  }

  TEST(Geometry_PolytopeTraitsExt, Tetrahedron_IsNotTensorProduct)
  {
    Polytope::Traits traits(Polytope::Type::Tetrahedron);
    EXPECT_TRUE(traits.isSimplex());
    EXPECT_FALSE(traits.isTensorProduct());
  }

  TEST(Geometry_PolytopeTraitsExt, Hexahedron_IsTensorProduct)
  {
    Polytope::Traits traits(Polytope::Type::Hexahedron);
    EXPECT_FALSE(traits.isSimplex());
    EXPECT_TRUE(traits.isTensorProduct());
  }

  TEST(Geometry_PolytopeTraitsExt, Wedge_IsTensorProduct)
  {
    // Wedge = triangle × segment, so it is a tensor product element
    Polytope::Traits traits(Polytope::Type::Wedge);
    EXPECT_FALSE(traits.isSimplex());
    EXPECT_TRUE(traits.isTensorProduct());
  }

  TEST(Geometry_PolytopeTraitsExt, Hexahedron_Dimension)
  {
    Polytope::Traits traits(Polytope::Type::Hexahedron);
    EXPECT_EQ(traits.getDimension(), 3);
  }

  TEST(Geometry_PolytopeTraitsExt, Hexahedron_VertexCount)
  {
    Polytope::Traits traits(Polytope::Type::Hexahedron);
    EXPECT_EQ(traits.getVertexCount(), 8);
  }

  // ---- getCentroid ----

  TEST(Geometry_PolytopeTraitsExt, PointCentroid)
  {
    Polytope::Traits traits(Polytope::Type::Point);
    const auto& c = traits.getCentroid();
    // Point centroid should be at the origin or the single vertex
    EXPECT_EQ(c.size(), 0);  // 0-dimensional
  }

  TEST(Geometry_PolytopeTraitsExt, SegmentCentroid)
  {
    Polytope::Traits traits(Polytope::Type::Segment);
    const auto& c = traits.getCentroid();
    EXPECT_EQ(c.size(), 1);
    // Reference segment [0,1]: centroid at 0.5
    EXPECT_NEAR(c(0), 0.5, 1e-14);
  }

  TEST(Geometry_PolytopeTraitsExt, TriangleCentroid)
  {
    Polytope::Traits traits(Polytope::Type::Triangle);
    const auto& c = traits.getCentroid();
    EXPECT_EQ(c.size(), 2);
    // Reference triangle with vertices (0,0), (1,0), (0,1): centroid at (1/3, 1/3)
    EXPECT_NEAR(c(0), 1.0 / 3.0, 1e-14);
    EXPECT_NEAR(c(1), 1.0 / 3.0, 1e-14);
  }

  TEST(Geometry_PolytopeTraitsExt, QuadrilateralCentroid)
  {
    Polytope::Traits traits(Polytope::Type::Quadrilateral);
    const auto& c = traits.getCentroid();
    EXPECT_EQ(c.size(), 2);
    // Reference quad [0,1]^2: centroid at (0.5, 0.5)
    EXPECT_NEAR(c(0), 0.5, 1e-14);
    EXPECT_NEAR(c(1), 0.5, 1e-14);
  }

  TEST(Geometry_PolytopeTraitsExt, TetrahedronCentroid)
  {
    Polytope::Traits traits(Polytope::Type::Tetrahedron);
    const auto& c = traits.getCentroid();
    EXPECT_EQ(c.size(), 3);
    // Reference tet: centroid at (1/4, 1/4, 1/4)
    EXPECT_NEAR(c(0), 0.25, 1e-14);
    EXPECT_NEAR(c(1), 0.25, 1e-14);
    EXPECT_NEAR(c(2), 0.25, 1e-14);
  }

  TEST(Geometry_PolytopeTraitsExt, HexahedronCentroid)
  {
    Polytope::Traits traits(Polytope::Type::Hexahedron);
    const auto& c = traits.getCentroid();
    EXPECT_EQ(c.size(), 3);
    // Reference hex [0,1]^3: centroid at (0.5, 0.5, 0.5)
    EXPECT_NEAR(c(0), 0.5, 1e-14);
    EXPECT_NEAR(c(1), 0.5, 1e-14);
    EXPECT_NEAR(c(2), 0.5, 1e-14);
  }

  // ---- getHalfSpace value validation ----

  TEST(Geometry_PolytopeTraitsExt, TriangleHalfSpace_Values)
  {
    Polytope::Traits traits(Polytope::Type::Triangle);
    const auto& hs = traits.getHalfSpace();
    // Reference triangle: 3 constraints (x>=0, y>=0, x+y<=1)
    EXPECT_EQ(hs.matrix.rows(), 3);
    EXPECT_EQ(hs.matrix.cols(), 2);
    EXPECT_EQ(hs.vector.size(), 3);
  }

  TEST(Geometry_PolytopeTraitsExt, QuadrilateralHalfSpace_Values)
  {
    Polytope::Traits traits(Polytope::Type::Quadrilateral);
    const auto& hs = traits.getHalfSpace();
    // Reference quad: 4 constraints (x>=0, x<=1, y>=0, y<=1)
    EXPECT_EQ(hs.matrix.rows(), 4);
    EXPECT_EQ(hs.matrix.cols(), 2);
    EXPECT_EQ(hs.vector.size(), 4);
  }

  TEST(Geometry_PolytopeTraitsExt, TetrahedronHalfSpace_Values)
  {
    Polytope::Traits traits(Polytope::Type::Tetrahedron);
    const auto& hs = traits.getHalfSpace();
    // Reference tet: 4 constraints
    EXPECT_EQ(hs.matrix.rows(), 4);
    EXPECT_EQ(hs.matrix.cols(), 3);
    EXPECT_EQ(hs.vector.size(), 4);
  }

  TEST(Geometry_PolytopeTraitsExt, HexahedronHalfSpace_Values)
  {
    Polytope::Traits traits(Polytope::Type::Hexahedron);
    const auto& hs = traits.getHalfSpace();
    // Reference hex: 6 constraints
    EXPECT_EQ(hs.matrix.rows(), 6);
    EXPECT_EQ(hs.matrix.cols(), 3);
    EXPECT_EQ(hs.vector.size(), 6);
  }

  // ---- Centroid is inside reference element (half-space test) ----

  TEST(Geometry_PolytopeTraitsExt, CentroidInsideReferenceElement_Triangle)
  {
    Polytope::Traits traits(Polytope::Type::Triangle);
    const auto& c = traits.getCentroid();
    const auto& hs = traits.getHalfSpace();
    // Convert centroid to Eigen vector for matrix multiplication
    Eigen::VectorXd cv(c.size());
    for (int j = 0; j < cv.size(); ++j)
      cv(j) = c(j);
    // Ax <= b should hold for the centroid
    Eigen::VectorXd result = hs.matrix * cv;
    for (int i = 0; i < result.size(); ++i)
      EXPECT_LE(result(i), hs.vector(i) + 1e-14);
  }

  TEST(Geometry_PolytopeTraitsExt, CentroidInsideReferenceElement_Quad)
  {
    Polytope::Traits traits(Polytope::Type::Quadrilateral);
    const auto& c = traits.getCentroid();
    const auto& hs = traits.getHalfSpace();
    Eigen::VectorXd cv(c.size());
    for (int j = 0; j < cv.size(); ++j)
      cv(j) = c(j);
    Eigen::VectorXd result = hs.matrix * cv;
    for (int i = 0; i < result.size(); ++i)
      EXPECT_LE(result(i), hs.vector(i) + 1e-14);
  }

  TEST(Geometry_PolytopeTraitsExt, CentroidInsideReferenceElement_Tet)
  {
    Polytope::Traits traits(Polytope::Type::Tetrahedron);
    const auto& c = traits.getCentroid();
    const auto& hs = traits.getHalfSpace();
    Eigen::VectorXd cv(c.size());
    for (int j = 0; j < cv.size(); ++j)
      cv(j) = c(j);
    Eigen::VectorXd result = hs.matrix * cv;
    for (int i = 0; i < result.size(); ++i)
      EXPECT_LE(result(i), hs.vector(i) + 1e-14);
  }

  TEST(Geometry_PolytopeTraitsExt, CentroidInsideReferenceElement_Hex)
  {
    Polytope::Traits traits(Polytope::Type::Hexahedron);
    const auto& c = traits.getCentroid();
    const auto& hs = traits.getHalfSpace();
    Eigen::VectorXd cv(c.size());
    for (int j = 0; j < cv.size(); ++j)
      cv(j) = c(j);
    Eigen::VectorXd result = hs.matrix * cv;
    for (int i = 0; i < result.size(); ++i)
      EXPECT_LE(result(i), hs.vector(i) + 1e-14);
  }
}
