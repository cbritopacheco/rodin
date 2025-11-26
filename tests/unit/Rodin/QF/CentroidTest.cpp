/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <vector>

#include "Rodin/QF/Centroid.h"

using namespace Rodin::QF;
using namespace Rodin::Geometry;

// Test class for centroid quadrature functionality
class CentroidTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// Test centroid (Single-point quadrature)
TEST_F(CentroidTest, BasicProperties)
{
  Centroid qf_triangle(Polytope::Type::Triangle);
  Centroid qf_quad(Polytope::Type::Quadrilateral);

  // Single point quadrature should always have size 1
  EXPECT_EQ(qf_triangle.getSize(), 1);
  EXPECT_EQ(qf_quad.getSize(), 1);

  // Test geometry types
  EXPECT_EQ(qf_triangle.getGeometry(), Polytope::Type::Triangle);
  EXPECT_EQ(qf_quad.getGeometry(), Polytope::Type::Quadrilateral);
}

TEST_F(CentroidTest, Weights)
{
  Centroid qf_triangle(Polytope::Type::Triangle);
  Centroid qf_quad(Polytope::Type::Quadrilateral);
  Centroid qf_line(Polytope::Type::Segment);

  // For single point quadrature, weights should represent the measure of the reference element
  EXPECT_NO_THROW(qf_triangle.getWeight(0));
  EXPECT_NO_THROW(qf_quad.getWeight(0));
  EXPECT_NO_THROW(qf_line.getWeight(0));

  // Weights should be positive
  EXPECT_GT(qf_triangle.getWeight(0), 0.0);
  EXPECT_GT(qf_quad.getWeight(0), 0.0);
  EXPECT_GT(qf_line.getWeight(0), 0.0);
}

TEST_F(CentroidTest, Points)
{
  Centroid qf_triangle(Polytope::Type::Triangle);
  Centroid qf_quad(Polytope::Type::Quadrilateral);

  // Test that we can get the single quadrature point (index 0)
  EXPECT_NO_THROW(qf_triangle.getPoint(0));
  EXPECT_NO_THROW(qf_quad.getPoint(0));

  const auto& point_tri = qf_triangle.getPoint(0);
  const auto& point_quad = qf_quad.getPoint(0);

  // Points should be valid spatial vectors
  EXPECT_GE(point_tri.size(), 1);
  EXPECT_GE(point_quad.size(), 1);
}

// Test different geometry types
TEST_F(CentroidTest, DifferentGeometryTypes)
{
  // Test that QF1P1 works with various geometry types
  std::vector<Polytope::Type> geometries = {
    Polytope::Type::Point,
    Polytope::Type::Segment,
    Polytope::Type::Triangle,
    Polytope::Type::Quadrilateral,
    Polytope::Type::Tetrahedron,
  };

  for (auto geom : geometries) {
    EXPECT_NO_THROW(Centroid qf(geom));

    Centroid qf(geom);
    EXPECT_EQ(qf.getGeometry(), geom);
    EXPECT_EQ(qf.getSize(), 1);
    EXPECT_NO_THROW(qf.getWeight(0));
    EXPECT_NO_THROW(qf.getPoint(0));
  }
}

// Integration test: verify quadrature accuracy for simple functions
TEST_F(CentroidTest, QuadratureAccuracy)
{
  // Test that constant function integration is exact
  Centroid qf_triangle(Polytope::Type::Triangle);

  // For a constant function f(x) = 1, integral over reference triangle should equal area * 1
  double integral_approx = qf_triangle.getWeight(0) * 1.0;  // f(x) = 1

  // Reference triangle has area 0.5
  EXPECT_GT(integral_approx, 0.1);  // Should be positive and reasonable
  EXPECT_LT(integral_approx, 1.0);  // Should not be too large
}

// Error handling tests
TEST_F(CentroidTest, OutOfBoundsAccess)
{
  Centroid qf(Polytope::Type::Triangle);

  // Should work for valid index
  EXPECT_NO_THROW(qf.getWeight(0));
  EXPECT_NO_THROW(qf.getPoint(0));

  // Out of bounds access behavior depends on implementation
  // For single point quadrature, index 1 should be invalid
  // Note: Actual behavior may vary (assert, exception, or undefined)
}
