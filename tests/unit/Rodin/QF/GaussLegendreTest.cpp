/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <cmath>

#include "Rodin/QF/GaussLegendre.h"

using namespace Rodin;
using namespace Rodin::QF;
using namespace Rodin::Geometry;

// Test class for GaussLegendre functionality
class GaussLegendreTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}

    // Helper function to test polynomial integration accuracy
    double integrate_polynomial_1d(const GaussLegendre& qf, int degree)
    {
      double result = 0.0;
      for (size_t i = 0; i < qf.getSize(); ++i)
      {
        const auto& point = qf.getPoint(i);
        double weight = qf.getWeight(i);
        double x = point[0];
        double value = std::pow(x, degree);
        result += weight * value;
      }
      return result;
    }

    // Helper function to compute exact integral of x^n from 0 to 1
    double exact_integral_1d(int n)
    {
      return 1.0 / (n + 1.0);
    }
};

// Test basic construction and properties
TEST_F(GaussLegendreTest, BasicConstruction)
{
  // Test default constructor (order 2)
  GaussLegendre gl_default(Polytope::Type::Segment);
  EXPECT_EQ(gl_default.getSize(), 2);
  EXPECT_EQ(gl_default.getGeometry(), Polytope::Type::Segment);

  // Test explicit order constructor
  GaussLegendre gl_order3(Polytope::Type::Segment, 3);
  EXPECT_EQ(gl_order3.getSize(), 3);
  EXPECT_EQ(gl_order3.getGeometry(), Polytope::Type::Segment);
}

// Test Point geometry
TEST_F(GaussLegendreTest, PointGeometry)
{
  GaussLegendre gl_point(Polytope::Type::Point);

  EXPECT_EQ(gl_point.getSize(), 1);
  EXPECT_EQ(gl_point.getGeometry(), Polytope::Type::Point);

  const auto& point = gl_point.getPoint(0);
  EXPECT_EQ(point.size(), 1);
  EXPECT_DOUBLE_EQ(point[0], 0.0);
  EXPECT_DOUBLE_EQ(gl_point.getWeight(0), 1.0);
}

// Test Segment geometry with different orders
TEST_F(GaussLegendreTest, SegmentGeometry)
{
  // Test order 1 (should have 1 point)
  GaussLegendre gl1(Polytope::Type::Segment, 1);
  EXPECT_EQ(gl1.getSize(), 1);
  EXPECT_DOUBLE_EQ(gl1.getPoint(0)[0], 0.5);  // Midpoint
  EXPECT_DOUBLE_EQ(gl1.getWeight(0), 1.0);

  // Test order 2 (should have 2 points)
  GaussLegendre gl2(Polytope::Type::Segment, 2);
  EXPECT_EQ(gl2.getSize(), 2);

  // Points should be symmetric around 0.5
  double p1 = gl2.getPoint(0)[0];
  double p2 = gl2.getPoint(1)[0];
  EXPECT_NEAR(p1 + p2, 1.0, 1e-14);

  // Weights should be equal
  EXPECT_NEAR(gl2.getWeight(0), gl2.getWeight(1), 1e-14);
  EXPECT_NEAR(gl2.getWeight(0) + gl2.getWeight(1), 1.0, 1e-14);
}

// Test quadrature accuracy for segments
TEST_F(GaussLegendreTest, SegmentQuadratureAccuracy)
{
  // n-point Gauss-Legendre quadrature should integrate polynomials of degree 2n-1 exactly

  // Test 1-point rule (integrates degree 1 exactly)
  GaussLegendre gl1(Polytope::Type::Segment, 1);
  EXPECT_NEAR(integrate_polynomial_1d(gl1, 0), exact_integral_1d(0), 1e-14);
  EXPECT_NEAR(integrate_polynomial_1d(gl1, 1), exact_integral_1d(1), 1e-14);

  // Test 2-point rule (integrates degree 3 exactly)
  GaussLegendre gl2(Polytope::Type::Segment, 2);
  for (int degree = 0; degree <= 3; ++degree)
  {
    EXPECT_NEAR(integrate_polynomial_1d(gl2, degree), exact_integral_1d(degree), 1e-14)
      << "Failed for degree " << degree;
  }

  // Test 3-point rule (integrates degree 5 exactly)
  GaussLegendre gl3(Polytope::Type::Segment, 3);
  for (int degree = 0; degree <= 5; ++degree)
  {
    EXPECT_NEAR(integrate_polynomial_1d(gl3, degree), exact_integral_1d(degree), 1e-13)
      << "Failed for degree " << degree;
  }
}

// Test Quadrilateral geometry
TEST_F(GaussLegendreTest, QuadrilateralGeometry)
{
  // Test uniform order
  GaussLegendre gl_uniform(Polytope::Type::Quadrilateral, 2);
  EXPECT_EQ(gl_uniform.getSize(), 4);  // 2x2 points
  EXPECT_EQ(gl_uniform.getGeometry(), Polytope::Type::Quadrilateral);

  // Test different orders in x and y
  GaussLegendre gl_mixed(Polytope::Type::Quadrilateral, 2, 3);
  EXPECT_EQ(gl_mixed.getSize(), 6);  // 2x3 points

  // Check that all points are 2D
  for (size_t i = 0; i < gl_uniform.getSize(); ++i)
  {
    const auto& point = gl_uniform.getPoint(i);
    EXPECT_EQ(point.size(), 2);
    EXPECT_GE(point[0], 0.0);
    EXPECT_LE(point[0], 1.0);
    EXPECT_GE(point[1], 0.0);
    EXPECT_LE(point[1], 1.0);
  }

  // Check that weights sum to 1
  double total_weight = 0.0;
  for (size_t i = 0; i < gl_uniform.getSize(); ++i)
  {
    total_weight += gl_uniform.getWeight(i);
  }
  EXPECT_NEAR(total_weight, 1.0, 1e-14);
}

// Test Triangle geometry
TEST_F(GaussLegendreTest, TriangleGeometry)
{
  GaussLegendre gl_tri(Polytope::Type::Triangle, 2);
  EXPECT_EQ(gl_tri.getSize(), 4);  // 2x2 points
  EXPECT_EQ(gl_tri.getGeometry(), Polytope::Type::Triangle);

  // Check that all points are 2D and within unit triangle
  for (size_t i = 0; i < gl_tri.getSize(); ++i)
  {
    const auto& point = gl_tri.getPoint(i);
    EXPECT_EQ(point.size(), 2);
    EXPECT_GE(point[0], 0.0);
    EXPECT_GE(point[1], 0.0);
    EXPECT_LE(point[0] + point[1], 1.0 + 1e-14);  // Allow small numerical error
  }

  // Check that weights sum to 0.5 (area of unit triangle)
  double total_weight = 0.0;
  for (size_t i = 0; i < gl_tri.getSize(); ++i)
  {
    total_weight += gl_tri.getWeight(i);
  }
  EXPECT_NEAR(total_weight, 0.5, 1e-13);
}

// Test Tetrahedron geometry
TEST_F(GaussLegendreTest, TetrahedronGeometry)
{
  GaussLegendre gl_tet(Polytope::Type::Tetrahedron, 2, 2, 2);
  EXPECT_EQ(gl_tet.getSize(), 8);  // 2x2x2 points
  EXPECT_EQ(gl_tet.getGeometry(), Polytope::Type::Tetrahedron);

  // Check that all points are 3D and within unit tetrahedron
  for (size_t i = 0; i < gl_tet.getSize(); ++i)
  {
    const auto& point = gl_tet.getPoint(i);
    EXPECT_EQ(point.size(), 3);
    EXPECT_GE(point[0], 0.0);
    EXPECT_GE(point[1], 0.0);
    EXPECT_GE(point[2], 0.0);
    EXPECT_LE(point[0] + point[1] + point[2], 1.0 + 1e-14);  // Allow small numerical error
  }

  // Check that weights sum to 1/6 (volume of unit tetrahedron)
  double total_weight = 0.0;
  for (size_t i = 0; i < gl_tet.getSize(); ++i)
  {
    total_weight += gl_tet.getWeight(i);
  }
  EXPECT_NEAR(total_weight, 1.0/6.0, 1e-12);
}

// Test Wedge geometry
TEST_F(GaussLegendreTest, WedgeGeometry)
{
  GaussLegendre gl_wedge(Polytope::Type::Wedge, 2, 3);
  EXPECT_EQ(gl_wedge.getSize(), 12);  // 2x2x3 points (triangle base * height)
  EXPECT_EQ(gl_wedge.getGeometry(), Polytope::Type::Wedge);

  // Check that all points are 3D
  for (size_t i = 0; i < gl_wedge.getSize(); ++i)
  {
    const auto& point = gl_wedge.getPoint(i);
    EXPECT_EQ(point.size(), 3);
    EXPECT_GE(point[0], 0.0);
    EXPECT_GE(point[1], 0.0);
    EXPECT_GE(point[2], 0.0);
    EXPECT_LE(point[2], 1.0);
    EXPECT_LE(point[0] + point[1], 1.0 + 1e-14);  // Triangle constraint in xy
  }

  // Check that weights sum to 0.5 (volume of unit wedge)
  double total_weight = 0.0;
  for (size_t i = 0; i < gl_wedge.getSize(); ++i)
  {
    total_weight += gl_wedge.getWeight(i);
  }
  EXPECT_NEAR(total_weight, 0.5, 1e-12);
}

// Test copy functionality
TEST_F(GaussLegendreTest, CopyFunctionality)
{
  GaussLegendre original(Polytope::Type::Segment, 3);

  // Test copy constructor
  GaussLegendre copied(original);
  EXPECT_EQ(copied.getSize(), original.getSize());
  EXPECT_EQ(copied.getGeometry(), original.getGeometry());

  for (size_t i = 0; i < original.getSize(); ++i)
  {
    EXPECT_DOUBLE_EQ(copied.getWeight(i), original.getWeight(i));
    const auto& orig_point = original.getPoint(i);
    const auto& copy_point = copied.getPoint(i);
    EXPECT_EQ(copy_point.size(), orig_point.size());
    for (size_t j = 0; j < orig_point.size(); ++j)
    {
      EXPECT_DOUBLE_EQ(copy_point[j], orig_point[j]);
    }
  }

  // Test copy() method
  std::unique_ptr<GaussLegendre> cloned(original.copy());
  EXPECT_EQ(cloned->getSize(), original.getSize());
  EXPECT_EQ(cloned->getGeometry(), original.getGeometry());

  for (size_t i = 0; i < original.getSize(); ++i)
  {
    EXPECT_DOUBLE_EQ(cloned->getWeight(i), original.getWeight(i));
    const auto& orig_point = original.getPoint(i);
    const auto& clone_point = cloned->getPoint(i);
    EXPECT_EQ(clone_point.size(), orig_point.size());
    for (size_t j = 0; j < orig_point.size(); ++j)
    {
      EXPECT_DOUBLE_EQ(clone_point[j], orig_point[j]);
    }
  }
}

// Test symmetry properties
TEST_F(GaussLegendreTest, SymmetryProperties)
{
  // Test that Gauss-Legendre points are symmetric on segments
  GaussLegendre gl(Polytope::Type::Segment, 4);

  for (size_t i = 0; i < gl.getSize() / 2; ++i)
  {
    size_t mirror_i = gl.getSize() - 1 - i;
    double x1 = gl.getPoint(i)[0];
    double x2 = gl.getPoint(mirror_i)[0];
    EXPECT_NEAR(x1 + x2, 1.0, 1e-14);  // Symmetric around 0.5
    EXPECT_NEAR(gl.getWeight(i), gl.getWeight(mirror_i), 1e-14);  // Equal weights
  }
}

// Test edge cases and parameter validation
TEST_F(GaussLegendreTest, EdgeCases)
{
  // Test with minimum order (1)
  GaussLegendre gl_min(Polytope::Type::Segment, 1);
  EXPECT_EQ(gl_min.getSize(), 1);
  EXPECT_DOUBLE_EQ(gl_min.getWeight(0), 1.0);
}
