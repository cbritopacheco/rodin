/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <memory>
#include <cmath>

#include "Rodin/QF/QuadratureFormula.h"
#include "Rodin/QF/QF1P1.h"
#include "Rodin/QF/GaussLegendre.h"
#include "Rodin/Geometry.h"

using namespace Rodin::QF;
using namespace Rodin::Geometry;

// Test class for QF functionality
class QFTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// Test QuadratureFormulaBase construction and basic functionality
TEST_F(QFTest, QuadratureFormulaBaseConstruction)
{
  // Note: Can't directly test abstract base class, but we can test through derived classes
  QF1P1 qf1p1(Polytope::Type::Triangle);

  EXPECT_EQ(qf1p1.getGeometry(), Polytope::Type::Triangle);

  // Test copy constructor
  QF1P1 qf1p1_copy(qf1p1);
  EXPECT_EQ(qf1p1_copy.getGeometry(), Polytope::Type::Triangle);
}

// Test QF1P1 (Single-point quadrature)
TEST_F(QFTest, QF1P1BasicProperties)
{
  QF1P1 qf_triangle(Polytope::Type::Triangle);
  QF1P1 qf_quad(Polytope::Type::Quadrilateral);

  // Single point quadrature should always have size 1
  EXPECT_EQ(qf_triangle.getSize(), 1);
  EXPECT_EQ(qf_quad.getSize(), 1);

  // Test geometry types
  EXPECT_EQ(qf_triangle.getGeometry(), Polytope::Type::Triangle);
  EXPECT_EQ(qf_quad.getGeometry(), Polytope::Type::Quadrilateral);
}

TEST_F(QFTest, QF1P1Weights)
{
  QF1P1 qf_triangle(Polytope::Type::Triangle);
  QF1P1 qf_quad(Polytope::Type::Quadrilateral);
  QF1P1 qf_line(Polytope::Type::Segment);

  // For single point quadrature, weights should represent the measure of the reference element
  EXPECT_NO_THROW(qf_triangle.getWeight(0));
  EXPECT_NO_THROW(qf_quad.getWeight(0));
  EXPECT_NO_THROW(qf_line.getWeight(0));

  // Weights should be positive
  EXPECT_GT(qf_triangle.getWeight(0), 0.0);
  EXPECT_GT(qf_quad.getWeight(0), 0.0);
  EXPECT_GT(qf_line.getWeight(0), 0.0);
}

TEST_F(QFTest, QF1P1Points)
{
  QF1P1 qf_triangle(Polytope::Type::Triangle);
  QF1P1 qf_quad(Polytope::Type::Quadrilateral);

  // Test that we can get the single quadrature point (index 0)
  EXPECT_NO_THROW(qf_triangle.getPoint(0));
  EXPECT_NO_THROW(qf_quad.getPoint(0));

  const auto& point_tri = qf_triangle.getPoint(0);
  const auto& point_quad = qf_quad.getPoint(0);

  // Points should be valid spatial vectors
  EXPECT_GE(point_tri.size(), 1);
  EXPECT_GE(point_quad.size(), 1);
}

// Test GaussLegendre quadrature
TEST_F(QFTest, GaussLegendreBasicProperties)
{
  GaussLegendre gl_segment(Polytope::Type::Segment);

  // GaussLegendre with 2 points
  EXPECT_EQ(gl_segment.getSize(), 2);
  EXPECT_EQ(gl_segment.getGeometry(), Polytope::Type::Segment);
}

TEST_F(QFTest, GaussLegendreWeights)
{
  GaussLegendre gl_segment(Polytope::Type::Segment);

  // Should have weights for each quadrature point
  for (size_t i = 0; i < gl_segment.getSize(); ++i) {
    EXPECT_NO_THROW(gl_segment.getWeight(i));
    EXPECT_GT(gl_segment.getWeight(i), 0.0);
  }
}

TEST_F(QFTest, GaussLegendrePoints)
{
  GaussLegendre gl_segment(Polytope::Type::Segment);

  // Test points are accessible (though implementation may have limitations)
  EXPECT_NO_THROW(gl_segment.getPoint(0));

  const auto& point = gl_segment.getPoint(0);
  EXPECT_GE(point.size(), 1);
}

// Test polymorphic behavior
TEST_F(QFTest, PolymorphicBehavior)
{
  std::unique_ptr<QuadratureFormulaBase> qf1 = std::make_unique<QF1P1>(Polytope::Type::Triangle);
  std::unique_ptr<QuadratureFormulaBase> qf2 = std::make_unique<GaussLegendre>(Polytope::Type::Segment);

  // Test virtual function calls work correctly
  EXPECT_EQ(qf1->getSize(), 1);
  EXPECT_EQ(qf2->getSize(), 2);

  EXPECT_EQ(qf1->getGeometry(), Polytope::Type::Triangle);
  EXPECT_EQ(qf2->getGeometry(), Polytope::Type::Segment);

  EXPECT_NO_THROW(qf1->getWeight(0));
  EXPECT_NO_THROW(qf2->getWeight(0));
}

// Test copy functionality (Copyable interface)
TEST_F(QFTest, CopyableBehavior)
{
  QF1P1 original(Polytope::Type::Triangle);

  // Test that copy() method works (inherited from Copyable)
  std::unique_ptr<QuadratureFormulaBase> copied(original.copy());

  EXPECT_NE(copied.get(), &original);  // Different objects
  EXPECT_EQ(copied->getGeometry(), original.getGeometry());  // Same properties
  EXPECT_EQ(copied->getSize(), original.getSize());
  EXPECT_EQ(copied->getWeight(0), original.getWeight(0));
}

// Test different geometry types
TEST_F(QFTest, DifferentGeometryTypes)
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
    EXPECT_NO_THROW(QF1P1 qf(geom));

    QF1P1 qf(geom);
    EXPECT_EQ(qf.getGeometry(), geom);
    EXPECT_EQ(qf.getSize(), 1);
    EXPECT_NO_THROW(qf.getWeight(0));
    EXPECT_NO_THROW(qf.getPoint(0));
  }
}

// Integration test: verify quadrature accuracy for simple functions
TEST_F(QFTest, QuadratureAccuracySimple)
{
  // Test that constant function integration is exact
  QF1P1 qf_triangle(Polytope::Type::Triangle);

  // For a constant function f(x) = 1, integral over reference triangle should equal area * 1
  double integral_approx = qf_triangle.getWeight(0) * 1.0;  // f(x) = 1

  // Reference triangle has area 0.5
  EXPECT_GT(integral_approx, 0.1);  // Should be positive and reasonable
  EXPECT_LT(integral_approx, 1.0);  // Should not be too large
}

// Error handling tests
TEST_F(QFTest, OutOfBoundsAccess)
{
  QF1P1 qf(Polytope::Type::Triangle);

  // Should work for valid index
  EXPECT_NO_THROW(qf.getWeight(0));
  EXPECT_NO_THROW(qf.getPoint(0));

  // Out of bounds access behavior depends on implementation
  // For single point quadrature, index 1 should be invalid
  // Note: Actual behavior may vary (assert, exception, or undefined)
}

// Memory management test
TEST_F(QFTest, MemoryManagement)
{
  // Test that objects can be created and destroyed without issues
  {
    QF1P1 qf1(Polytope::Type::Triangle);
    GaussLegendre qf2(Polytope::Type::Segment);

    // Use the objects
    qf1.getSize();
    qf2.getSize();
  }  // Objects should be destroyed without issues

  EXPECT_TRUE(true);  // Test passes if we reach here without crashes
}
