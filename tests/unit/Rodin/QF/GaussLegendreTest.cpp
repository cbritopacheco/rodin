/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/QF/GaussLegendre.h"
#include "Rodin/Geometry.h"

using namespace Rodin::QF;
using namespace Rodin::Geometry;

// Test class for GaussLegendre functionality
class GaussLegendreTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// Test GaussLegendre quadrature
TEST_F(GaussLegendreTest, BasicProperties)
{
  GaussLegendre gl_segment(Polytope::Type::Segment);

  // GaussLegendre with 2 points
  EXPECT_EQ(gl_segment.getSize(), 2);
  EXPECT_EQ(gl_segment.getGeometry(), Polytope::Type::Segment);
}

TEST_F(GaussLegendreTest, Weights)
{
  GaussLegendre gl_segment(Polytope::Type::Segment);

  // Should have weights for each quadrature point
  for (size_t i = 0; i < gl_segment.getSize(); ++i) {
    EXPECT_NO_THROW(gl_segment.getWeight(i));
    EXPECT_GT(gl_segment.getWeight(i), 0.0);
  }
}

TEST_F(GaussLegendreTest, Points)
{
  GaussLegendre gl_segment(Polytope::Type::Segment);

  // Test points are accessible (though implementation may have limitations)
  EXPECT_NO_THROW(gl_segment.getPoint(0));

  const auto& point = gl_segment.getPoint(0);
  EXPECT_GE(point.size(), 1);
}