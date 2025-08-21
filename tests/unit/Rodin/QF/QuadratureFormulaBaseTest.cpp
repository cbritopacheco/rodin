/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <memory>

#include "Rodin/QF/QuadratureFormula.h"
#include "Rodin/QF/QF1P1.h"
#include "Rodin/QF/GaussLegendre.h"
#include "Rodin/Geometry.h"

using namespace Rodin::QF;
using namespace Rodin::Geometry;

// Test class for QuadratureFormulaBase functionality
class QuadratureFormulaBaseTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// Test QuadratureFormulaBase construction and basic functionality
TEST_F(QuadratureFormulaBaseTest, Construction)
{
  // Note: Can't directly test abstract base class, but we can test through derived classes
  QF1P1 qf1p1(Polytope::Type::Triangle);

  EXPECT_EQ(qf1p1.getGeometry(), Polytope::Type::Triangle);

  // Test copy constructor
  QF1P1 qf1p1_copy(qf1p1);
  EXPECT_EQ(qf1p1_copy.getGeometry(), Polytope::Type::Triangle);
}

// Test polymorphic behavior
TEST_F(QuadratureFormulaBaseTest, PolymorphicBehavior)
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
TEST_F(QuadratureFormulaBaseTest, CopyableBehavior)
{
  QF1P1 original(Polytope::Type::Triangle);

  // Test that copy() method works (inherited from Copyable)
  std::unique_ptr<QuadratureFormulaBase> copied(original.copy());

  EXPECT_NE(copied.get(), &original);  // Different objects
  EXPECT_EQ(copied->getGeometry(), original.getGeometry());  // Same properties
  EXPECT_EQ(copied->getSize(), original.getSize());
  EXPECT_EQ(copied->getWeight(0), original.getWeight(0));
}

// Memory management test
TEST_F(QuadratureFormulaBaseTest, MemoryManagement)
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