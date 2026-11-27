/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <gtest/gtest.h>
#include <limits>
#include "Rodin/Adaptation/SWIFT/DirectionalNewton.h"
#include "Rodin/Adaptation/SWIFT/Loss.h"

using namespace Rodin;
namespace SWIFT = Rodin::Adaptation::SWIFT;

TEST(SWIFTDirectionalNewton, NonpositiveCurvatureUsesMotionBound)
{
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(1, 0, 4, 8), 2);
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(1, -1, 4, 8), 2);
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(1, 0, 4, 0), 0);
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(1, -1, 4, 0), 0);
}

TEST(SWIFTDirectionalNewton, PositiveCurvatureSelectsNewtonStep)
{
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(4, 2, 1, 100), 2);
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(1, 2, 1, 100), 0.5);
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(400, 2, 1, 100), 100);
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(4, 2, 4, 2), 0.5);
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(1e-9, 1, 1, 100), 1e-9);
}

TEST(SWIFTDirectionalNewton, BoundedStepMinimizesScalarModel)
{
  for (const Real action : {Real(0.1), Real(1), Real(10)})
  {
    for (const Real curvature : {Real(-2), Real(0), Real(0.5), Real(4)})
    {
      const Real norm = Real(2);
      const Real bound = Real(3);
      const Real step = SWIFT::getDirectionalNewtonStep(action, curvature, norm, bound);
      ASSERT_GT(step, Real(0));
      ASSERT_LE(norm * step, bound);
      const Real value = -action * step + Real(0.5) * curvature * step * step;
      for (Integer i = 0; i <= 100; ++i)
      {
        const Real candidate = bound / norm * i / Real(100);
        const Real candidateValue =
          -action * candidate + Real(0.5) * curvature * candidate * candidate;
        EXPECT_LE(value, candidateValue + Real(1e-14));
      }
    }
  }
}

TEST(SWIFTDirectionalNewton, RejectsInvalidInputs)
{
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(-1, 1, 1, 100), 0);
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(0, 1, 1, 100), 0);
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(1, 1, 0, 100), 0);
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(1, 1, -1, 100), 0);
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(1, 1, 1, -1), 0);
  for (const Real value : {std::numeric_limits<Real>::infinity(),
         -std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()})
  {
    EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(value, 1, 1, 100), 0);
    EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(1, value, 1, 100), 0);
    EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(1, 1, value, 100), 0);
    EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(1, 1, 1, value), 0);
  }
}

TEST(SWIFTDirectionalNewton, PhysicalStepIsIndependentOfDirectionScale)
{
  for (const Real curvature : {Real(-2), Real(0), Real(2)})
  {
    for (const Real bound : {Real(0), Real(0.1), Real(100)})
    {
      const Real reference = SWIFT::getDirectionalNewtonStep(1, curvature, 2, bound);
      for (const Real scale : {Real(1e-8), Real(1e-3), Real(1e4)})
      {
        EXPECT_NEAR(scale *
            SWIFT::getDirectionalNewtonStep(
              scale, scale * scale * curvature, scale * 2, bound),
          reference, 1e-14);
      }
    }
  }
}

TEST(SWIFTDirectionalNewton, ZeroMotionBoundIsUnrestricted)
{
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(400, 2, 1, 0), 200);
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(400, 2, 1e6, 0), 200);
}

TEST(SWIFTDirectionalNewton, AffineResidualCurvatureMatchesForceDifference)
{
  const Adaptation::SWIFT::Loss loss(Real(0.3));
  for (Real residual : {Real(0), Real(0.1), Real(0.3), Real(-0.7)})
  {
    for (Real projection : {Real(-2), Real(0.4)})
    {
      const Real epsilon = Real(1e-6);
      const Real curvature = loss.getCurvature(residual) * projection * projection;
      const Real difference = (loss.getInfluence(residual + epsilon * projection) -
                                loss.getInfluence(residual - epsilon * projection)) *
        projection / (Real(2) * epsilon);
      EXPECT_NEAR(curvature, difference, 1e-8);
    }
  }
}
