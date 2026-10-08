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

TEST(SWIFTDirectionalNewton, StepSelectionAndFallback)
{
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(4, 2, 3, 1, 100), 2);
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(1, 2, 3, 1, 100), 0.5);
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(400, 2, 3, 1, 100), 100);
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(1, 0, 2, 1, 100), 0.5);
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(1, -1, 2, 1, 100), 0.5);
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(1e-9, 1, 1, 1, 100), 1e-9);
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(-1, 1, 1, 1, 100), 0);
  EXPECT_DOUBLE_EQ(
    SWIFT::getDirectionalNewtonStep(1, std::numeric_limits<Real>::infinity(), 2, 1, 100), 0.5);
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(1, 0, 0, 1, 100), 0);
}

TEST(SWIFTDirectionalNewton, PhysicalStepIsIndependentOfDirectionScale)
{
  for (const Real curvature : {Real(-2), Real(2)})
    for (const Real bound : {Real(0), Real(0.1), Real(100)})
    {
      const Real reference = SWIFT::getDirectionalNewtonStep(1, curvature, 3, 2, bound);
      for (const Real scale : {Real(1e-8), Real(1e-3), Real(1e4)})
        EXPECT_NEAR(scale *
            SWIFT::getDirectionalNewtonStep(
              scale, scale * scale * curvature, scale * scale * 3, scale * 2, bound),
          reference, 1e-14);
    }
}

TEST(SWIFTDirectionalNewton, ZeroMotionBoundIsUnrestricted)
{
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(400, 2, 3, 1, 0), 200);
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(400, 2, 3, 1e6, 0), 200);
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(400, -1, 4, 1e6, 0), 100);
  EXPECT_DOUBLE_EQ(SWIFT::getDirectionalNewtonStep(1, 2, 3, 1, -1), 0);
  EXPECT_DOUBLE_EQ(
    SWIFT::getDirectionalNewtonStep(1, 2, 3, 1, std::numeric_limits<Real>::quiet_NaN()), 0);
}

TEST(SWIFTDirectionalNewton, AffineResidualCurvatureMatchesForceDifference)
{
  const Adaptation::SWIFT::Loss loss(Real(0.3));
  for (Real residual : {Real(0), Real(0.1), Real(0.3), Real(-0.7)})
    for (Real projection : {Real(-2), Real(0.4)})
    {
      const Real epsilon = Real(1e-6);
      const Real curvature = loss.getWeight(residual) *
        (Real(1) - Real(2) * residual * residual / loss.getScaleSquared()) * projection * projection;
      const Real difference = (loss.getInfluence(residual + epsilon * projection) -
        loss.getInfluence(residual - epsilon * projection)) * projection / (Real(2) * epsilon);
      EXPECT_NEAR(curvature, difference, 1e-8);
    }
}
