/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <gtest/gtest.h>
#include <limits>
#include "Rodin/Adaptation/WNGIR/DirectionalNewton.h"
#include "Rodin/Adaptation/WNGIR/Loss.h"

using namespace Rodin;

TEST(WNGIRDirectionalNewton, StepSelectionAndFallback)
{
  using Adaptation::wngirDirectionalNewtonStep;
  EXPECT_DOUBLE_EQ(wngirDirectionalNewtonStep(4, 2, 3, 1, 100), 2);
  EXPECT_DOUBLE_EQ(wngirDirectionalNewtonStep(1, 2, 3, 1, 100), 0.5);
  EXPECT_DOUBLE_EQ(wngirDirectionalNewtonStep(400, 2, 3, 1, 100), 100);
  EXPECT_DOUBLE_EQ(wngirDirectionalNewtonStep(1, 0, 2, 1, 100), 0.5);
  EXPECT_DOUBLE_EQ(wngirDirectionalNewtonStep(1, -1, 2, 1, 100), 0.5);
  EXPECT_DOUBLE_EQ(wngirDirectionalNewtonStep(1e-9, 1, 1, 1, 100), 1e-9);
  EXPECT_DOUBLE_EQ(wngirDirectionalNewtonStep(-1, 1, 1, 1, 100), 0);
  EXPECT_DOUBLE_EQ(
    wngirDirectionalNewtonStep(1, std::numeric_limits<Real>::infinity(), 2, 1, 100), 0.5);
  EXPECT_DOUBLE_EQ(wngirDirectionalNewtonStep(1, 0, 0, 1, 100), 0);
}

TEST(WNGIRDirectionalNewton, PhysicalStepIsIndependentOfDirectionScale)
{
  using Adaptation::wngirDirectionalNewtonStep;
  for (const Real curvature : {Real(-2), Real(2)})
    for (const Real bound : {Real(0), Real(0.1), Real(100)})
    {
      const Real reference = wngirDirectionalNewtonStep(1, curvature, 3, 2, bound);
      for (const Real scale : {Real(1e-8), Real(1e-3), Real(1e4)})
        EXPECT_NEAR(scale *
            wngirDirectionalNewtonStep(
              scale, scale * scale * curvature, scale * scale * 3, scale * 2, bound),
          reference, 1e-14);
    }
}

TEST(WNGIRDirectionalNewton, ZeroMotionBoundIsUnrestricted)
{
  using Adaptation::wngirDirectionalNewtonStep;
  EXPECT_DOUBLE_EQ(wngirDirectionalNewtonStep(400, 2, 3, 1, 0), 200);
  EXPECT_DOUBLE_EQ(wngirDirectionalNewtonStep(400, 2, 3, 1e6, 0), 200);
  EXPECT_DOUBLE_EQ(wngirDirectionalNewtonStep(400, -1, 4, 1e6, 0), 100);
  EXPECT_DOUBLE_EQ(wngirDirectionalNewtonStep(1, 2, 3, 1, -1), 0);
  EXPECT_DOUBLE_EQ(
    wngirDirectionalNewtonStep(1, 2, 3, 1, std::numeric_limits<Real>::quiet_NaN()), 0);
}

TEST(WNGIRDirectionalNewton, AffineResidualCurvatureMatchesForceDifference)
{
  const Adaptation::WNGIRLoss loss(Real(0.3));
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
