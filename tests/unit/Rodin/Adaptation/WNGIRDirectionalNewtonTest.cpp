/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <gtest/gtest.h>
#include <limits>
#include "Rodin/Adaptation/WNGIRDirectionalNewton.h"
#include "Rodin/Adaptation/WNGIRLoss.h"

using namespace Rodin;

TEST(WNGIRDirectionalNewton, StepSelectionAndFallback)
{
  using Adaptation::Detail::wngirDirectionalNewtonStep;
  EXPECT_DOUBLE_EQ(wngirDirectionalNewtonStep(4, 2, 1e-6, 100), 2);
  EXPECT_DOUBLE_EQ(wngirDirectionalNewtonStep(1, 2, 1e-6, 100), 0.5);
  EXPECT_DOUBLE_EQ(wngirDirectionalNewtonStep(400, 2, 1e-6, 100), 100);
  EXPECT_DOUBLE_EQ(wngirDirectionalNewtonStep(1, 0, 1e-6, 100), 1);
  EXPECT_DOUBLE_EQ(wngirDirectionalNewtonStep(1, -1, 1e-6, 100), 1);
  EXPECT_DOUBLE_EQ(wngirDirectionalNewtonStep(1e-9, 1, 1e-6, 100), 1);
  EXPECT_DOUBLE_EQ(wngirDirectionalNewtonStep(-1, 1, 1e-6, 100), 1);
  EXPECT_DOUBLE_EQ(wngirDirectionalNewtonStep(1, std::numeric_limits<Real>::infinity(), 1e-6, 100), 1);
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
