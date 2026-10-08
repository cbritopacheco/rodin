/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <gtest/gtest.h>

#include "Rodin/Adaptation/SWIFT/Loss.h"

namespace Rodin::Adaptation
{
  TEST(Rodin_Adaptation_SWIFTLoss, SmallWelschResidualDoesNotCancel)
  {
    const SWIFT::Loss loss(1);
    EXPECT_NEAR(loss.getValue(Real(1e-10)), Real(0.5e-20), Real(1e-35));
  }

  TEST(Rodin_Adaptation_SWIFTLoss, DroppedLevelSetHessianCurvature)
  {
    const SWIFT::Loss loss(Real(0.3));
    constexpr Real eps = Real(1e-6);
    for (const Real r : {Real(0), Real(0.1), Real(0.3), Real(-0.7)})
    {
      EXPECT_NEAR(loss.getCurvature(r),
        (loss.getInfluence(r + eps) - loss.getInfluence(r - eps)) / (Real(2) * eps),
        Real(1e-9));
    }
  }
  /// @brief Verifies that the analytic influence equals the loss derivative.
  TEST(Rodin_Adaptation_SWIFTLoss, InfluenceMatchesFiniteDifference)
  {
    const SWIFT::Loss loss(Real(0.3));
    for (const Real residual : {Real(-0.7), Real(-0.1), Real(0.1), Real(0.7)})
    {
      const Real epsilon = Real(1e-7);
      const Real finiteDifference =
        (loss.getValue(residual + epsilon) - loss.getValue(residual - epsilon)) /
        (Real(2) * epsilon);
      EXPECT_NEAR(loss.getInfluence(residual), finiteDifference, Real(1e-9));
      EXPECT_GT(loss.getWeight(residual), Real(0));
    }
  }

  /// @brief Verifies bounded Welsch energy.
  TEST(Rodin_Adaptation_SWIFTLoss, SaturatesAtHalfScaleSquared)
  {
    const SWIFT::Loss loss(Real(0.3));
    EXPECT_NEAR(loss.getValue(Real(100)), Real(0.045), Real(1e-12));
    EXPECT_LT(loss.getWeight(Real(100)), Real(1e-12));
  }
}
