#include <gtest/gtest.h>

#include "Rodin/Adaptation/WNGIRLoss.h"

namespace Rodin::Adaptation
{
  /// @brief Verifies that the analytic influence equals the loss derivative.
  TEST(Rodin_Adaptation_WNGIRLoss, InfluenceMatchesFiniteDifference)
  {
    const WNGIRLoss loss(Real(0.3));
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
  TEST(Rodin_Adaptation_WNGIRLoss, SaturatesAtHalfScaleSquared)
  {
    const WNGIRLoss loss(Real(0.3));
    EXPECT_NEAR(loss.getValue(Real(100)), Real(0.045), Real(1e-12));
    EXPECT_LT(loss.getWeight(Real(100)), Real(1e-12));
  }
}
