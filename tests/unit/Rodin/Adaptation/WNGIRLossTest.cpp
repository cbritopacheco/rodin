#include <gtest/gtest.h>

#include "Rodin/Adaptation/WNGIRLoss.h"

namespace Rodin::Adaptation
{
  TEST(Rodin_Adaptation_WNGIRLoss, InfluenceMatchesFiniteDifference)
  {
    for (const auto type :
      {WNGIRLossType::Welsch, WNGIRLossType::Cauchy, WNGIRLossType::PseudoHuber})
    {
      const WNGIRLoss loss(type, Real(0.3));
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
  }

  TEST(Rodin_Adaptation_WNGIRLoss, DistinguishesBoundedAndUnboundedLosses)
  {
    const WNGIRLoss welsch(WNGIRLossType::Welsch, Real(0.3));
    const WNGIRLoss cauchy(WNGIRLossType::Cauchy, Real(0.3));
    EXPECT_NEAR(welsch.getValue(Real(100)), Real(0.045), Real(1e-12));
    EXPECT_GT(cauchy.getValue(Real(100)), cauchy.getValue(Real(10)));
  }
}
