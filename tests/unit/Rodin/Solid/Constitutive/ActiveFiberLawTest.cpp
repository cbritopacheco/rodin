/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/// @file @brief Starling-envelope properties and steady active-state invariants.
#include <gtest/gtest.h>
#include <cmath>
#include <utility>

#include "Rodin/Solid/Constitutive/ActiveFiberLaw.h"

using namespace Rodin;
using namespace Rodin::Solid;

/// @brief The calibrated activation envelope is bounded, continuous and has a unit plateau.
TEST(Solid_ActiveFiberLaw, StarlingEnvelopeAndContinuity)
{
  constexpr Real step = Real(0.01);
  Real previous = ActiveFiberLaw::starling(Real(-0.6));
  for (Index i = 0; i <= 320; ++i)
  {
    const Real extension = Real(-0.6) + Real(i) * step;
    const Real gain = ActiveFiberLaw::starling(extension);
    SCOPED_TRACE(extension);
    ASSERT_TRUE(std::isfinite(gain));
    EXPECT_GE(gain, Real(0));
    EXPECT_LE(gain, Real(1));
    // Every calibrated segment has slope magnitude below one.
    EXPECT_LE(std::abs(gain - previous), step + Real(1e-14));
    if (extension < Real(-0.4) || extension > Real(2.4))
    {
      EXPECT_EQ(gain, Real(0));
    }
    if (extension >= Real(1) && extension <= Real(1.3))
    {
      EXPECT_NEAR(gain, Real(1), Real(1e-14));
    }
    previous = gain;
  }
}

using GainPoint = std::pair<Real, Real>;
class Solid_ActiveFiberEquilibrium : public testing::TestWithParam<GainPoint>
{};

/// @brief Positive constant activation preserves the analytic cross-bridge equilibrium.
TEST_P(Solid_ActiveFiberEquilibrium, ConstantExtensionPreservesSteadyState)
{
  const auto [extension, gain] = GetParam();
  constexpr Real crossBridgeStiffness = Real(7.5), contractility = Real(3.25);
  ActiveFiberLaw::Parameters parameters;
  parameters.crossBridgeStiffness = crossBridgeStiffness;
  parameters.contractility = contractility;
  parameters.destructionRate = Real(2.5);
  const ActiveFiberLaw law(parameters);
  const ActiveFiberLaw::State equilibrium{std::sqrt(gain * crossBridgeStiffness),
    gain * contractility / std::sqrt(gain * crossBridgeStiffness)};
  for (Real activation : {Real(0.2), Real(1.3)})
  {
    for (Real dt : {Real(0), Real(0.001), Real(0.5), Real(2)})
    {
      SCOPED_TRACE(::testing::Message() << "activation=" << activation << " dt=" << dt);
      const auto updated = law.update(dt, equilibrium, extension, extension, activation);
      ASSERT_TRUE(std::isfinite(updated.gamma));
      ASSERT_TRUE(std::isfinite(updated.beta));
      EXPECT_NEAR(updated.gamma, equilibrium.gamma, Real(1e-12));
      EXPECT_NEAR(updated.beta, equilibrium.beta, Real(1e-12));
      EXPECT_NEAR(
        updated.gamma * updated.gamma, gain * crossBridgeStiffness, Real(1e-12));
      EXPECT_NEAR(updated.activeStress(), gain * contractility, Real(1e-12));
    }
  }
}

// Interior calibration points exercise each linear segment and the plateau.
INSTANTIATE_TEST_SUITE_P(StarlingSegments, Solid_ActiveFiberEquilibrium,
  testing::Values(GainPoint{Real(-0.05), Real(0.19)}, GainPoint{Real(0.515), Real(0.56)},
    GainPoint{Real(0.865), Real(0.87)}, GainPoint{Real(1.15), Real(1)},
    GainPoint{Real(1.85), Real(0.5)}));
