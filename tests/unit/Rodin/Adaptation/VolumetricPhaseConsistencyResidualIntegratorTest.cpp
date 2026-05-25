/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "TMOPFDTestHelpers.h"

TEST(Rodin_Adaptation_TargetMatrixOptimization_PerClass, VolumetricPhaseConsistencyResidualIntegratorFiniteDifferenceConsistencyP1)
{
  const auto result = Rodin::Tests::Unit::TMOPFD::phaseFdSweep(false);
  Rodin::Tests::Unit::TMOPFD::printFDSweep("VolumetricPhaseConsistencyResidualIntegrator", "P1", result);
  EXPECT_LT(result.bestError(), Rodin::Real(3e-5));
  EXPECT_TRUE(Rodin::Tests::Unit::TMOPFD::hasEpsilonScalingTrend(result));
}

TEST(Rodin_Adaptation_TargetMatrixOptimization_PerClass, VolumetricPhaseConsistencyResidualIntegratorFiniteDifferenceConsistencyP2)
{
  const auto result = Rodin::Tests::Unit::TMOPFD::phaseFdSweep(true);
  Rodin::Tests::Unit::TMOPFD::printFDSweep("VolumetricPhaseConsistencyResidualIntegrator", "P2", result);
  EXPECT_LT(result.bestError(), Rodin::Real(5e-5));
  EXPECT_TRUE(Rodin::Tests::Unit::TMOPFD::hasEpsilonScalingTrend(result));
}
