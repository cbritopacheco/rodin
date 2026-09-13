/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Tests for convergence-study infrastructure shared by all strategies.
 */

#include <cmath>

#include <gtest/gtest.h>

#include "Convergence.h"

namespace Rodin::Tests::Convergence
{
  /**
   * @brief Verifies that algebraic rates use the actual scale ratio.
   *
   * Here @f$h@f$ decreases by four rather than two. Errors proportional to
   * @f$h^2@f$ and @f$h@f$ must still report rates two and one respectively.
   */
  TEST(ErrorHistoryTest, ComputesAlgebraicRatesForNonDyadicSpacing)
  {
    ErrorHistory history;
    history.append(0.5, ErrorNorms(0.25, 0.5))
           .append(0.125, ErrorNorms(0.015625, 0.125));

    const auto rates = history.getAlgebraicRates(1);
    EXPECT_DOUBLE_EQ(rates.getL2(), 2);
    EXPECT_DOUBLE_EQ(rates.getH1Seminorm(), 1);
  }

  /**
   * @brief Verifies exponential rates for non-unit degree spacing.
   *
   * Errors @f$e^{-2p}@f$ and @f$e^{-p}@f$ must report decay constants two
   * and one even when two polynomial degrees separate the samples.
   */
  TEST(ErrorHistoryTest, ComputesExponentialRatesForNonUnitSpacing)
  {
    ErrorHistory history;
    history.append(1, ErrorNorms(std::exp(-2), std::exp(-1)))
           .append(3, ErrorNorms(std::exp(-6), std::exp(-3)));

    const auto rates = history.getExponentialRates(1);
    EXPECT_DOUBLE_EQ(rates.getL2(), 2);
    EXPECT_DOUBLE_EQ(rates.getH1Seminorm(), 1);
  }
}
