/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief SWIFT termination tokens remain distinct and usable after report destruction. */
#include <gtest/gtest.h>
#include <array>
#include <limits>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include "Rodin/Adaptation/SWIFT/Report.h"

namespace Rodin::Tests::Unit
{
  /// @brief Log consumers can distinguish every outcome and retain its static token.
  TEST(Adaptation_SWIFT_Report, TokensAreDistinctAndOutliveReports)
  {
    using Report = Adaptation::SWIFT::Report;
    using Reason = Report::Reason;
    constexpr std::array Reasons{Reason::IterationLimit, Reason::MissingInterface,
      Reason::EmptyInterface, Reason::DegenerateGradient, Reason::InvalidInitialGeometry,
      Reason::InvalidGeometry, Reason::GeometricTarget, Reason::PredictorFailure,
      Reason::NonDescentPredictor, Reason::InvalidScaling, Reason::LinearFailure,
      Reason::InnerLineSearchFailure, Reason::InnerIterationLimit,
      Reason::NonfiniteDirection, Reason::NonDescentDirection, Reason::LineSearchFailure,
      Reason::SmallAcceptedSteps, Reason::SmallEnergyChanges};
    std::set<std::string> distinct;
    std::vector<std::pair<const char*, std::string>> retained;
    for (const Reason reason : Reasons)
    {
      SCOPED_TRACE(static_cast<int>(reason));
      Report report;
      report.reason = reason;
      const char* token = report.getReasonString();
      ASSERT_NE(token, nullptr);
      const std::string value(token);
      ASSERT_FALSE(value.empty());
      EXPECT_NE(value, "unknown");
      EXPECT_TRUE(distinct.insert(value).second);
      for (const char character : value)
      {
        EXPECT_TRUE((character >= 'a' && character <= 'z') || character == '-' ||
          (character >= '0' && character <= '9'));
      }
      retained.emplace_back(token, value);
    }
    for (const auto& [token, value] : retained)
      EXPECT_EQ(std::string(token), value);
    EXPECT_EQ(distinct.size(), Reasons.size());
    EXPECT_STREQ(Report().getReasonString(), "iter-budget");

    // Values from an unrecognized/newer producer still have a safe fallback.
    for (const int unknown :
      {-1, static_cast<int>(Reasons.size()), std::numeric_limits<int>::max()})
    {
      Report report;
      report.reason = static_cast<Reason>(unknown);
      EXPECT_STREQ(report.getReasonString(), "unknown");
    }
  }
}
