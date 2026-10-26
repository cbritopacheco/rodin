/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Measures H1 basis consumption with certified reference-table readiness.
 *
 * Preparation is outside the measured region. Every selected case checks
 * reference-table presence before and after timing. Binding-primed controls
 * remain in RodinBasisScalingBenchmarks and are not relabelled as prepared.
 */
#include "../Mult.h"
#include "PreparedBasis.h"

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  registerH1PreparedBasisLow();
  registerH1PreparedBasisHigh();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  benchmark::AddCustomContext("scope",
    "Prepared-reference public basis consumption; preparation, pair contraction and "
    "assembly excluded");
  const auto cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return !cases || basisScalingFailed ? 1 : 0;
}
