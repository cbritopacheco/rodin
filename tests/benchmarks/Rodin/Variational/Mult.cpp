/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Public scalar-weighted basis production benchmark entry point.
 */
#include "Mult.h"

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  registerBasisScalingH1();
  registerBasisScalingHigherOrder();
  registerBasisScalingFixedDegree();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  benchmark::AddCustomContext("scope",
    "Warmed basis-value production; no pair contraction, integration or assembly");
  const auto cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return !cases || basisScalingFailed ? 1 : 0;
}
