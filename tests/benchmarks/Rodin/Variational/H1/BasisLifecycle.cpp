/*
 * Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE or https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief H1 first-binding, identity-working-set and sample-sweep driver.
 */
#include "BasisLifecycle.h"

namespace Rodin::Tests::Benchmarks
{
  void registerH1BasisLifecycleLow();
  void registerH1BasisLifecycleMiddle();
  void registerH1BasisLifecycleHigh();
}

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  registerH1BasisLifecycleLow();
  registerH1BasisLifecycleMiddle();
  registerH1BasisLifecycleHigh();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return h1BasisLifecycleFailed ? 1 : 0;
}
