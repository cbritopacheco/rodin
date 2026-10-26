/*
 * Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE or https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief High-order H1 basis-lifecycle registrations.
 */
#include "BasisLifecycle.h"

namespace Rodin::Tests::Benchmarks
{
  /// @brief Registers degree-five and degree-six lifecycle protocols.
  void registerH1BasisLifecycleHigh()
  {
    registerH1BasisLifecycleDegree<5>();
    registerH1BasisLifecycleDegree<6>();
  }
}
