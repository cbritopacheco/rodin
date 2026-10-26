/*
 * Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE or https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Intermediate-order H1 basis-lifecycle registrations.
 */
#include "BasisLifecycle.h"

namespace Rodin::Tests::Benchmarks
{
  /// @brief Registers degree-three and degree-four lifecycle protocols.
  void registerH1BasisLifecycleMiddle()
  {
    registerH1BasisLifecycleDegree<3>();
    registerH1BasisLifecycleDegree<4>();
  }
}
