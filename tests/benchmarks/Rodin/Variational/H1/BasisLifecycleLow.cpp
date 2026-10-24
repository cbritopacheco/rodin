/*
 * Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE or https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Low-order H1 basis-lifecycle registrations.
 */
#include "BasisLifecycle.h"

namespace Rodin::Tests::Benchmarks
{
  /// @brief Registers degree-one and degree-two lifecycle protocols.
  void registerH1BasisLifecycleLow()
  {
    registerH1BasisLifecycleDegree<1>();
    registerH1BasisLifecycleDegree<2>();
  }
}
