/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief P1 and constant-space basis-production registrations.
 */
#include "Mult.h"

namespace Rodin::Tests::Benchmarks
{
  void registerBasisScalingFixedDegree()
  {
    registerBasisScaling<1, BasisFamily::P1>();
    registerBasisScaling<0, BasisFamily::P0>();
    registerBasisScaling<0, BasisFamily::P0g>();
  }
}
