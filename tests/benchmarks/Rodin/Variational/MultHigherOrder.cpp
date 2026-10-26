/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Degree-four through degree-six H1 basis-production registrations.
 */
#include "Mult.h"

namespace Rodin::Tests::Benchmarks
{
  void registerBasisScalingHigherOrder()
  {
    registerBasisScaling<4, BasisFamily::H1>();
    registerBasisScaling<5, BasisFamily::H1>();
    registerBasisScaling<6, BasisFamily::H1>();
  }
}
