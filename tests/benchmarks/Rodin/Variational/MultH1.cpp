/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Degree-one through degree-three H1 basis-production registrations.
 */
#include "Mult.h"

namespace Rodin::Tests::Benchmarks
{
  void registerBasisScalingH1()
  {
    registerBasisScaling<1, BasisFamily::H1>();
    registerBasisScaling<2, BasisFamily::H1>();
    registerBasisScaling<3, BasisFamily::H1>();
  }
}
