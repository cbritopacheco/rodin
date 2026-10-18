/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Registers high-order prepared-reference basis consumers.
 */
#include "../Mult.h"
#include "PreparedBasis.h"

namespace Rodin::Tests::Benchmarks
{
  void registerH1PreparedBasisHigh()
  {
    registerBasisScaling<4, BasisFamily::H1, true>();
    registerBasisScaling<5, BasisFamily::H1, true>();
    registerBasisScaling<6, BasisFamily::H1, true>();
  }
}
