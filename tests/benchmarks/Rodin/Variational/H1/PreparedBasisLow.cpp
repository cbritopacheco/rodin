/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Registers low-order prepared-reference basis consumers.
 */
#include "../Mult.h"
#include "PreparedBasis.h"

namespace Rodin::Tests::Benchmarks
{
  void registerH1PreparedBasisLow()
  {
    registerBasisScaling<1, BasisFamily::H1, true>();
    registerBasisScaling<2, BasisFamily::H1, true>();
    registerBasisScaling<3, BasisFamily::H1, true>();
  }
}
