/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_HEART_FORWARDDECLS_H
#define RODIN_HEART_FORWARDDECLS_H

#include "Rodin/Types.h"
#include "Rodin/Heart/CCMLC2014/HolzapfelReducedLaw.h"
#include "Rodin/Heart/CCMLC2014/PassiveLaw.h"

namespace Rodin::Heart
{
  template <class PassiveEnergyLaw, class PassiveLaw>
  class CCMLC2014T;

  using CCMLC2014 = CCMLC2014T<HolzapfelReducedLaw<Real>, CCMLC2014PassiveLaw<Real>>;
}

#endif
