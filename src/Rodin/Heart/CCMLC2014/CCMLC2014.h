/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file CCMLC2014.h
 * @brief User-facing alias for the CCMLC2014 0D ventricular model stepper.
 */
#ifndef RODIN_HEART_CCMLC2014_CCMLC2014_H
#define RODIN_HEART_CCMLC2014_CCMLC2014_H

#include "Rodin/Heart/CCMLC2014/HolzapfelReducedLaw.h"
#include "Rodin/Heart/CCMLC2014/PassiveLaw.h"
#include "Rodin/Heart/CCMLC2014/Solver/Stepper.h"

namespace Rodin::Heart
{
  /**
   * @brief Default CCMLC2014 0D model type.
   *
   * @tparam PassiveEnergyLaw Reduced passive energy law.
   * @tparam PassiveLaw Passive stress operator from the reduced invariants.
   */
  template <
    class PassiveEnergyLaw = HolzapfelReducedLaw<Real>,
    class PassiveLaw = CCMLC2014PassiveLaw<Real>>
  using CCMLC2014T = CCMLC2014::Solver::StepperT<PassiveEnergyLaw, PassiveLaw>;
}

#endif
