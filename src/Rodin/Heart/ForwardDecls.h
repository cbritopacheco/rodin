/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ForwardDecls.h
 * @brief Forward declarations for the Heart module.
 */
#ifndef RODIN_HEART_FORWARDDECLS_H
#define RODIN_HEART_FORWARDDECLS_H

#include "Rodin/Types.h"
#include "Rodin/Heart/CCMLC2014/HolzapfelReducedLaw.h"
#include "Rodin/Heart/CCMLC2014/PassiveLaw.h"

namespace Rodin::Heart
{
  namespace CCMLC2014::Solver
  {
    template <class PassiveEnergyLaw, class PassiveLaw>
    class StepperT;
  }

  /**
   * @brief Default CCMLC2014 0D ventricular model type.
   *
   * The default template arguments are declared here, and only here — an
   * alias template cannot be redeclared elsewhere with defaults.
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
