/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_LEVELSETREGISTRATION_H
#define RODIN_ADAPTATION_LEVELSETREGISTRATION_H

#include "LSRRegistration.h"

namespace Rodin::Adaptation
{
  template <class Trial, class Test, class State>
  class LevelSetRegistration
    : public LSRRegistration<Trial, Test, State>
  {
    public:
      using LSRRegistration<Trial, Test, State>::LSRRegistration;
  };

  template <class Trial, class Test, class State>
  LevelSetRegistration(const Trial&, const Test&, const State&)
    -> LevelSetRegistration<Trial, Test, State>;
}

#endif
