/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_SIGNEDDISTANCEREGISTRATION_H
#define RODIN_ADAPTATION_SIGNEDDISTANCEREGISTRATION_H

#include "SDFRRegistration.h"

namespace Rodin::Adaptation
{
  template <class Trial, class Test, class State>
  class SignedDistanceRegistration
    : public SDFRRegistration<Trial, Test, State>
  {
    public:
      using SDFRRegistration<Trial, Test, State>::SDFRRegistration;
  };

  template <class Trial, class Test, class State>
  SignedDistanceRegistration(const Trial&, const Test&, const State&)
    -> SignedDistanceRegistration<Trial, Test, State>;
}

#endif
