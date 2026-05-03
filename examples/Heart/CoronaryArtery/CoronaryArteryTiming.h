// CoronaryArteryTiming.h
#ifndef EXAMPLES_HEART_CORONARYARTERY_CORONARYARTERYTIMING_H
#define EXAMPLES_HEART_CORONARYARTERY_CORONARYARTERYTIMING_H

#include <algorithm>
#include <chrono>

#include <Rodin/Types.h>

namespace Rodin::Examples::Heart
{
  using CoronaryClock = std::chrono::steady_clock;

  struct StepTiming
  {
    Real total = 0.0;
    Real advance0D = 0.0;
    Real setup3DForm = 0.0;
    Real assemble3D = 0.0;
    Real fieldSplits = 0.0;
    Real solve3D = 0.0;
    Real fluxes = 0.0;
    Real outletRCR = 0.0;
    Real csv = 0.0;
    Real history = 0.0;
    Real output = 0.0;
  };

  struct MaxReal
  {
    Real operator()(Real lhs, Real rhs) const
    {
      return std::max(lhs, rhs);
    }
  };

  inline Real secondsSince(CoronaryClock::time_point start)
  {
    return std::chrono::duration<Real>(CoronaryClock::now() - start).count();
  }
}

#endif
