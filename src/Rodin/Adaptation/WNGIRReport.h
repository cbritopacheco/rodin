/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_WNGIRREPORT_H
#define RODIN_ADAPTATION_WNGIRREPORT_H

#include "WNGIRCommon.h"

namespace Rodin::Adaptation
{
  struct WNGIRReport
  {
    std::size_t iterations = 0;
    Real sigma = 0;
    Real lastAlpha = 0;
    Real acceptedStep = 0;
    Real minJ = 1;
    Real maxQRel = 1;
    Real activeRMS = 0;
    Real activeSup = 0;
    Real activeFraction = 0;
    Real energy = 0;
    const char* exitReason = "iter-budget";
    // Wall-clock breakdown (seconds, accumulated over iterations).
    Real tAssembly = 0;   ///< WNGIR variational problem assembly.
    Real tFactor = 0;     ///< CG setup/preconditioner.
    Real tSolve = 0;      ///< CG iterations.
    Real tLineSearch = 0; ///< true-geometry admissibility + energy LS.
    std::size_t linearIterations = 0;
    Real linearError = 0;
    std::size_t andersonTried = 0;
    std::size_t andersonAccepted = 0;
    Real lastAndersonTheta = 0;
  };
}

#endif
