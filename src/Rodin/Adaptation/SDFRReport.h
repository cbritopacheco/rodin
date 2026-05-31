/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_SDFRREPORT_H
#define RODIN_ADAPTATION_SDFRREPORT_H

#include <cstddef>

#include "Rodin/Types.h"

namespace Rodin::Adaptation
{
  struct SDFRReport
  {
    bool converged = false;
    bool lineSearchFailed = false;

    std::size_t iterations = 0;
    std::size_t totalBacktracks = 0;
    std::size_t inadmissibleCount = 0;

    Real initialResidual = 0;
    Real finalResidual = 0;
    Real finalStepNorm = 0;
    Real minJRatio = 0;
    Real jLineSearchRatio = 0;
    Real lastAcceptedAlpha = 0;
  };
}

#endif
