/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ATRIUM_RUNPARAMETERS_H
#define RODIN_ATRIUM_RUNPARAMETERS_H

#include <string>
#include <Rodin/Types.h>

namespace Rodin::Examples::Heart::Atrium
{
  /**
   * @brief Common time and output configuration passed unchanged to both components.
   *
   * The time step is in seconds. The step count measures physical advances, excluding the
   * initial state. An output interval of zero retains only the initial and final
   * states (fluid cycle boundaries are also exported). CSV diagnostics are
   * recorded on every accepted step independently of the export interval.
   */
  struct RunParameters
  {
    Real timeStep = 1e-3;
    size_t steps = 24;
    size_t outputEvery = 1;
    std::string outputDirectory = "output/atrium";
  };
}

#endif
