/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Common.h"

namespace Rodin::External::MMG
{
  const char* getISCDMshdistExecutable()
  {
   return ISCD_MSHDIST_EXECUTABLE;
  }

  const char* getISCDAdvectExecutable()
  {
   return ISCD_ADVECTION_EXECUTABLE;
  }

  int getMMGVerbosityLevel()
  {
   return VERBOSITY_LEVEL;
  }
}
