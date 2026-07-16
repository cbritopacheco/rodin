/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Common.h"

namespace Rodin::MMG
{
  const char* getISCDMshdistExecutable()
  {
    return iscdMshdistExecutable;
  }

  const char* getISCDAdvectExecutable()
  {
    return iscdAdvectionExecutable;
  }

  int getMMGVerbosityLevel()
  {
    return verbosityLevel;
  }
}
