/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_RODINEXTERNAL_MMG_CONFIGURE_H
#define RODIN_RODINEXTERNAL_MMG_CONFIGURE_H

namespace Rodin::External::MMG
{
  /**
   * @internal
   * @brief Path to the ISCD Mshdist executable.
   */
  inline constexpr const char* ISCD_MSHDIST_EXECUTABLE =
    "/Users/carlos/Projects/rodin/build2/third-party/ISCD/Mshdist/mshdist";

  inline constexpr const char* MSHDIST_EXECUTABLE =
    "/Users/carlos/Projects/rodin/build2/third-party/ISCD/Mshdist/mshdist";

  /**
   * @internal
   * @brief Path to the ISCD Mshdist executable.
   */
  inline constexpr const char* ISCD_ADVECTION_EXECUTABLE =
    "/Users/carlos/Projects/rodin/build2/third-party/ISCD/Advection/Advection";

  inline constexpr const char* ADVECTION_EXECUTABLE =
    "/Users/carlos/Projects/rodin/build2/third-party/ISCD/Advection/Advection";

  /**
   * @brief Verbosity level for MMG console output.
   *
   * Ranges from -1 to INT_MAX. A higher number means more verbose output.
   */
  inline constexpr const int VERBOSITY_LEVEL = 1;
}

#endif
