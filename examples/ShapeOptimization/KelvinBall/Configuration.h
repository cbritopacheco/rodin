/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef KELVIN_BALL_CONFIGURATION_H
#define KELVIN_BALL_CONFIGURATION_H

#include <string_view>

#include <Rodin/Types.h>

#include "Common.h"

namespace KelvinBall
{
  /** Common chamber, mesh, and Stokes parameters. */
  struct Configuration
  {
      size_t points = 13;
      Real requestedH = 0;
      Real outerRadius = 2;
      Real nitschePenalty = DefaultNitschePenalty;
      Real stabilizationFactor = DefaultStabilizationFactor;
      bool pointsSpecified = false;
      bool hSpecified = false;

      bool parse(std::string_view option);

      void finalize();

      Real getH() const;
  };
}

#endif
