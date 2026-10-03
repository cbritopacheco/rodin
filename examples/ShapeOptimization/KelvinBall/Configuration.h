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
      /// Only controls construction of the initial uniform chamber.
      Real requestedGridSpacing = 0;
      Real hminFactor = 0.1;
      Real hmaxFactor = 10;
      Real hmin = 0;
      Real hmax = 0;
      Real outerRadius = 2;
      Real nitschePenalty = DefaultNitschePenalty;
      Real stabilizationFactor = DefaultStabilizationFactor;
      /// Settings of the initial WNGIR background preparation.
      Real backgroundHausdorff = 0.05;
      Real backgroundGradation = 2;
      /// Optional MMG adaptation after each MMG cut or WNGIR fit, and during
      /// initial background preparation. Its size map ranges from hmin to hmax.
      bool adapt = false;
      Real adaptGradation = 1.3;
      /// Crossing fraction below which the MMG cut snaps the near endpoint of
      /// a crossed edge onto the level set; zero disables the snapping.
      Real mmgSnap = 0;
      /// Number of times a failed MMG reconstruction is retried with the MMG
      /// sizes computed from half the previous scale.
      size_t mmgRetries = 2;
      bool pointsSpecified = false;
      bool gridSpacingSpecified = false;

      bool parse(std::string_view option);

      void finalize();

      Real getGridSpacing() const;
  };
}

#endif
