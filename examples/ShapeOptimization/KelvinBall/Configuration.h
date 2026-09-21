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
      /// Settings of the single MMG optimization of the WNGIR background, the
      /// sizes in multiples of h.
      Real backgroundHMin = 0.1;
      Real backgroundHMax = 1;
      Real backgroundHausdorff = 0.05;
      Real backgroundGradation = 2;
      /// Optional MMG adaptation that replaces the optimization pass after each
      /// level-set cut. The size grows linearly from the interface size on
      /// Gamma to the far size at the given width; all in multiples of h.
      bool adapt = false;
      Real adaptInterfaceSize = 1;
      Real adaptFarSize = 1;
      Real adaptWidth = 3;
      Real adaptGradation = 1.3;
      bool pointsSpecified = false;
      bool hSpecified = false;

      bool parse(std::string_view option);

      void finalize();

      Real getH() const;
  };
}

#endif
