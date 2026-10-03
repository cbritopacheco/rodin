/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Configuration.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace KelvinBall
{
  bool Configuration::parse(std::string_view option)
  {
    if (option.rfind("--n=", 0) == 0)
    {
      points = std::stoul(std::string(option.substr(4)));
      pointsSpecified = true;
    }
    else if (option.rfind("--h=", 0) == 0)
    {
      requestedGridSpacing = std::stod(std::string(option.substr(4)));
      gridSpacingSpecified = true;
    }
    else if (option.rfind("--hmin-factor=", 0) == 0)
      hminFactor = std::stod(std::string(option.substr(14)));
    else if (option.rfind("--hmax-factor=", 0) == 0)
      hmaxFactor = std::stod(std::string(option.substr(14)));
    else if (option.rfind("--outer-radius=", 0) == 0)
      outerRadius = std::stod(std::string(option.substr(15)));
    else if (option.rfind("--penalty=", 0) == 0)
      nitschePenalty = std::stod(std::string(option.substr(10)));
    else if (option.rfind("--stabilization=", 0) == 0)
      stabilizationFactor = std::stod(std::string(option.substr(16)));
    else if (option.rfind("--background-hausdorff=", 0) == 0)
      backgroundHausdorff = std::stod(std::string(option.substr(23)));
    else if (option.rfind("--background-gradation=", 0) == 0)
      backgroundGradation = std::stod(std::string(option.substr(23)));
    else if (option == "--mmg-adapt")
      adapt = true;
    else if (option.rfind("--mmg-adapt-gradation=", 0) == 0)
      adaptGradation = std::stod(std::string(option.substr(22)));
    else if (option.rfind("--mmg-snap=", 0) == 0)
      mmgSnap = std::stod(std::string(option.substr(11)));
    else if (option.rfind("--mmg-retries=", 0) == 0)
      mmgRetries = std::stoul(std::string(option.substr(14)));
    else
      return false;
    return true;
  }

  void Configuration::finalize()
  {
    if (pointsSpecified && gridSpacingSpecified)
      throw std::runtime_error("Use either --n or --h, not both.");
    if (!(outerRadius > 1))
      throw std::runtime_error(
        "The chamber outer radius must exceed the unit sphere radius.");
    if (gridSpacingSpecified)
    {
      if (!std::isfinite(requestedGridSpacing) || !(requestedGridSpacing > 0))
        throw std::runtime_error("The requested initial grid spacing must be positive.");
      points = static_cast<size_t>(std::ceil(outerRadius / requestedGridSpacing)) + 1;
    }
    if (points < 5)
      throw std::runtime_error("The chamber uniform grid requires at least five points per edge.");
    if (!std::isfinite(hminFactor) || !std::isfinite(hmaxFactor) ||
        !(hminFactor > 0) || !(hmaxFactor >= hminFactor))
      throw std::runtime_error(
        "The size factors must satisfy 0 < --hmin-factor <= --hmax-factor.");
    hmin = hminFactor * getGridSpacing();
    hmax = hmaxFactor * getGridSpacing();
    if (!(nitschePenalty > 0))
      throw std::runtime_error("The Nitsche penalty must be positive.");
    if (stabilizationFactor < 0)
      throw std::runtime_error("The pressure stabilization must be nonnegative.");
    if (!(backgroundHausdorff > 0))
      throw std::runtime_error("The background Hausdorff tolerance must be positive.");
    if (!(backgroundGradation > 1))
      throw std::runtime_error("The background gradation must exceed one.");
    if (!(adaptGradation > 1))
      throw std::runtime_error("The adaptation gradation must exceed one.");
    if (!(mmgSnap >= 0) || !(mmgSnap < 0.5))
      throw std::runtime_error(
        "The snapping fraction must satisfy 0 <= --mmg-snap < 0.5.");
  }

  Real Configuration::getGridSpacing() const
  {
    return outerRadius / static_cast<Real>(points - 1);
  }
}
