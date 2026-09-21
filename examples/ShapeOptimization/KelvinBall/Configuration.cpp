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
      requestedH = std::stod(std::string(option.substr(4)));
      hSpecified = true;
    }
    else if (option.rfind("--outer-radius=", 0) == 0)
      outerRadius = std::stod(std::string(option.substr(15)));
    else if (option.rfind("--penalty=", 0) == 0)
      nitschePenalty = std::stod(std::string(option.substr(10)));
    else if (option.rfind("--stabilization=", 0) == 0)
      stabilizationFactor = std::stod(std::string(option.substr(16)));
    else if (option.rfind("--background-hmin=", 0) == 0)
      backgroundHMin = std::stod(std::string(option.substr(18)));
    else if (option.rfind("--background-hmax=", 0) == 0)
      backgroundHMax = std::stod(std::string(option.substr(18)));
    else if (option.rfind("--background-hausdorff=", 0) == 0)
      backgroundHausdorff = std::stod(std::string(option.substr(23)));
    else if (option.rfind("--background-gradation=", 0) == 0)
      backgroundGradation = std::stod(std::string(option.substr(23)));
    else
      return false;
    return true;
  }

  void Configuration::finalize()
  {
    if (pointsSpecified && hSpecified)
      throw std::runtime_error("Use either --n or --h, not both.");
    if (!(outerRadius > 1))
      throw std::runtime_error(
        "The chamber outer radius must exceed the unit sphere radius.");
    if (hSpecified)
    {
      if (!(requestedH > 0))
        throw std::runtime_error("The requested mesh size must be positive.");
      points = static_cast<size_t>(std::ceil(outerRadius / requestedH)) + 1;
    }
    if (points < 5)
      throw std::runtime_error("The chamber uniform grid requires --n >= 5.");
    if (!(nitschePenalty > 0))
      throw std::runtime_error("The Nitsche penalty must be positive.");
    if (stabilizationFactor < 0)
      throw std::runtime_error("The pressure stabilization must be nonnegative.");
    if (!(backgroundHMin > 0) || !(backgroundHMax >= backgroundHMin))
      throw std::runtime_error(
        "The background sizes must satisfy 0 < --background-hmin <= --background-hmax.");
    if (!(backgroundHausdorff > 0))
      throw std::runtime_error("The background Hausdorff tolerance must be positive.");
    if (!(backgroundGradation > 1))
      throw std::runtime_error("The background gradation must exceed one.");
  }

  Real Configuration::getH() const
  {
    return outerRadius / static_cast<Real>(points - 1);
  }
}
