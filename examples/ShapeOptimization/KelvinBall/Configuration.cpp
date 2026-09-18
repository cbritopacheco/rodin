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
  }

  Real Configuration::getH() const
  {
    return outerRadius / static_cast<Real>(points - 1);
  }
}
