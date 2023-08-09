/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <iostream>
#include <rang.hpp>

#include "Rodin/Configure.h"

#include "Warning.h"

namespace Rodin::Alert
{
  Warning::Warning(const std::string& what)
    : Alert(what, RODIN_ALERT_WARNING_PREFIX_LENGTH)
  {}

  void Warning::raise() const noexcept
  {
#ifdef RODIN_SILENCE_WARNINGS
#else
    std::cerr << rang::fg::yellow
              << RODIN_ALERT_WARNING_PREFIX
              << rang::fg::reset
              << what()
              << std::endl;
#endif
  }
}

