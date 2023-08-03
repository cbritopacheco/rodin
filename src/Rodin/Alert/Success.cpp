/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <iostream>
#include <rang.hpp>

#include "Success.h"

namespace Rodin::Alert
{
  Success::Success()
    : Alert(RODIN_ALERT_SUCCESS_PREFIX_LENGTH)
  {}

  void Success::raise() const noexcept
  {
    std::cout << rang::fg::green
           << RODIN_ALERT_SUCCESS_PREFIX
           << rang::fg::reset
           << what()
           << std::endl;
  }
}

