/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <iostream>
#include <rang.hpp>

#include "Info.h"

namespace Rodin::Alert
{
  Info::Info()
    : Alert(6)
  {}

  void Info::raise() const noexcept
  {
    std::cout << rang::fg::blue
           << "Info: "
           << rang::fg::reset
           << what()
           << std::endl;
  }
}
