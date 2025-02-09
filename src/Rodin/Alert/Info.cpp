/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cstdint>
#include <cstdlib>
#include <iostream>

#include "Info.h"

namespace Rodin::Alert
{
  Info::Info()
    : Info(std::cout)
  {}

  Info::Info(std::ostream& os)
    : Parent(os, InfoPrefix())
  {}

  Info::Info(const Info& other)
    : Parent(other)
  {}

  Info::Info(Info&& other)
    : Parent(std::move(other))
  {}
}
