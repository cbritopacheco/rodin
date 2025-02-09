/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cstdint>
#include <cstdlib>
#include <iostream>

#include "Success.h"

namespace Rodin::Alert
{
  Success::Success()
    : Success(std::cout)
  {}

  Success::Success(std::ostream& os)
    : Parent(os, SuccessPrefix())
  {}

  Success::Success(const Success& other)
    : Parent(other)
  {}

  Success::Success(Success&& other)
    : Parent(std::move(other))
  {}
}

