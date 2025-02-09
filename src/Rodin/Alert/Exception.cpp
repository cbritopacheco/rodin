/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cstdint>
#include <cstdlib>
#include <iostream>

#include "Rodin/Configure.h"

#include "Exception.h"

namespace Rodin::Alert
{
  Exception::Exception()
    : Exception(std::cerr)
  {}

  Exception::Exception(std::ostream& os)
    : Parent(os, ExceptionPrefix())
  {}

  Exception::Exception(const Exception& other)
    : std::exception(other),
      Parent(other)
  {}

  Exception::Exception(Exception&& other)
    : std::exception(std::move(other)),
      Parent(std::move(other))
  {}

  const char* Exception::what() const noexcept
  {
    return Parent::what();
  }

  void Exception::raise() const
  {
    Parent::raise();
    throw *this;
  }
}

