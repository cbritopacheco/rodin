/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <iostream>
#include <rang.hpp>

#include "Rodin/Configure.h"

#include "Exception.h"

namespace Rodin::Alert
{
  Exception::Exception(const std::string& what)
    : Alert(what, 0)
  {}

  void Exception::raise() const
  {
#ifdef RODIN_SILENCE_EXCEPTIONS
#else
    std::cerr << rang::fg::red
              << RODIN_ALERT_WARNING_PREFIX
              << rang::fg::reset
              << what()
              << std::endl;
#endif
    throw *this;
  }
}

