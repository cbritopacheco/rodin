/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_STYLIZE_H
#define RODIN_ALERT_STYLIZE_H

#include <ostream>
#include <termcolor/termcolor.hpp>

namespace Rodin::Alert
{
  struct StylizeT {};

  static constexpr StylizeT Stylize;

  inline
  std::ostream& operator<<(std::ostream& os, const StylizeT&)
  {
    os << termcolor::colorize;
    return os;
  }
}

#endif
