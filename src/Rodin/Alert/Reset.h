/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_RESET_H
#define RODIN_ALERT_RESET_H

#include <ostream>
#include <termcolor/termcolor.hpp>

namespace Rodin::Alert
{
  struct ResetT
  {
    constexpr
    ResetT() = default;

    constexpr
    ResetT(const ResetT&) = default;

    constexpr
    ResetT(ResetT&&) = default;

    constexpr
    ResetT& operator=(const ResetT&) = default;

    constexpr
    ResetT& operator=(ResetT&&) = default;
  };

  static constexpr ResetT Reset;

  inline
  std::ostream& operator<<(std::ostream& os, const ResetT&)
  {
    os << termcolor::reset;
    return os;
  }
}

#endif
