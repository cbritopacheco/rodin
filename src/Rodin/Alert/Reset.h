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
  /**
   * @brief Empty tag type for resetting terminal formatting.
   *
   * Tag type used to reset all terminal text formatting and colors
   * back to default. This includes clearing any colors, bold, italic,
   * underline, and other attributes previously applied.
   */
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

  /**
   * @brief Instance of ResetT tag type.
   *
   * Constant instance of the ResetT tag type for convenient usage.
   * Use this to reset terminal formatting to default.
   */
  static constexpr ResetT Reset;

  /**
   * @brief Stream insertion operator for ResetT.
   * @param os The output stream to write to.
   * @return Reference to the output stream.
   *
   * Resets all terminal formatting and colors to default using the
   * termcolor library.
   */
  inline
  std::ostream& operator<<(std::ostream& os, const ResetT&)
  {
    os << termcolor::reset;
    return os;
  }
}

#endif
