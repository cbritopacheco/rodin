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
  /**
   * @brief Empty tag type for enabling terminal color output.
   *
   * Tag type used to enable colored output in terminal streams.
   * When streamed to an output stream, it activates the termcolor
   * library's colorization for subsequent output.
   */
  struct StylizeT {};

  /**
   * @brief Instance of StylizeT tag type.
   *
   * Constant instance of the StylizeT tag type for convenient usage.
   * Use this to enable colored terminal output.
   */
  static constexpr StylizeT Stylize;

  /**
   * @brief Stream insertion operator for StylizeT.
   * @param os The output stream to write to.
   * @return Reference to the output stream.
   *
   * Enables terminal colorization for the output stream using the
   * termcolor library.
   */
  inline
  std::ostream& operator<<(std::ostream& os, const StylizeT&)
  {
    os << termcolor::colorize;
    return os;
  }
}

#endif
