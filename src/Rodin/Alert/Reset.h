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

    /**
     * @brief Copy constructor.
     * @param other Object to copy from.
     */
    constexpr ResetT(const ResetT& other) = default;

    /**
     * @brief Move constructor.
     * @param other Object to move from.
     */
    constexpr ResetT(ResetT&& other) = default;

    /**
     * @brief Copy assignment operator.
     * @returns Reference to this object after the operation.
     * @param other Object to copy from.
     */
    constexpr ResetT& operator=(const ResetT& other) = default;

    /**
     * @brief Move assignment operator.
     * @returns Reference to this object after the operation.
     * @param other Object to move from.
     */
    constexpr ResetT& operator=(ResetT&& other) = default;
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
   * @param tag Formatting or action tag selected through its type.
   */
  inline std::ostream& operator<<(std::ostream& os, [[maybe_unused]] const ResetT& tag)
  {
    os << termcolor::reset;
    return os;
  }
}

#endif
