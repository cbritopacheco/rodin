/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_NEWLINE_H
#define RODIN_ALERT_NEWLINE_H

#include <ostream>

namespace Rodin::Alert
{
  /**
   * @brief Empty tag type for inserting newline characters.
   * @ingroup AlertModule
   *
   * Tag type used to insert newline characters into Alert messages.
   * When streamed to an output stream or Message object, it produces
   * a newline and triggers proper indentation for subsequent content.
   */
  struct NewLineT
  {
    explicit constexpr NewLineT() = default;
  };

  /**
   * @brief Stream insertion operator for NewLineT.
   * @ingroup AlertModule
   * @param os The output stream to write to.
   * @return Reference to the output stream.
   *
   * Inserts a newline character into the output stream.
   */
  inline
  std::ostream& operator<<(std::ostream& os, const NewLineT&)
  {
    os << '\n';
    return os;
  }

  /**
   * @brief Instance of NewLineT tag type.
   * @ingroup AlertModule
   *
   * Constant instance of the NewLineT tag type for convenient usage.
   * Use this to insert newlines in Alert messages with proper indentation.
   */
  static constexpr NewLineT NewLine;
}

#endif

