/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_WARNING_H
#define RODIN_ALERT_WARNING_H

/**
 * @file
 * @brief Warning messages with yellow-colored output.
 */

#define RODIN_ALERT_WARNING_PREFIX "Warning"

#include "Message.h"

namespace Rodin::Alert
{
  /**
   * @brief Prefix class for warning messages.
   * @ingroup AlertModule
   *
   * Provides yellow-colored prefixing for warning messages using "Warning" as
   * the default prefix text.
   */
  class WarningPrefix : public MessagePrefix<YellowT>
  {
    public:
      /// @brief Parent class type alias
      using Parent = MessagePrefix<YellowT>;

      /**
       * @brief Default constructor.
       *
       * Initializes the warning prefix with the default "Warning" text.
       */
      WarningPrefix()
        : Parent(RODIN_ALERT_WARNING_PREFIX)
      {}
  };

  /**
   * @brief Warning message class with formatted output.
   * @ingroup AlertModule
   *
   * A specialized message class for displaying warnings with yellow-colored
   * formatting. Warnings do not terminate the program and are typically used
   * to alert users about potential issues or non-critical problems.
   *
   * The Warning class provides stream-like insertion operators for building
   * detailed warning messages with formatting capabilities.
   */
  class Warning : public Message<WarningPrefix>
  {
    public:
      /// @brief Parent class type alias
      using Parent = Message<WarningPrefix>;

      /**
       * @brief Constructs a Warning alert with an empty message.
       *
       * Creates a warning message that outputs to std::cout by default.
       */
      Warning();

      /**
       * @brief Constructs a Warning alert with a specific output stream.
       * @param os The output stream to write the warning message to
       */
      Warning(std::ostream& os);

      /**
       * @brief Copy constructor.
       * @param other The Warning object to copy from
       */
      Warning(const Warning& other);

      /**
       * @brief Move constructor.
       * @param other The Warning object to move from
       */
      Warning(Warning&& other);
  };
}

#endif
