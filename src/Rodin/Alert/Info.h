/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_INFO_H
#define RODIN_ALERT_INFO_H

#define RODIN_ALERT_INFO_PREFIX "Info"

#include "Message.h"

namespace Rodin::Alert
{
  /**
   * @brief Prefix class for informational messages.
   * @ingroup AlertModule
   *
   * Provides blue-colored prefixing for informational messages using "Info" as
   * the default prefix text.
   */
  class InfoPrefix : public MessagePrefix<BlueT>
  {
    public:
      /// @brief Parent class type alias.
      using Parent = MessagePrefix<BlueT>;

      /**
       * @brief Default constructor.
       *
       * Initializes the info prefix with the default "Info" text.
       */
      InfoPrefix()
        : Parent(RODIN_ALERT_INFO_PREFIX)
      {}
  };

  /**
   * @brief Informational message class with formatted output.
   * @ingroup AlertModule
   *
   * A specialized message class for displaying informational content with
   * blue-colored formatting. Unlike exceptions, info messages do not terminate
   * the program and are typically used for logging and user notifications.
   *
   * The Info class provides stream-like insertion operators for building
   * detailed informational messages with formatting capabilities.
   */
  class Info : public Message<InfoPrefix>
  {
    public:
      /// @brief Parent class type alias.
      using Parent = Message<InfoPrefix>;

      /**
       * @brief Constructs an Info alert with an empty message.
       *
       * Creates an informational message that outputs to std::cout by default.
       */
      Info();

      /**
       * @brief Constructs an Info alert with a specific output stream.
       * @param os The output stream to write the informational message to.
       */
      Info(std::ostream& os);

      /**
       * @brief Copy constructor.
       * @param other The Info object to copy from.
       */
      Info(const Info& other);

      /**
       * @brief Move constructor.
       * @param other The Info object to move from.
       */
      Info(Info&& other);
  };
}

#endif
