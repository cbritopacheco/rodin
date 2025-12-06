/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_SUCCESS_H
#define RODIN_ALERT_SUCCESS_H

#define RODIN_ALERT_SUCCESS_PREFIX "Success"

#include "Message.h"

namespace Rodin::Alert
{
  /**
   * @brief Prefix class for success messages.
   *
   * Provides green-colored prefixing for success messages using "Success" as
   * the default prefix text.
   */
  class SuccessPrefix : public MessagePrefix<GreenT>
  {
    public:
      /// @brief Parent class type alias.
      using Parent = MessagePrefix<GreenT>;

      /**
       * @brief Default constructor.
       *
       * Initializes the success prefix with the default "Success" text.
       */
      SuccessPrefix()
        : Parent(RODIN_ALERT_SUCCESS_PREFIX)
      {}
  };

  /**
   * @brief Success message class with formatted output.
   *
   * A specialized message class for displaying success notifications with
   * green-colored formatting. Success messages indicate successful completion
   * of operations and are typically used for user feedback and logging.
   *
   * The Success class provides stream-like insertion operators for building
   * detailed success messages with formatting capabilities.
   */
  class Success : public Message<SuccessPrefix>
  {
    public:
      /// @brief Parent class type alias.
      using Parent = Message<SuccessPrefix>;

      /**
       * @brief Constructs a Success alert with an empty message.
       *
       * Creates a success message that outputs to std::cout by default.
       */
      Success();

      /**
       * @brief Constructs a Success alert with a specific output stream.
       * @param os The output stream to write the success message to.
       */
      Success(std::ostream& os);

      /**
       * @brief Copy constructor.
       * @param other The Success object to copy from.
       */
      Success(const Success& other);

      /**
       * @brief Move constructor.
       * @param other The Success object to move from.
       */
      Success(Success&& other);
  };
}

#endif

