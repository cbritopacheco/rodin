/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_EXCEPTION_H
#define RODIN_ALERT_EXCEPTION_H

/**
 * @defgroup AlertModule Alert
 * @brief Alert and exception handling system for Rodin.
 *
 * The Alert module provides a comprehensive system for handling exceptions,
 * warnings, informational messages, and success notifications with colored
 * terminal output and formatted messaging capabilities.
 */

#define RODIN_ALERT_EXCEPTION_PREFIX "Error"

#include <exception>

#include "Message.h"

namespace Rodin::Alert
{
  /**
   * @brief Prefix class for exception messages.
   * @ingroup AlertModule
   *
   * Provides red-colored prefixing for exception messages using "Error" as
   * the default prefix text.
   */
  class ExceptionPrefix : public MessagePrefix<RedT>
  {
    public:
      /// @brief Parent class type alias.
      using Parent = MessagePrefix<RedT>;

      /**
       * @brief Default constructor.
       *
       * Initializes the exception prefix with the default "Error" text.
       */
      ExceptionPrefix()
        : Parent(RODIN_ALERT_EXCEPTION_PREFIX)
      {}
  };

  /**
   * @brief Exception class with formatted output capabilities.
   * @ingroup AlertModule
   *
   * A specialized exception class that extends both std::exception and the 
   * Alert messaging system. When raised, it outputs a formatted error message
   * with red coloring and terminates the program.
   *
   * The Exception class provides stream-like insertion operators for building
   * detailed error messages with formatting capabilities.
   */
  class Exception : public std::exception, public Message<ExceptionPrefix>
  {
    public:
      /// @brief Parent class type alias.
      using Parent = Message<ExceptionPrefix>;

      /**
       * @brief Constructs an Exception with an empty message.
       *
       * Creates an exception that outputs to std::cerr by default.
       */
      Exception();

      /**
       * @brief Constructs an Exception with a specific output stream.
       * @param os The output stream to write the error message to.
       */
      Exception(std::ostream& os);

      /**
       * @brief Copy constructor.
       * @param other The Exception object to copy from.
       */
      Exception(const Exception& other);

      /**
       * @brief Move constructor.
       * @param other The Exception object to move from.
       */
      Exception(Exception&& other);

      /**
       * @brief Raises the exception to the user.
       *
       * Outputs a formatted error message to the configured stream and
       * terminates the program by calling std::abort(). The message includes
       * the exception prefix and any accumulated message content.
       */
      virtual void raise() const override;

      /**
       * @brief Returns the exception message.
       * @return A C-style string containing the exception message.
       *
       * This method provides compatibility with std::exception::what().
       */
      const char* what() const noexcept override;
  };
}

#endif
