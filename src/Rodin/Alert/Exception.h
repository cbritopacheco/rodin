/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_EXCEPTION_H
#define RODIN_ALERT_EXCEPTION_H

#define RODIN_ALERT_EXCEPTION_PREFIX "Error"

#include <exception>

#include "Message.h"

namespace Rodin::Alert
{
  class ExceptionPrefix : public MessagePrefix<RedT>
  {
    public:
      using Parent = MessagePrefix<RedT>;

      ExceptionPrefix()
        : Parent(RODIN_ALERT_EXCEPTION_PREFIX)
      {}
  };

  /**
   * @brief Derived Alert class representing an exception.
   *
   * Represents an alert which, when raised, will terminate the program after
   * outputting a message with the reason why.
   */
  class Exception : public std::exception, public Message<ExceptionPrefix>
  {
    public:
      using Parent = Message<ExceptionPrefix>;

      /**
       * @brief Constructs an Exception with an empty message.
       */
      Exception();

      Exception(std::ostream& os);

      /**
       * @brief Copy constructor.
       */
      Exception(const Exception& other);

      /**
       * @brief Move constructor.
       */
      Exception(Exception&& other);

      /**
       * @brief Raises the exception to the user.
       *
       * Default behaviour is to output a formatted error message and call
       * std::abort.
       */
      virtual void raise() const override;

      const char* what() const noexcept override;
  };
}

#endif
