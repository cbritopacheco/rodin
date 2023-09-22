/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_WARNING_H
#define RODIN_ALERT_WARNING_H

#define RODIN_ALERT_WARNING_PREFIX "Warning"

#include "Message.h"

namespace Rodin::Alert
{
  class WarningPrefix : public MessagePrefix<YellowT>
  {
    public:
      using Parent = MessagePrefix<YellowT>;

      WarningPrefix()
        : Parent(RODIN_ALERT_WARNING_PREFIX)
      {}
  };

  /**
   * @brief Derived Alert class representing a warning.
   *
   * Represents an alert which does not terminate the program when raised, yet
   * still has visible effects when raised. For example, showing a formatted
   * message on the screen.
   */
  class Warning : public Message<WarningPrefix>
  {
    public:
      using Parent = Message<WarningPrefix>;

      /**
       * @brief Constructs a Warning alert with an empty message.
       */
      Warning();

      Warning(std::ostream& os);

      /**
       * @brief Copy constructor.
       */
      Warning(const Warning& other);

      /**
       * @brief Move constructor.
       */
      Warning(Warning&& other);
  };
}

#endif
