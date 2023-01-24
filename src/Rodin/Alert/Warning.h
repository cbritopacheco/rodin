/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_WARNING_H
#define RODIN_ALERT_WARNING_H

#include "Alert.h"

namespace Rodin::Alert
{
  /**
   * @brief Derived Alert class representing a warning.
   *
   * Represents an alert which does not terminate the program when raised, yet
   * still has visible effects when raised. For example, showing a formatted
   * message on the screen.
   */
  class Warning : public Alert
  {
    public:
      /**
       * @brief Constructs a Warning with an empty message.
       */
      Warning() = default;

      /**
       * @brief Constructs a Warning with the given message.
       * @param[in] what Description or reason for the Warning being raised.
       */
      Warning(const std::string& what);

      /**
       * @brief Copies the Warning message.
       */
      Warning(const Warning&) = default;

      /**
       * @brief Raises the warning to the user.
       *
       * Default behaviour is to output a formatted warning message to
       * stderr.
       */
      virtual void raise() const noexcept override;
  };
}

#endif
