/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_SUCCESS_H
#define RODIN_ALERT_SUCCESS_H

#define RODIN_ALERT_SUCCESS_PREFIX "Success: "
#define RODIN_ALERT_SUCCESS_PREFIX_LENGTH (sizeof(RODIN_ALERT_SUCCESS_PREFIX) - 1)

#include "Alert.h"

namespace Rodin::Alert
{
  class Success : public Alert
  {
    public:
      /**
       * @brief Constructs a Success with an empty message.
       */
      Success();

      /**
       * @brief Constructs a Success with the given message.
       * @param[in] what Description or reason for the Success being raised.
       */
      Success(const std::string& what)
        : Alert(what, RODIN_ALERT_SUCCESS_PREFIX_LENGTH)
      {}

      /**
       * @brief Copies the Success message.
       */
      Success(const Success&) = default;

      /**
       * @brief Raises the Success to the user.
       *
       * Default behaviour is to output a formatted Success message to
       * stderr.
       */
      virtual void raise() const noexcept override;
  };
}

#endif

