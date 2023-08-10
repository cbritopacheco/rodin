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
  class SuccessPrefix : public MessagePrefix<GreenT>
  {
    public:
      using Parent = MessagePrefix<GreenT>;

      SuccessPrefix()
        : Parent(RODIN_ALERT_SUCCESS_PREFIX)
      {}
  };

  class Success : public Message<SuccessPrefix>
  {
    public:
      using Parent = Message<SuccessPrefix>;

      /**
       * @brief Constructs a Success alert with an empty message.
       */
      Success();

      Success(std::ostream& os);

      /**
       * @brief Copy constructor.
       */
      Success(const Success& other);

      /**
       * @brief Move constructor.
       */
      Success(Success&& other);
  };
}

#endif

