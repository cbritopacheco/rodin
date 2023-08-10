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
  class InfoPrefix : public MessagePrefix<BlueT>
  {
    public:
      using Parent = MessagePrefix<BlueT>;

      InfoPrefix()
        : Parent(RODIN_ALERT_INFO_PREFIX)
      {}
  };

  class Info : public Message<InfoPrefix>
  {
    public:
      using Parent = Message<InfoPrefix>;

      /**
       * @brief Constructs an Info alert with an empty message.
       */
      Info();

      Info(std::ostream& os);

      /**
       * @brief Copy constructor.
       */
      Info(const Info& other);

      /**
       * @brief Move constructor.
       */
      Info(Info&& other);
  };
}

#endif
