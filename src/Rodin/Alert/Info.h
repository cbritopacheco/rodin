/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_INFO_H
#define RODIN_ALERT_INFO_H

#include "Alert.h"

namespace Rodin::Alert
{
   class Info : public Alert
   {
      public:
         /**
          * @brief Constructs a Info with an empty message.
          */
         Info() = default;

         /**
          * @brief Constructs a Info with the given message.
          * @param[in] what Description or reason for the Info being raised.
          */
         Info(const std::string& what)
            : Alert(what)
         {}

         /**
          * @brief Copies the Info message.
          */
         Info(const Info&) = default;

         /**
          * @brief Raises the Info to the user.
          *
          * Default behaviour is to output a formatted Info message to
          * stderr.
          */
         virtual void raise() const noexcept override;
   };
}

#endif
