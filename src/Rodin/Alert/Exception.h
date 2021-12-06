/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_EXCEPTION_H
#define RODIN_ALERT_EXCEPTION_H

#include "Alert.h"

namespace Rodin::Alert
{
   class Exception : public Alert
   {
      public:
         Exception() = default;

         Exception(const Exception& other) = default;

         Exception(const std::string& what);

         virtual void raise() override;

         virtual Exception& operator<<(const std::string& ss) override;
   };
}

#endif
