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
   class Warning : public Alert
   {
      public:
         Warning() = default;

         Warning(const std::string& what);

         Warning(const Warning& other) = default;

         virtual void raise() override;

         void operator<<(const RaiseT&)
         {
            this->raise();
         }

         template <class T>
         Warning& operator<<(T&& v)
         {
            return static_cast<Warning&>(Alert::operator<<(std::forward<T>(v)));
         }
   };
}

#endif