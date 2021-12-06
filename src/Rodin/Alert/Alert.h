/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_ALERT_H
#define RODIN_ALERT_ALERT_H

#include <string>
#include <sstream>

namespace Rodin::Alert
{
   class Alert
   {
      public:
         Alert() = default;

         Alert(const std::string& what);

         virtual void raise() = 0;

         std::string what() const;

         Alert& operator<<(const std::string& ss);

      private:
         std::stringstream m_what;
   };
}

#endif
