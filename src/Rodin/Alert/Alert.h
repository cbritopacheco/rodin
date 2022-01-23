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
   struct RaiseT
   {};

   static constexpr RaiseT Raise;

   class Alert
   {
      public:
         Alert() = default;

         Alert(const std::string& what);

         Alert(const Alert& other);

         virtual ~Alert() = default;

         std::string what() const;


         virtual void raise() = 0;

         template <class T>
         Alert& operator<<(T&& v)
         {
            m_what << std::forward<T>(v);
            return *this;
         }

      private:
         std::stringstream m_what;
   };
}

#endif
