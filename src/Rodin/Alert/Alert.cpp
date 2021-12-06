/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <iostream>

#include "Alert.h"

namespace Rodin::Alert
{
   Alert::Alert(const std::string& what)
   {
      m_what << what;
   }

   std::string Alert::what() const
   {
      return m_what.str();
   }

   Alert& Alert::operator<<(const std::string& ss)
   {
      m_what << ss;
      return *this;
   }
}
