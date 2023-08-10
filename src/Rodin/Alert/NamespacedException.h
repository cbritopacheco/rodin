/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_NAMESPACEDEXCEPTION_H
#define RODIN_ALERT_NAMESPACEDEXCEPTION_H

#include "Exception.h"
#include "Identifier.h"

namespace Rodin::Alert
{
  class NamespacedException : public Exception
  {
    public:
      NamespacedException(const std::string& namesp)
      {
        *this << Identifier::Namespace(namesp) << ". ";
      }
  };
}

#endif
