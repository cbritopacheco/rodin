/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_CLASSEXCEPTION_H
#define RODIN_ALERT_CLASSEXCEPTION_H

#include <boost/type_index.hpp>

#include "Exception.h"
#include "Identifier.h"

namespace Rodin::Alert
{
  template <class T>
  class ClassException : public Exception
  {
    public:
      ClassException(const T&)
      {
        *this << "In class " << Identifier::Class(
            boost::typeindex::type_id_with_cvr<T>().pretty_name()) << ". ";
      }
  };
}

#endif

