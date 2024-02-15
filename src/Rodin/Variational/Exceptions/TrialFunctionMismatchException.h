/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_EXCEPTIONS_TRIALFUNCTIONMISMATCHEXCEPTION_H
#define RODIN_VARIATIONAL_EXCEPTIONS_TRIALFUNCTIONMISMATCHEXCEPTION_H

#include "Rodin/Alert/Exception.h"
#include "Rodin/Alert/Identifier.h"

namespace Rodin::Variational
{
  class TrialFunctionMismatchException : public Alert::Exception
  {
    public:
      using Parent = Alert::Exception;

      template <class UPb>
      TrialFunctionMismatchException(const UPb& uPb)
      {
        *this << "Bad "
              << Alert::Identifier::Class("Problem")
              << " definition. Mismatching "
              << Alert::Identifier::Class("TrialFunction")
              << " in the problem body: ";
        if (uPb.getName())
        {
          *this << Alert::Notation("\"")
                << Alert::Notation(uPb.getName())
                << Alert::Notation("\"");
        }
        else
        {
          *this << Alert::Notation("TrialFunction[ UUID = ")
                << Alert::Notation(uPb.getUUID())
                << Alert::Notation(" ]");
        }
        *this << " does not appear in the declaration.";
      }
  };
}

#endif




