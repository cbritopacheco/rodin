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
  /**
   * @brief Exception class that automatically includes the namespace name.
   * @ingroup AlertModule
   *
   * Specialized exception class that prepends a namespace identifier to the
   * error message. The namespace name is formatted with appropriate coloring
   * and emphasis for clear visual identification in error output.
   *
   * Example usage:
   * @code{.cpp}
   * NamespacedException("Rodin::Geometry") << "Invalid operation" << Raise;
   * @endcode
   *
   * This produces output like:
   * "Error: Rodin::Geometry. Invalid operation"
   */
  class NamespacedException : public Exception
  {
    public:
      /**
       * @brief Constructs a NamespacedException with the given namespace name.
       * @param namesp The namespace name to include in the error message.
       *
       * Automatically formats the namespace name with appropriate styling and
       * prepends it to the exception message.
       */
      NamespacedException(const std::string& namesp)
      {
        *this << Identifier::Namespace(namesp) << ". ";
      }
  };
}

#endif
