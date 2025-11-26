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
  /**
   * @brief Exception class that automatically includes the class name.
   * @ingroup AlertModule
   * @tparam T The class type where the exception occurred.
   *
   * Specialized exception class that prepends the class name to the
   * error message. Uses Boost.TypeIndex to automatically extract and
   * format the class name with proper coloring and emphasis.
   *
   * Example usage:
   * @code{.cpp}
   * ClassException(*this) << "Invalid state detected" << Raise;
   * @endcode
   *
   * This produces output like:
   * "Error: In class MyClass. Invalid state detected"
   */
  template <class T>
  class ClassException : public Exception
  {
    public:
      /**
       * @brief Constructs a ClassException for the given class instance.
       * @param An instance of the class (used for type deduction).
       *
       * Automatically extracts the class name using Boost.TypeIndex and
       * prepends it to the exception message with appropriate formatting.
       */
      ClassException(const T&)
      {
        *this << "In class " << Identifier::Class(
            boost::typeindex::type_id_with_cvr<T>().pretty_name()) << ". ";
      }
  };
}

#endif

