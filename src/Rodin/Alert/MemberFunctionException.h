/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_MEMBERFUNCTIONEXCEPTION_H
#define RODIN_ALERT_MEMBERFUNCTIONEXCEPTION_H

#include <boost/type_index.hpp>
#include <boost/current_function.hpp>

#include "Exception.h"
#include "Identifier.h"

namespace Rodin::Alert
{
  /**
   * @brief Exception type for errors occurring in member functions.
   *
   * This exception class is used to capture errors in member functions,
   * providing context about the class and function in which the exception
   * occurred.
   *
   * @tparam T The type of the class where the member function resides.
   * @tparam FuncName The type used to represent the function name.
   */
  template <class T, class FuncName>
  class MemberFunctionException : public Exception
  {
    public:
      /// Alias for the parent exception type.
      using Parent = Exception;

      /**
       * @brief Constructs a MemberFunctionException.
       *
       * This constructor formats an error message that includes the name of
       * the member function and the class where the exception occurred.
       *
       * @param cls An instance of the class where the error occurred. (Unused
       * directly, but used to deduce the class type.)
       * @param funcName The name of the member function in which the error
       * occurred.
       */
      MemberFunctionException(const T&, const FuncName& funcName)
      {
        const auto& className = boost::typeindex::type_id_with_cvr<T>().pretty_name();
        *this << "In member function " << Identifier::Function(funcName)
              << " of class " << Identifier::Class(className) << ": ";
      }

      /**
       * @brief Raises the exception.
       *
       * Calls the parent exception's raise() method to propagate the
       * exception.
       */
      virtual void raise() const override
      {
        Parent::raise();
      }
  };
}

#endif


