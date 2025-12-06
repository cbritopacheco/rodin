/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_MEMBERFUNCTIONwarning_H
#define RODIN_ALERT_MEMBERFUNCTIONwarning_H

#include <boost/type_index.hpp>
#include <boost/current_function.hpp>

#include "Warning.h"
#include "Identifier.h"

namespace Rodin::Alert
{
  /**
   * @brief Warning type for issues occurring in member functions.
   * @tparam T The type of the class where the member function resides.
   * @tparam FuncName The type used to represent the function name.
   *
   * This warning class is used to report issues in member functions,
   * providing context about the class and function where the warning
   * originated. Unlike MemberFunctionException, this does not terminate
   * the program but instead displays a formatted warning message.
   *
   * Example usage:
   * @code{.cpp}
   * MemberFunctionWarning(*this, __func__) << "Deprecated feature used" << Raise;
   * @endcode
   *
   * This produces output like:
   * "Warning: In member function myFunction of class MyClass: Deprecated feature used"
   */
  template <class T, class FuncName>
  class MemberFunctionWarning : public Warning
  {
    public:
      /// @brief Alias for the parent warning type.
      using Parent = Warning;

      /**
       * @brief Constructs a MemberFunctionWarning.
       * @param funcName The name of the member function where the warning occurred.
       *
       * This constructor formats a warning message that includes the name of
       * the member function and the class where the warning was generated.
       * It uses Boost.TypeIndex to automatically extract the class name.
       */
      MemberFunctionWarning(const T&, const FuncName& funcName)
      {
        const auto& className = boost::typeindex::type_id_with_cvr<T>().pretty_name();
        *this << "In member function " << Identifier::Function(funcName)
              << " of class " << Identifier::Class(className) << ": ";
      }
  };
}

#endif


