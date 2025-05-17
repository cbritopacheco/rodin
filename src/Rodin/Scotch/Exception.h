#ifndef RODIN_SCOTCH_EXCEPTION_H
#define RODIN_SCOTCH_EXCEPTION_H

#include "Rodin/Alert/MemberFunctionException.h"

namespace Rodin::Scotch
{
  template <class T, class FuncName>
  class Exception : public Alert::MemberFunctionException<T, FuncName>
  {
    public:
      /// Alias for the parent exception class.
      using Parent = Alert::MemberFunctionException<T, FuncName>;

      /**
       * @brief Constructs an UnsafeAccessException.
       *
       * Initializes the exception with the given object and function name.
       *
       * @param cls The object (resource) being accessed unsafely.
       * @param funcName The name of the member function where the unsafe
       * access occurred.
       */
      Exception(const char* msg, const T& cls, const FuncName& funcName)
        : Parent(cls, funcName)
      {
        *this << Alert::Text(msg).setBold().setUnderline() << ".";
      }

      virtual void raise() const override
      {
        Parent::raise();
      }
  };
}

#endif
