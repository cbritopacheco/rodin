/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_EXCEPTIONS_UNDETERMINEDTRACEDOMAINEXCEPTION_H
#define RODIN_VARIATIONAL_EXCEPTIONS_UNDETERMINEDTRACEDOMAINEXCEPTION_H

#include <boost/type_index.hpp>
#include <boost/current_function.hpp>

#include "Rodin/Alert/MemberFunctionException.h"

namespace Rodin::Variational
{
  template <class T, class FuncName>
  class UndeterminedTraceDomainException : public Alert::MemberFunctionException<T, FuncName>
  {
    public:
      using Parent = Alert::MemberFunctionException<T, FuncName>;

      template <class Iterator>
      UndeterminedTraceDomainException(const T& cls, const FuncName& funcName,
          const std::pair<size_t, Index>& p, Iterator begin, Iterator end)
        : Parent(cls, funcName)
      {
        *this << "Could not determine the trace polytope for the interface "
              << Alert::Notation::Polytope(p.first, p.second)
              << " with the provided trace domain "
              << Alert::Notation::Set(begin, end)
              << '.'
              << Alert::Raise;
      }
  };

  template <class T, class FuncName, class Iterator>
  UndeterminedTraceDomainException(const T&, const FuncName&, const std::pair<size_t, Index>&, Iterator, Iterator)
    -> UndeterminedTraceDomainException<T, FuncName>;

}

#endif



