/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_EXCEPTIONS_H
#define RODIN_VARIATIONAL_EXCEPTIONS_H

#include <set>

#include "Rodin/Alert/Exception.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
   class RangeShapeMismatchException : public Alert::Exception
   {
      public:
         RangeShapeMismatchException(
               const RangeShape& a, const RangeShape& b)
         {
            *this << "RangeShape mismatch " << a << " != " << b;
         }
   };

   class IncompatibleShapeException : public Alert::Exception
   {
      public:
         IncompatibleShapeException(
               const RangeShape& a, const RangeShape& b)
         {
            *this << "Incompatible RangeShape " << a << " and " << b;
         }
   };

   class NotSquareRangeShapeException : public Alert::Exception
   {
      public:
         NotSquareRangeShapeException(const RangeShape& r)
         {
            *this << "RangeShape " << r << " must be square.";
         }
   };

   class UnexpectedRangeTypeException : public Alert::Exception
   {
      public:
         UnexpectedRangeTypeException(
               const RangeType& expected, const RangeType& actual)
            : UnexpectedRangeTypeException(std::set<RangeType>{expected}, actual)
         {}

         UnexpectedRangeTypeException(
               const std::set<RangeType>& expected, const RangeType& actual)
         {
            *this << "Unexpected \"" << actual << "\" RangeType. "
                  << "Expected one of the following: [";
            for (const auto& r : expected)
               *this << r << ", ";
            *this << "]";
         }
   };
}

#endif
