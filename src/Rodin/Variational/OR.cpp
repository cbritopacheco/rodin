/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "OR.h"

namespace Rodin::Variational
{
   OR<BooleanFunctionBase, BooleanFunctionBase>
   operator||(const BooleanFunctionBase& lhs, const BooleanFunctionBase& rhs)
   {
      return OR(lhs, rhs);
   }

   OR<BooleanFunctionBase, BooleanFunctionBase>
   operator||(bool lhs, const BooleanFunctionBase& rhs)
   {
      return OR(BooleanFunction(lhs), rhs);
   }

   OR<BooleanFunctionBase, BooleanFunctionBase>
   operator||(const BooleanFunctionBase& lhs, bool rhs)
   {
      return OR(lhs, BooleanFunction(rhs));
   }
}


