/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "LEQ.h"

namespace Rodin::Variational
{
   LEQ<FunctionBase, FunctionBase>
   operator<=(const FunctionBase& lhs, const FunctionBase& rhs)
   {
      return LEQ(lhs, rhs);
   }
}
