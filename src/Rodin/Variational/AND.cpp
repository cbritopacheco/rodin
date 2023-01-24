/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "AND.h"

namespace Rodin::Variational
{
  AND<BooleanFunctionBase, BooleanFunctionBase>
  operator&&(const BooleanFunctionBase& lhs, const BooleanFunctionBase& rhs)
  {
    return AND(lhs, rhs);
  }

  AND<BooleanFunctionBase, BooleanFunctionBase>
  operator&&(bool lhs, const BooleanFunctionBase& rhs)
  {
    return AND(BooleanFunction(lhs), rhs);
  }

  AND<BooleanFunctionBase, BooleanFunctionBase>
  operator&&(const BooleanFunctionBase& lhs, bool rhs)
  {
    return AND(lhs, BooleanFunction(rhs));
  }
}

