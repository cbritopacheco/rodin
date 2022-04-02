/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Division.h"

namespace Rodin::Variational
{
   Division<VectorFunctionBase, ScalarFunctionBase>
   operator/(const VectorFunctionBase& lhs, const ScalarFunctionBase& rhs)
   {
      return Division(lhs, rhs);
   }
}

