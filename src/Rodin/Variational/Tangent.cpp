/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Tangent.h"

namespace Rodin::Variational
{
  Tangent<FunctionBase> tan(const FunctionBase& op)
  {
    return Tangent<FunctionBase>(op);
  }
}
