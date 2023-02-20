/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_RANGETYPE_H
#define RODIN_VARIATIONAL_RANGETYPE_H

#include <iostream>
#include "ForwardDecls.h"

namespace Rodin::Variational
{
  enum class RangeType
  {
    Boolean,
    Scalar,
    Vector,
    Matrix
  };

  std::ostream& operator<<(std::ostream& os, const RangeType& obj);
}

#endif
