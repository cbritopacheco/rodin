/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_CORE_CONSTANTS_H
#define RODIN_CORE_CONSTANTS_H

#include <cmath>
#include <type_traits>

#include "Types.h"

namespace Rodin::Math::Constants
{
  /**
   * @brief Computes the number @f$ \pi @f$ to machine precision.
   */
  inline
  constexpr
  Scalar pi()
  {
    return M_PI;
  }

  inline
  constexpr
  Scalar nan()
  {
    return NAN;
  }

  inline
  constexpr
  Scalar epsilon()
  {
    return std::numeric_limits<Scalar>::epsilon();
  }

  inline
  constexpr
  Scalar zero()
  {
    return Scalar(0);
  }

  inline
  constexpr
  Boolean isZero(Scalar x)
  {
    return x == Scalar(0);
  }

  inline
  constexpr
  Scalar one()
  {
    return Scalar(1);
  }

  inline
  constexpr
  Boolean isOne(Scalar x)
  {
    return x == Scalar(1);
  }
}

#endif