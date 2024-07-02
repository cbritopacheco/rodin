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
  Real pi()
  {
    return M_PI;
  }

  inline
  constexpr
  Real nan()
  {
    return NAN;
  }

  inline
  constexpr
  Real epsilon()
  {
    return std::numeric_limits<Real>::epsilon();
  }

  inline
  constexpr
  Real zero()
  {
    return Real(0);
  }

  inline
  constexpr
  Boolean isZero(Real x)
  {
    return x == Real(0);
  }

  inline
  constexpr
  Real one()
  {
    return Real(1);
  }

  inline
  constexpr
  Boolean isOne(Real x)
  {
    return x == Real(1);
  }
}

#endif
