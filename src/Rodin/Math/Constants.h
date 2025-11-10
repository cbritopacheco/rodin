/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Constants.h
 * @brief Mathematical constants and related utility functions.
 *
 * This file provides commonly used mathematical constants with machine precision,
 * as well as utility functions for checking special values.
 */
#ifndef RODIN_CORE_CONSTANTS_H
#define RODIN_CORE_CONSTANTS_H

#include <cmath>
#include <limits>
#include <type_traits>

#include "Types.h"

namespace Rodin::Math::Constants
{
  /**
   * @brief Returns the mathematical constant @f$ \pi @f$ to machine precision.
   *
   * The value of @f$ \pi \approx 3.14159265358979323846 @f$ is the ratio of a
   * circle's circumference to its diameter.
   *
   * @return The value of @f$ \pi @f$ as a Real
   */
  inline
  constexpr
  Real pi()
  {
    return M_PI;
  }

  /**
   * @brief Returns the machine epsilon for Real type.
   *
   * Machine epsilon is the smallest positive number @f$ \varepsilon @f$ such that
   * @f$ 1 + \varepsilon \neq 1 @f$ in floating-point arithmetic. It represents
   * the upper bound on relative rounding error.
   *
   * @return The machine epsilon value
   */
  inline
  constexpr
  Real epsilon()
  {
    return std::numeric_limits<Real>::epsilon();
  }

  /**
   * @brief Returns the value zero.
   *
   * @return The Real value @f$ 0 @f$
   */
  inline
  constexpr
  Real zero()
  {
    return Real(0);
  }

  /**
   * @brief Checks if a value is exactly zero.
   *
   * Performs an exact comparison with zero. For approximate comparisons
   * near zero, consider using a tolerance based on epsilon().
   *
   * @param[in] x Value to check
   * @return True if @f$ x = 0 @f$, false otherwise
   */
  inline
  constexpr
  Boolean isZero(Real x)
  {
    return x == Real(0);
  }

  /**
   * @brief Returns the value one.
   *
   * @return The Real value @f$ 1 @f$
   */
  inline
  constexpr
  Real one()
  {
    return Real(1);
  }

  /**
   * @brief Checks if a value is exactly one.
   *
   * Performs an exact comparison with one. For approximate comparisons,
   * consider using a tolerance based on epsilon().
   *
   * @param[in] x Value to check
   * @return True if @f$ x = 1 @f$, false otherwise
   */
  inline
  constexpr
  Boolean isOne(Real x)
  {
    return x == Real(1);
  }
}

#endif
