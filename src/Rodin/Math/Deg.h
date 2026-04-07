/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Deg.h
 * @brief Degree angle unit type with radian conversion.
 *
 * This file defines the Deg class for representing angles in degrees with
 * type-safe unit semantics and conversion to radians.
 */
#ifndef RODIN_MATH_DEG_H
#define RODIN_MATH_DEG_H

#include "Rodin/Types.h"
#include "Unit.h"
#include "Constants.h"
#include "Rad.h"

namespace Rodin::Math
{
  /**
   * @brief Represents an angle in degrees.
   *
   * This class provides a type-safe wrapper for angles measured in degrees,
   * inheriting all the unit arithmetic operations from the Unit base class.
   * Provides a conversion method to radians.
   *
   * ## Example Usage
   * ```cpp
   * Deg angle(180.0);  // 180 degrees
   * Rad radians = angle.toRad();  // π radians
   * ```
   *
   * @see Rad, Unit
   */
  class Deg : public Unit<Deg, Real>
  {
    public:
      using Parent = Unit<Deg, Real>;
      using Parent::Parent;

      /**
       * @brief Converts the angle from degrees to radians.
       *
       * Applies the conversion @f$ \text{rad} = \text{deg} \times \frac{\pi}{180} @f$.
       *
       * @return The angle in radians
       */
      Rad toRad() const
      {
        return Rad(Real(*this) * Constants::pi() / Real(180));
      }
  };
}

#endif
