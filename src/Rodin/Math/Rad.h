/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Rad.h
 * @brief Radian angle unit type.
 *
 * This file defines the Rad class for representing angles in radians with
 * type-safe unit semantics.
 */
#ifndef RODIN_MATH_RAD_H
#define RODIN_MATH_RAD_H

#include "Rodin/Types.h"
#include "Unit.h"

namespace Rodin::Math
{
  /**
   * @brief Represents an angle in radians.
   *
   * This class provides a type-safe wrapper for angles measured in radians,
   * inheriting all the unit arithmetic operations from the Unit base class.
   *
   * ## Example Usage
   * ```cpp
   * Rad angle1(Math::Constants::pi());  // π radians (180 degrees)
   * Rad angle2(Math::Constants::pi() / 2);  // π/2 radians (90 degrees)
   * Rad sum = angle1 + angle2;  // 3π/2 radians
   * ```
   *
   * @see Unit
   */
  class Rad : public Unit<Rad, Real>
  {
    public:
      using Parent = Unit<Rad, Real>;
      using Parent::Parent;
  };
}

#endif
