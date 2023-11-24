/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_F_H
#define RODIN_VARIATIONAL_F_H

#include "ScalarFunction.h"

/**
 * @brief Contains built-in functions.
 */
namespace Rodin::Variational::F
{
  /**
   * @brief Represents the first coordinate of the point.
   *
   * Represents the function @f$ f : \mathcal{T}_h \rightarrow \mathbb{R} @f$
   * defined by:
   * @f[
   *   f(x_1, \cdots, x_d) = x_1 \: .
   * @f]
   */
  static ScalarFunction x([](const Geometry::Point& p) { return p.x(); });

  /**
   * @brief Represents the first coordinate of the point.
   *
   * Represents the function @f$ f : \mathcal{T}_h \rightarrow \mathbb{R} @f$
   * defined by:
   * @f[
   *   f(x_1, \cdots, x_d) = x_2 \: .
   * @f]
   */
  static ScalarFunction y([](const Geometry::Point& p) { return p.y(); });

  /**
   * @brief Represents the first coordinate of the point.
   *
   * Represents the function @f$ f : \mathcal{T}_h \rightarrow \mathbb{R} @f$
   * defined by:
   * @f[
   *   f(x_1, \cdots, x_d) = x_3 \: .
   * @f]
   */
  static ScalarFunction z([](const Geometry::Point& p) { return p.z(); });
}

#endif
