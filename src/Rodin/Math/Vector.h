/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MATH_VECTOR_H
#define RODIN_MATH_VECTOR_H

#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "Rodin/Types.h"

namespace Rodin::Math
{
  using Vector = Eigen::VectorX<Scalar>;

  template <size_t Size>
  using FixedSizeVector = Eigen::Vector<Scalar, Size>;

  using Vector2 = FixedSizeVector<2>;
  using Vector3 = FixedSizeVector<3>;
  using Vector4 = FixedSizeVector<4>;
  using Vector8 = FixedSizeVector<8>;
  using Vector16 = FixedSizeVector<16>;
  using Vector32 = FixedSizeVector<32>;
  using Vector128 = FixedSizeVector<128>;
}

#endif
