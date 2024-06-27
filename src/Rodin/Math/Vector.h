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
  /**
   * @brief Dense vector type.
   */
  template <class NumberType>
  using Vector = Eigen::VectorX<NumberType>;

  using ComplexVector = Vector<Complex>;

  template <class NumberType>
  using SpatialVector =
    Eigen::Matrix<NumberType, Eigen::Dynamic, 1, 0, RODIN_MAXIMAL_SPACE_DIMENSION, 1>;

  template <class NumberType, size_t Size>
  using FixedSizeVector = Eigen::Vector<NumberType, Size>;

  template <class NumberType>
  using Vector2 = FixedSizeVector<NumberType, 2>;

  template <class NumberType>
  using Vector3 = FixedSizeVector<NumberType, 3>;

  template <class NumberType>
  using Vector4 = FixedSizeVector<NumberType, 4>;

  template <class NumberType>
  using Vector8 = FixedSizeVector<NumberType, 8>;

  template <class NumberType>
  using Vector16 = FixedSizeVector<NumberType, 16>;

  template <class NumberType>
  using Vector32 = FixedSizeVector<NumberType, 32>;

  template <class NumberType>
  using Vector64 = FixedSizeVector<NumberType, 64>;

  template <class NumberType>
  using Vector128 = FixedSizeVector<NumberType, 128>;
}

#endif
