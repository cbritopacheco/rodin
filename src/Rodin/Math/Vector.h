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

#include "Rodin/FormLanguage/Traits.h"

namespace Rodin::Math
{
  /**
   * @brief Dense vector type.
   */
  template <class ScalarType>
  using Vector = Eigen::VectorX<ScalarType>;

  using ComplexVector = Vector<Complex>;

  template <class ScalarType>
  using SpatialVector =
    Eigen::Matrix<ScalarType, Eigen::Dynamic, 1, 0, RODIN_MAXIMAL_SPACE_DIMENSION, 1>;

  using PointVector = SpatialVector<Real>;

  template <class ScalarType, size_t Size>
  using FixedSizeVector = Eigen::Vector<ScalarType, Size>;

  template <class ScalarType>
  using Vector2 = FixedSizeVector<ScalarType, 2>;

  template <class ScalarType>
  using Vector3 = FixedSizeVector<ScalarType, 3>;

  template <class ScalarType>
  using Vector4 = FixedSizeVector<ScalarType, 4>;

  template <class ScalarType>
  using Vector8 = FixedSizeVector<ScalarType, 8>;

  template <class ScalarType>
  using Vector16 = FixedSizeVector<ScalarType, 16>;

  template <class ScalarType>
  using Vector32 = FixedSizeVector<ScalarType, 32>;

  template <class ScalarType>
  using Vector64 = FixedSizeVector<ScalarType, 64>;

  template <class ScalarType>
  using Vector128 = FixedSizeVector<ScalarType, 128>;
}

namespace Rodin::FormLanguage
{
  template <class Number>
  struct Traits<Math::Vector<Number>>
  {
    using ScalarType = Number;
  };

  template <class Number>
  struct Traits<Math::SpatialVector<Number>>
  {
    using ScalarType = Number;
  };

  template <class Number, size_t S>
  struct Traits<Math::FixedSizeVector<Number, S>>
  {
    using ScalarType = Number;
    static constexpr size_t Size = S;
  };
}

#endif
