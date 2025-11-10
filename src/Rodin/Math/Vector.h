/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Vector.h
 * @brief Dense vector type aliases and definitions.
 *
 * This file provides type aliases for various vector types built on Eigen.
 * Vectors are column-oriented by default and support both dynamic and fixed sizes.
 */
#ifndef RODIN_MATH_VECTOR_H
#define RODIN_MATH_VECTOR_H

#include <Eigen/Core>
#include <Eigen/Dense>

#include "Rodin/Types.h"

#include "Rodin/FormLanguage/Traits.h"

namespace Rodin::Math
{
  /**
   * @brief Dynamic-size dense vector type.
   *
   * A column vector with dynamic size determined at runtime. The vector
   * supports standard linear algebra operations.
   *
   * @tparam ScalarType The element type (e.g., Real, Complex)
   */
  template <class ScalarType>
  using Vector = Eigen::VectorX<ScalarType>;

  /**
   * @brief Dynamic-size complex-valued vector.
   *
   * Convenience alias for Vector<Complex>.
   */
  using ComplexVector = Vector<Complex>;

  /**
   * @brief Spatial vector with bounded maximum size.
   *
   * A dynamic-size vector with maximum size bounded by RODIN_MAXIMAL_SPACE_DIMENSION.
   * Used for geometric quantities in 2D or 3D space to optimize memory allocation.
   *
   * @tparam ScalarType The element type
   */
  template <class ScalarType>
  using SpatialVector =
    Eigen::Matrix<ScalarType, Eigen::Dynamic, 1, 0, RODIN_MAXIMAL_SPACE_DIMENSION, 1>;

  /**
   * @brief Real-valued spatial vector for point coordinates.
   *
   * Convenience alias for SpatialVector<Real>, commonly used to represent
   * points in 2D or 3D space.
   */
  using SpatialPoint = SpatialVector<Real>;

  /**
   * @brief Fixed-size vector type.
   *
   * A compile-time fixed-size column vector. The size is known at compile time,
   * enabling optimizations and stack allocation.
   *
   * @tparam ScalarType The element type
   * @tparam Size The number of elements (must be known at compile time)
   */
  template <class ScalarType, size_t Size>
  using FixedSizeVector = Eigen::Vector<ScalarType, Size>;

  /**
   * @brief 2D fixed-size vector.
   *
   * A vector with exactly 2 elements.
   *
   * @tparam ScalarType The element type
   */
  template <class ScalarType>
  using Vector2 = FixedSizeVector<ScalarType, 2>;

  /**
   * @brief 3D fixed-size vector.
   *
   * A vector with exactly 3 elements.
   *
   * @tparam ScalarType The element type
   */
  template <class ScalarType>
  using Vector3 = FixedSizeVector<ScalarType, 3>;

  /**
   * @brief 4D fixed-size vector.
   *
   * A vector with exactly 4 elements.
   *
   * @tparam ScalarType The element type
   */
  template <class ScalarType>
  using Vector4 = FixedSizeVector<ScalarType, 4>;

  /**
   * @brief 8D fixed-size vector.
   *
   * A vector with exactly 8 elements.
   *
   * @tparam ScalarType The element type
   */
  template <class ScalarType>
  using Vector8 = FixedSizeVector<ScalarType, 8>;

  /**
   * @brief 16D fixed-size vector.
   *
   * A vector with exactly 16 elements.
   *
   * @tparam ScalarType The element type
   */
  template <class ScalarType>
  using Vector16 = FixedSizeVector<ScalarType, 16>;

  /**
   * @brief 32D fixed-size vector.
   *
   * A vector with exactly 32 elements.
   *
   * @tparam ScalarType The element type
   */
  template <class ScalarType>
  using Vector32 = FixedSizeVector<ScalarType, 32>;

  /**
   * @brief 64D fixed-size vector.
   *
   * A vector with exactly 64 elements.
   *
   * @tparam ScalarType The element type
   */
  template <class ScalarType>
  using Vector64 = FixedSizeVector<ScalarType, 64>;

  /**
   * @brief 128D fixed-size vector.
   *
   * A vector with exactly 128 elements.
   *
   * @tparam ScalarType The element type
   */
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
