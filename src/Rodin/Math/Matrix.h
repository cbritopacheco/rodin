/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Matrix.h
 * @brief Dense matrix type aliases and definitions.
 *
 * This file provides type aliases for various dense matrix types built on Eigen.
 * Matrices are stored in column-major order by default and support both dynamic
 * and fixed sizes.
 */
#ifndef RODIN_MATH_MATRIX_H
#define RODIN_MATH_MATRIX_H

#include <Eigen/Core>
#include <Eigen/Dense>

#include "Rodin/Configure.h"

#include "Rodin/FormLanguage/Traits.h"

#include "Rodin/Types.h"

namespace Rodin::Math
{
  /**
   * @brief Dynamic-size dense matrix type.
   *
   * A dense matrix with both dimensions determined at runtime. The matrix
   * is stored in column-major order (Fortran style) and supports standard
   * linear algebra operations.
   *
   * @tparam ScalarType The element type (e.g., Real, Complex)
   */
  template <class ScalarType>
  using Matrix = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;

  /**
   * @brief Real-valued dense matrix.
   *
   * Convenience alias for Matrix<Real>.
   */
  using RealMatrix = Matrix<Real>;

  /**
   * @brief Complex-valued dense matrix.
   *
   * Convenience alias for Matrix<Complex>, commonly used in eigenvalue
   * problems, frequency-domain analysis, and other applications requiring
   * complex arithmetic.
   */
  using ComplexMatrix = Matrix<Complex>;

  /**
   * @brief Spatial matrix with bounded maximum dimensions.
   *
   * A dynamic-size matrix with maximum dimensions bounded by
   * RODIN_MAXIMAL_SPACE_DIMENSION. Used for geometric transformations,
   * Jacobians, and other spatial operators to optimize memory allocation.
   *
   * Typical uses:
   * - Jacobian matrices: @f$ J = \frac{\partial \mathbf{f}}{\partial \mathbf{x}} @f$
   * - Metric tensors
   * - Local coordinate transformations
   *
   * @tparam ScalarType The element type
   */
  template <class ScalarType>
  using SpatialMatrix = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, 0,
        RODIN_MAXIMAL_SPACE_DIMENSION, RODIN_MAXIMAL_SPACE_DIMENSION>;

  /**
   * @brief Point matrix with bounded row dimension.
   *
   * A dynamic-size matrix where the number of rows is bounded by
   * RODIN_MAXIMAL_SPACE_DIMENSION. Commonly used to store collections of
   * spatial points as columns.
   *
   * Example: storing vertex coordinates of a mesh element where each column
   * represents one vertex's coordinates.
   */
  using PointMatrix =
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, 0, RODIN_MAXIMAL_SPACE_DIMENSION, Eigen::Dynamic>;

  /**
   * @brief Fixed-size matrix type.
   *
   * A compile-time fixed-size matrix. Both dimensions are known at compile time,
   * enabling optimizations and stack allocation.
   *
   * Example uses:
   * - 2×2 rotation matrices
   * - 3×3 transformation matrices
   * - Small dense blocks in larger systems
   *
   * @tparam ScalarType The element type
   * @tparam Rows The number of rows (must be known at compile time)
   * @tparam Cols The number of columns (must be known at compile time)
   */
  template <class ScalarType, size_t Rows, size_t Cols>
  using FixedSizeMatrix = Eigen::Matrix<ScalarType, Rows, Cols>;
}

namespace Rodin::FormLanguage
{
  template <class Number>
  struct Traits<Math::Matrix<Number>>
  {
    using ScalarType = Number;
  };

  template <class Number>
  struct Traits<Math::SpatialMatrix<Number>>
  {
    using ScalarType = Number;
  };
}

#endif

