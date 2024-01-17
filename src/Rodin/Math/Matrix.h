/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MATH_MATRIX_H
#define RODIN_MATH_MATRIX_H

#include <Eigen/Core>
#include <Eigen/Dense>

#include "Rodin/Configure.h"

#include "Rodin/Types.h"
#include "Rodin/Array.h"

namespace Rodin::Math
{
  /**
   * @brief Dense matrix type
   */
  using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

  /**
   * @brief Spatial matrix
   */
  using SpatialMatrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, 0,
        RODIN_MAXIMAL_SPACE_DIMENSION, RODIN_MAXIMAL_SPACE_DIMENSION>;

  /**
   * @brief Point matrix
   */
  using PointMatrix =
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, 0, RODIN_MAXIMAL_SPACE_DIMENSION, Eigen::Dynamic>;

  /**
   * @brief Represents a fixed size matrix
   */
  template <size_t Rows, size_t Cols>
  using FixedSizeMatrix = Eigen::Matrix<Scalar, Rows, Cols>;
}

#endif

