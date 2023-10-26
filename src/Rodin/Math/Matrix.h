/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MATH_DENSEMATRIX_H
#define RODIN_MATH_DENSEMATRIX_H

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
  using Matrix = Eigen::MatrixX<Scalar>;

  using SpatialMatrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, 0,
        RODIN_MAXIMAL_SPACE_DIMENSION, Eigen::Dynamic>;

  using PointMatrix =
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, 0, RODIN_MAXIMAL_SPACE_DIMENSION, Eigen::Dynamic>;

  template <size_t Rows, size_t Cols>
  using FixedSizeMatrix = Eigen::Matrix<Scalar, Rows, Cols>;
}

#endif

