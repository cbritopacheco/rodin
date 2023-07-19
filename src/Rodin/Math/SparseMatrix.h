/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MATH_SPARSEMATRIX_H
#define RODIN_MATH_SPARSEMATRIX_H

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "Rodin/Types.h"
#include "Rodin/Math/Matrix.h"

namespace Rodin::Math
{
  /**
   * @brief Sparse matrix type
   */
  using SparseMatrix = Eigen::SparseMatrix<Scalar>;
}

#endif

