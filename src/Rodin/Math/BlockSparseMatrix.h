/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MATH_BLOCKSPARSEMATRIX_H
#define RODIN_MATH_BLOCKSPARSEMATRIX_H

#include <Eigen/Core>
#include <unsupported/Eigen/SparseExtra>
#include <unsupported/Eigen/src/SparseExtra/BlockSparseMatrix.h>

#include "Rodin/Types.h"

#include "ForwardDecls.h"

namespace Rodin::Math
{
  /**
   * @brief Sparse matrix type
   */
  using BlockSparseMatrix = Eigen::BlockSparseMatrix<Scalar>;
}

#endif


