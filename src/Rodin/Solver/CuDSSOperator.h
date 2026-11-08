/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file CuDSSOperator.h
 * @brief Compressed arrays of a host sparse matrix for the GPU direct solvers.
 *
 * The specialization needs no GPU header, so that CuDSS and RocSolverRF
 * share it.
 */
#ifndef RODIN_SOLVER_CUDSSOPERATOR_H
#define RODIN_SOLVER_CUDSSOPERATOR_H

#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"

namespace Rodin::Solver::Internal
{
  /// @brief Host CSC arrays, read as the CSR arrays of the transpose.
  template <>
  struct CuDSSOperator<Math::SparseMatrix<Real>>
  {
    /// @brief Matrix type.
    using MatrixType = Math::SparseMatrix<Real>;

    /// @brief Whether the arrays are those of the transpose.
    static constexpr bool Transposed = true;

    /// @brief Whether the arrays are in device memory.
    static constexpr bool Resident = false;

    /**
     * @brief Compresses the matrix.
     * @param matrix Host sparse matrix, compressed in place.
     */
    static void prepare(MatrixType& matrix)
    {
      matrix.makeCompressed();
    }

    /**
     * @brief Number of rows.
     * @param matrix Compressed host sparse matrix.
     * @return Row count of @p matrix.
     */
    static size_t getRowCount(const MatrixType& matrix)
    {
      return static_cast<size_t>(matrix.rows());
    }

    /**
     * @brief Number of columns.
     * @param matrix Compressed host sparse matrix.
     * @return Column count of @p matrix.
     */
    static size_t getColumnCount(const MatrixType& matrix)
    {
      return static_cast<size_t>(matrix.cols());
    }

    /**
     * @brief Number of stored coefficients.
     * @param matrix Compressed host sparse matrix.
     * @return Stored coefficient count of @p matrix.
     */
    static size_t getNonZeroCount(const MatrixType& matrix)
    {
      return static_cast<size_t>(matrix.nonZeros());
    }

    /**
     * @brief Start of every compressed row, followed by the coefficient count.
     * @param matrix Compressed host sparse matrix.
     * @return Pointer to the outer index array of @p matrix.
     */
    static const Math::SparseIndex* getStarts(const MatrixType& matrix)
    {
      return matrix.outerIndexPtr();
    }

    /**
     * @brief Index of every coefficient in its compressed row.
     * @param matrix Compressed host sparse matrix.
     * @return Pointer to the inner index array of @p matrix.
     */
    static const Math::SparseIndex* getIndices(const MatrixType& matrix)
    {
      return matrix.innerIndexPtr();
    }

    /**
     * @brief Value of every coefficient.
     * @param matrix Compressed host sparse matrix.
     * @return Pointer to the coefficient array of @p matrix.
     */
    static const Real* getValues(const MatrixType& matrix)
    {
      return matrix.valuePtr();
    }
  };
}

#endif
