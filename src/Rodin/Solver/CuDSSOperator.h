/*
 *          Copyright Jean-Baptiste KECK 2026.
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

#include <algorithm>
#include <vector>

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
     * @brief Owned compressed pattern, compared exactly in host memory.
     *
     * Values and array addresses do not participate in equality. The snapshot
     * owns its indices so that in-place edits of the matrix are detected.
     */
    class Pattern
    {
      public:
        /**
         * @brief Retains the compressed pattern of a prepared matrix.
         * @param matrix Compressed host matrix.
         */
        void assign(const MatrixType& matrix)
        {
          m_rows = getRowCount(matrix);
          m_columns = getColumnCount(matrix);
          const size_t nonzeros = getNonZeroCount(matrix);
          m_indices.resize(m_rows + 1 + nonzeros);
          std::copy_n(getStarts(matrix), m_rows + 1, m_indices.data());
          std::copy_n(getIndices(matrix), nonzeros, m_indices.data() + m_rows + 1);
        }

        /**
         * @brief Compares all offsets and indices, without comparing values.
         * @param matrix Compressed host matrix.
         * @return Whether the matrix has the retained pattern.
         */
        bool matches(const MatrixType& matrix) const
        {
          const size_t nonzeros = getNonZeroCount(matrix);
          return m_rows == getRowCount(matrix) && m_columns == getColumnCount(matrix) &&
            m_indices.size() == m_rows + 1 + nonzeros &&
            std::equal(
              getStarts(matrix), getStarts(matrix) + m_rows + 1, m_indices.begin()) &&
            std::equal(getIndices(matrix), getIndices(matrix) + nonzeros,
              m_indices.begin() + m_rows + 1);
        }

        /// @brief Releases the retained pattern.
        void clear()
        {
          m_indices = {};
          m_indices.shrink_to_fit();
          m_rows = 0;
          m_columns = 0;
        }

      private:
        size_t m_rows = 0;
        size_t m_columns = 0;
        std::vector<Math::SparseIndex> m_indices;
    };

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
