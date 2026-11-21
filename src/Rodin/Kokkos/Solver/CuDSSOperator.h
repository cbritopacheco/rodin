/*
 *          Copyright Carlos BRITO PACHECO 2026.
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_SOLVER_CUDSSOPERATOR_H
#define RODIN_KOKKOS_SOLVER_CUDSSOPERATOR_H

/**
 * @file
 * @brief Direct-solver access to backend sparsematrix storage.
 */

#include "Rodin/Solver/CuDSSOperator.h"
#include "Rodin/Kokkos/Math/SparseMatrix.h"

namespace Rodin::Solver
{
  /**
   * @brief CSR arrays of a Kokkos::Math::SparseMatrix, read in place by the GPU
   * direct solvers when they are in device memory and copied to the device
   * otherwise.
   */
  template <class Backend>
  struct CuDSSOperator<::Rodin::Kokkos::Math::SparseMatrix<Real, Backend>>
  {
    /// @brief Matrix type.
    using MatrixType = ::Rodin::Kokkos::Math::SparseMatrix<Real, Backend>;

    /// @brief Whether the arrays are those of the transpose.
    static constexpr bool Transposed = false;

    /// @brief Whether the arrays are in device memory.
    static constexpr bool Resident = !::Kokkos::SpaceAccessibility<::Kokkos::HostSpace,
      typename MatrixType::MemorySpace>::accessible;

    /**
     * @brief Owned CSR pattern in the matrix memory space.
     *
     * ## Architecture
     * Snapshot the actual compressed graph independently of the matrix arrays.
     * Compare every row offset and column index on the backend, reducing to
     * one mismatch flag on the host. Values and allocation addresses are
     * irrelevant; in-place structural edits are detected without a hash.
     */
    class Pattern
    {
      using IndexType = typename MatrixType::IndexType;
      using Array = typename MatrixType::template Array<IndexType>;
      using ExecutionSpace = typename MatrixType::ExecutionSpace;

      public:
        /**
         * @brief Copies a compressed graph into independent backend storage.
         * @param matrix Matrix whose pattern is retained.
         */
        void assign(const MatrixType& matrix)
        {
          m_rows = matrix.getRowCount();
          m_columns = matrix.getColumnCount();
          m_starts = Array(::Kokkos::view_alloc(::Kokkos::WithoutInitializing,
                             "Rodin cuDSS pattern offsets"),
            matrix.getRowStarts().extent(0));
          m_indices = Array(::Kokkos::view_alloc(::Kokkos::WithoutInitializing,
                              "Rodin cuDSS pattern indices"),
            matrix.getColumns().extent(0));
          ::Kokkos::deep_copy(m_starts, matrix.getRowStarts());
          ::Kokkos::deep_copy(m_indices, matrix.getColumns());
        }

        /**
         * @brief Compares every CSR offset and index on the backend.
         * @param matrix Matrix whose current pattern is compared.
         * @return Whether the exact compressed graph is unchanged.
         */
        bool matches(const MatrixType& matrix) const
        {
          const auto starts = matrix.getRowStarts();
          const auto indices = matrix.getColumns();
          if (m_rows != matrix.getRowCount() || m_columns != matrix.getColumnCount() ||
            m_starts.extent(0) == 0 || m_starts.extent(0) != starts.extent(0) ||
            m_indices.extent(0) != indices.extent(0))
          {
            return false;
          }
          int changed = 0;
          ::Kokkos::parallel_reduce("Rodin cuDSS compare pattern",
            ::Kokkos::RangePolicy<ExecutionSpace, ::Kokkos::IndexType<size_t>>(
              0, starts.extent(0) + indices.extent(0)),
            Compare{starts, indices, m_starts, m_indices, starts.extent(0)},
            ::Kokkos::Max<int>(changed));
          return changed == 0;
        }

        /// @brief Releases the retained backend arrays.
        void clear()
        {
          m_starts = Array{};
          m_indices = Array{};
          m_rows = 0;
          m_columns = 0;
        }

        /// @brief Reduces exact entry mismatches to a single flag.
        struct Compare
        {
          /// @brief Current row offsets.
          Array starts;

          /// @brief Current column indices.
          typename MatrixType::View::index_type indices;

          /// @brief Retained row offsets.
          Array previousStarts;

          /// @brief Retained column indices.
          Array previousIndices;

          /// @brief Number of row offsets preceding the column indices.
          size_t offsets;

          /**
           * @param i Offset or column position in the concatenated graph.
           * @param mismatch Maximum mismatch flag for the reduction.
           */
          KOKKOS_INLINE_FUNCTION void operator()(size_t i, int& mismatch) const
          {
            const bool different = i < offsets
              ? starts(i) != previousStarts(i)
              : indices(i - offsets) != previousIndices(i - offsets);
            if (mismatch < static_cast<int>(different))
              mismatch = static_cast<int>(different);
          }
        };

      private:
        size_t m_rows = 0;
        size_t m_columns = 0;
        Array m_starts;
        Array m_indices;
    };

    /**
     * @brief Leaves the matrix unchanged, being already compressed.
     * @param matrix Compressed device matrix.
     */
    static void prepare([[maybe_unused]] MatrixType& matrix) {}

    /**
     * @brief Number of rows.
     * @param matrix Compressed device matrix.
     * @return Row count of @p matrix.
     */
    static size_t getRowCount(const MatrixType& matrix)
    {
      return matrix.getRowCount();
    }

    /**
     * @brief Number of columns.
     * @param matrix Compressed device matrix.
     * @return Column count of @p matrix.
     */
    static size_t getColumnCount(const MatrixType& matrix)
    {
      return matrix.getColumnCount();
    }

    /**
     * @brief Number of stored coefficients.
     * @param matrix Compressed device matrix.
     * @return Stored coefficient count of @p matrix.
     */
    static size_t getNonZeroCount(const MatrixType& matrix)
    {
      return matrix.getNonZeroCount();
    }

    /**
     * @brief Start of every row, followed by the coefficient count.
     * @param matrix Compressed device matrix.
     * @return Pointer to the CSR row offsets of @p matrix.
     */
    static const ::Rodin::Math::SparseIndex* getStarts(const MatrixType& matrix)
    {
      return matrix.getRowStarts().data();
    }

    /**
     * @brief Borrows input storage from a non-const matrix for vendor APIs.
     * The vendor consumes these arrays as inputs; no allocation or copy occurs.
     * @param matrix Non-const owner of the compressed arrays.
     * @return Borrowed pointer, valid until the owner replaces its storage.
     */
    static ::Rodin::Math::SparseIndex* getStarts(MatrixType& matrix)
    {
      return matrix.getRowStarts().data();
    }

    /**
     * @brief Column of every coefficient.
     * @param matrix Compressed device matrix.
     * @return Pointer to the CSR column indices of @p matrix.
     */
    static const ::Rodin::Math::SparseIndex* getIndices(const MatrixType& matrix)
    {
      return matrix.getColumns().data();
    }

    /**
     * @brief Borrows input storage from a non-const matrix for vendor APIs.
     * The vendor consumes these arrays as inputs; no allocation or copy occurs.
     * @param matrix Non-const owner of the compressed arrays.
     * @return Borrowed pointer, valid until the owner replaces its storage.
     */
    static ::Rodin::Math::SparseIndex* getIndices(MatrixType& matrix)
    {
      return matrix.getColumns().data();
    }

    /**
     * @brief Value of every coefficient.
     * @param matrix Compressed device matrix.
     * @return Pointer to the CSR values of @p matrix.
     */
    static const Real* getValues(const MatrixType& matrix)
    {
      return matrix.getValues().data();
    }

    /**
     * @brief Borrows input storage from a non-const matrix for vendor APIs.
     * The vendor consumes these arrays as inputs; no allocation or copy occurs.
     * @param matrix Non-const owner of the compressed arrays.
     * @return Borrowed pointer, valid until the owner replaces its storage.
     */
    static Real* getValues(MatrixType& matrix)
    {
      return matrix.getValues().data();
    }
  };
}

#endif
