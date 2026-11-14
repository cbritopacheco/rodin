/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_MATH_SPARSEMATRIX_H
#define RODIN_KOKKOS_MATH_SPARSEMATRIX_H

/**
 * @file
 * @brief Sparse matrix stored in compressed rows in the memory of a Kokkos
 * backend.
 */

#include <string>
#include <utility>
#include <vector>

#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosSparse_spmv.hpp>
#include <KokkosBlas1_scal.hpp>
#include "Vector.h"

#include "Rodin/Types.h"
#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/FormLanguage/Traits.h"
#include "Rodin/Solver/ForwardDecls.h"

#include "Rodin/Kokkos/Context/Kokkos.h"

namespace Rodin::Kokkos::Math
{
  /**
   * @brief Sparse matrix in compressed row storage, resident in the memory
   * space of a Kokkos backend.
   *
   * ## Architecture
   * The owner holds a native KokkosSparse::CrsMatrix and a writable handle to
   * its row offsets for reassembly. Both handles refer to the same allocation.
   * getView() exposes the native matrix without copying it.
   * The device assembly of a Problem builds it without copying the matrix to
   * the host, so that a solver running on the device reads it in place. The
   * columns of a row are sorted and unique. Copies are deep, as for
   * ::Rodin::Math::SparseMatrix; moves keep the storage.
   *
   * @tparam Scalar Coefficient type.
   * @tparam Backend Kokkos backend, as for Context::Kokkos.
   */
  template <class Scalar, class Backend>
  class SparseMatrix
  {
    public:
      /// @brief Coefficient type.
      using ScalarType = Scalar;

      /// @brief Kokkos backend.
      using BackendType = Backend;

      /// @brief Context of the backend.
      using ContextType = ::Rodin::Context::Kokkos<::Rodin::Context::Local, Backend>;

      /// @brief Execution space of the backend.
      using ExecutionSpace = typename ContextType::ExecutionSpaceType;

      /// @brief Memory space holding the matrix.
      using MemorySpace = typename ContextType::MemorySpaceType;

      /// @brief Index type of the compressed arrays.
      using IndexType = ::Rodin::Math::SparseIndex;

      /// @brief Array in the memory space.
      template <class T>
      using Array = ::Kokkos::View<T*, MemorySpace>;

      /// @brief Native Kokkos-Kernels compressed row matrix.
      using View = ::KokkosSparse::CrsMatrix<Scalar, IndexType,
        ::Kokkos::Device<ExecutionSpace, MemorySpace>, void, IndexType>;

      /**
       * @brief Borrows the native compressed row matrix.
       * @return Native matrix sharing this owner's arrays.
       */
      const View& getView() const
      {
        return m_matrix;
      }

      /// @return Number of rows, for the shared algebra API.
      size_t rows() const
      {
        return getRowCount();
      }

      /// @return Number of columns, for the shared algebra API.
      size_t cols() const
      {
        return getColumnCount();
      }

      /// @return Number of stored coefficients.
      size_t nonZeros() const
      {
        return getNonZeroCount();
      }

      /// @return Frobenius norm of the stored coefficients.
      auto norm() const
      {
        return ::KokkosBlas::nrm2(m_matrix.values);
      }

      /// @param rhs Scalar multiplier. @return This matrix.
      SparseMatrix& operator*=(Scalar rhs)
      {
        ::KokkosBlas::scal(m_matrix.values, rhs, m_matrix.values);
        return *this;
      }

      /// @param rhs Scalar multiplier. @return Independent scaled matrix.
      SparseMatrix operator*(Scalar rhs) const
      {
        SparseMatrix result(*this);
        result *= rhs;
        return result;
      }

      /// @param rhs Backend vector. @return Backend matrix-vector product.
      Vector<Scalar, Backend> operator*(const Vector<Scalar, Backend>& rhs) const
      {
        assert(rhs.size() == getColumnCount());
        Vector<Scalar, Backend> result(getRowCount());
        ::KokkosSparse::spmv(
          "N", Scalar(1), m_matrix, rhs.getView(), Scalar(0), result.getView());
        return result;
      }

      /// @brief Empty matrix.
      SparseMatrix() = default;

      /**
       * @brief Matrix holding compressed arrays.
       * @param rows Number of rows.
       * @param cols Number of columns.
       * @param starts Position of the first coefficient of every row,
       * followed by the number of coefficients.
       * @param columns Column of every coefficient.
       * @param values Value of every coefficient.
       */
      SparseMatrix(size_t rows, size_t cols, Array<IndexType> starts,
        Array<IndexType> columns, Array<Scalar> values)
        : m_rows(rows),
          m_cols(cols),
          m_starts(starts),
          m_matrix("Rodin matrix", rows, cols, values.extent(0), values, starts, columns)
      {}

      /**
       * @brief Copies the coefficients.
       * @param other Matrix to copy.
       */
      SparseMatrix(const SparseMatrix& other)
        : SparseMatrix(other.m_rows, other.m_cols, clone(other.getRowStarts()),
            clone(other.getColumns()), clone(other.getValues()))
      {}

      /**
       * @brief Moves the storage.
       * @param other Matrix to move.
       */
      SparseMatrix(SparseMatrix&& other) noexcept
        : m_rows(std::exchange(other.m_rows, 0)),
          m_cols(std::exchange(other.m_cols, 0)),
          m_starts(std::move(other.m_starts)),
          m_matrix(std::move(other.m_matrix))
      {}

      /**
       * @brief Copies the coefficients.
       * @param other Matrix to copy.
       * @returns Reference to this matrix.
       */
      SparseMatrix& operator=(const SparseMatrix& other)
      {
        if (this != &other)
          *this = SparseMatrix(other);
        return *this;
      }

      /**
       * @brief Moves the storage.
       * @param other Matrix to move.
       * @returns Reference to this matrix.
       */
      SparseMatrix& operator=(SparseMatrix&& other) noexcept
      {
        m_rows = std::exchange(other.m_rows, 0);
        m_cols = std::exchange(other.m_cols, 0);
        m_starts = std::move(other.m_starts);
        m_matrix = std::move(other.m_matrix);
        return *this;
      }

      /**
       * @brief Number of rows.
       * @return Row count of the matrix.
       */
      size_t getRowCount() const
      {
        return m_rows;
      }

      /**
       * @brief Number of columns.
       * @return Column count of the matrix.
       */
      size_t getColumnCount() const
      {
        return m_cols;
      }

      /**
       * @brief Number of stored coefficients.
       * @return Length of the column and value arrays.
       */
      size_t getNonZeroCount() const
      {
        return m_matrix.values.extent(0);
      }

      /**
       * @brief Position of the first coefficient of every row, followed by
       * the number of coefficients. Empty for an empty matrix.
       * @return Resident CSR row offsets.
       */
      const Array<IndexType>& getRowStarts() const
      {
        return m_starts;
      }

      /**
       * @brief Column of every coefficient.
       * @return Resident CSR column indices.
       */
      const typename View::index_type& getColumns() const
      {
        return m_matrix.graph.entries;
      }

      /**
       * @brief Value of every coefficient.
       * @return Resident CSR values, read-only.
       */
      const typename View::values_type& getValues() const
      {
        return m_matrix.values;
      }

      /**
       * @brief Value of every coefficient, for an update in place.
       * @return Resident CSR values, writable.
       */
      typename View::values_type& getValues()
      {
        return m_matrix.values;
      }

      /**
       * @brief Copies the matrix to the host.
       * @returns Matrix in compressed column storage.
       */
      ::Rodin::Math::SparseMatrix<Scalar> toHost() const
      {
        ::Rodin::Math::SparseMatrix<Scalar> res(m_rows, m_cols);
        if (m_matrix.graph.row_map.extent(0) == 0)
          return res;
        const auto starts = ::Kokkos::create_mirror_view_and_copy(
          ::Kokkos::HostSpace(), m_matrix.graph.row_map);
        const auto columns = ::Kokkos::create_mirror_view_and_copy(
          ::Kokkos::HostSpace(), m_matrix.graph.entries);
        const auto values =
          ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace(), m_matrix.values);
        std::vector<::Rodin::Math::SparseTriplet<Scalar>> triplets;
        triplets.reserve(getNonZeroCount());
        for (size_t i = 0; i < m_rows; ++i)
        {
          for (IndexType k = starts(i); k < starts(i + 1); ++k)
            triplets.emplace_back(static_cast<IndexType>(i), columns(k), values(k));
        }
        res.setFromTriplets(triplets.begin(), triplets.end());
        return res;
      }

    private:
      /**
       * @brief Deep copy of a device array, with the same label.
       * @param a Array to copy.
       * @return New array holding a copy of @p a, empty when @p a is empty.
       */
      template <class V>
      static Array<typename V::non_const_value_type> clone(const V& a)
      {
        using T = typename V::non_const_value_type;
        if (a.extent(0) == 0)
          return Array<T>();
        Array<T> res(
          ::Kokkos::view_alloc(::Kokkos::WithoutInitializing, a.label()), a.extent(0));
        ::Kokkos::deep_copy(res, a);
        return res;
      }

      size_t m_rows = 0;
      size_t m_cols = 0;
      Array<IndexType> m_starts;
      View m_matrix;
  };
}

namespace Rodin::FormLanguage
{
  /// @brief Type traits for a Kokkos::Math::SparseMatrix, exposing the scalar type.
  template <class Scalar, class Backend>
  struct Traits<::Rodin::Kokkos::Math::SparseMatrix<Scalar, Backend>>
  {
    /// @brief Scalar value type.
    using ScalarType = Scalar;
  };
}

namespace Rodin::Solver::Internal
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
     * @brief Column of every coefficient.
     * @param matrix Compressed device matrix.
     * @return Pointer to the CSR column indices of @p matrix.
     */
    static const ::Rodin::Math::SparseIndex* getIndices(const MatrixType& matrix)
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
  };
}

#endif
