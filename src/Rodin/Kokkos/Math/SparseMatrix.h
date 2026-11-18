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
       * @brief Borrows writable CSR storage without reallocating it.
       * The caller must preserve the CSR invariants and explicitly refactorize
       * any solver after changing the matrix.
       * @return Resident array handle; valid until storage replacement.
       */
      Array<IndexType>& getRowStarts()
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
       * @brief Borrows writable CSR storage without reallocating it.
       * The caller must preserve the CSR invariants and explicitly refactorize
       * any solver after changing the matrix.
       * @return Resident array handle; valid until storage replacement.
       */
      typename View::index_type& getColumns()
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
       * @brief Replaces this backend matrix with an independent host snapshot.
       * @param source Host matrix, whose column storage is converted to rows.
       * @returns This matrix, after every transfer completes.
       */
      SparseMatrix& upload(const ::Rodin::Math::SparseMatrix<Scalar>& source)
      {
        Eigen::SparseMatrix<Scalar, Eigen::RowMajor, IndexType> rows(source);
        rows.makeCompressed();
        const ExecutionSpace execution;
        Array<IndexType> starts(::Kokkos::view_alloc(execution,
                                  ::Kokkos::WithoutInitializing, "Rodin uploaded rows"),
          rows.rows() + 1);
        Array<IndexType> columns(
          ::Kokkos::view_alloc(
            execution, ::Kokkos::WithoutInitializing, "Rodin uploaded columns"),
          rows.nonZeros());
        Array<Scalar> values(::Kokkos::view_alloc(execution,
                               ::Kokkos::WithoutInitializing, "Rodin uploaded values"),
          rows.nonZeros());
        using HostIndices = ::Kokkos::View<const IndexType*, ::Kokkos::HostSpace,
          ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>>;
        using HostValues = ::Kokkos::View<const Scalar*, ::Kokkos::HostSpace,
          ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>>;
        ::Kokkos::deep_copy(
          execution, starts, HostIndices(rows.outerIndexPtr(), rows.rows() + 1));
        ::Kokkos::deep_copy(
          execution, columns, HostIndices(rows.innerIndexPtr(), rows.nonZeros()));
        ::Kokkos::deep_copy(
          execution, values, HostValues(rows.valuePtr(), rows.nonZeros()));
        execution.fence("Rodin matrix upload");
        *this = SparseMatrix(rows.rows(), rows.cols(), starts, columns, values);
        return *this;
      }

      /**
       * @brief Replaces a host matrix with an independent backend snapshot.
       * @param res Host destination, rebuilt in compressed column storage.
       * The backend is unchanged and transfers complete before returning.
       */
      void download(::Rodin::Math::SparseMatrix<Scalar>& res) const
      {
        res.resize(m_rows, m_cols);
        res.setZero();
        if (m_matrix.graph.row_map.extent(0) == 0)
          return;
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
      }

      /// @returns Independent host snapshot, in compressed column storage.
      ::Rodin::Math::SparseMatrix<Scalar> download() const
      {
        ::Rodin::Math::SparseMatrix<Scalar> res;
        download(res);
        return res;
      }

      /// @returns Independent host snapshot; compatibility spelling of download().
      ::Rodin::Math::SparseMatrix<Scalar> toHost() const
      {
        return download();
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

#endif
