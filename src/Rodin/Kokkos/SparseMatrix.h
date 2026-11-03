/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_SPARSEMATRIX_H
#define RODIN_KOKKOS_SPARSEMATRIX_H

/**
 * @file
 * @brief Sparse matrix stored in compressed rows in the memory of a Kokkos
 * backend.
 */

#include <string>
#include <utility>
#include <vector>

#include <Kokkos_Core.hpp>

#include "Rodin/Types.h"
#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/FormLanguage/Traits.h"
#include "Rodin/Solver/ForwardDecls.h"

#include "Context/Kokkos.h"

namespace Rodin::Kokkos
{
  /**
   * @brief Sparse matrix in compressed row storage, resident in the memory
   * space of a Kokkos backend.
   *
   * The device assembly of a Problem builds it without copying the matrix to
   * the host, so that a solver running on the device reads it in place. The
   * columns of a row are sorted and unique. Copies are deep, as for
   * Math::SparseMatrix; moves keep the storage.
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
      using IndexType = Math::SparseIndex;

      /// @brief Array in the memory space.
      template <class T>
      using Array = ::Kokkos::View<T*, MemorySpace>;

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
          m_starts(std::move(starts)),
          m_columns(std::move(columns)),
          m_values(std::move(values))
      {}

      /**
       * @brief Copies the coefficients.
       * @param other Matrix to copy.
       */
      SparseMatrix(const SparseMatrix& other)
        : m_rows(other.m_rows),
          m_cols(other.m_cols),
          m_starts(clone(other.m_starts)),
          m_columns(clone(other.m_columns)),
          m_values(clone(other.m_values))
      {}

      /**
       * @brief Moves the storage.
       * @param other Matrix to move.
       */
      SparseMatrix(SparseMatrix&& other) noexcept
        : m_rows(std::exchange(other.m_rows, 0)),
          m_cols(std::exchange(other.m_cols, 0)),
          m_starts(std::move(other.m_starts)),
          m_columns(std::move(other.m_columns)),
          m_values(std::move(other.m_values))
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
        m_columns = std::move(other.m_columns);
        m_values = std::move(other.m_values);
        return *this;
      }

      /// @brief Number of rows.
      size_t getRowCount() const
      {
        return m_rows;
      }

      /// @brief Number of columns.
      size_t getColumnCount() const
      {
        return m_cols;
      }

      /// @brief Number of stored coefficients.
      size_t getNonZeroCount() const
      {
        return m_values.extent(0);
      }

      /// @brief Position of the first coefficient of every row, followed by
      /// the number of coefficients. Empty for an empty matrix.
      const Array<IndexType>& getRowStarts() const
      {
        return m_starts;
      }

      /// @brief Column of every coefficient.
      const Array<IndexType>& getColumns() const
      {
        return m_columns;
      }

      /// @brief Value of every coefficient.
      const Array<Scalar>& getValues() const
      {
        return m_values;
      }

      /// @brief Value of every coefficient, for an update in place.
      Array<Scalar>& getValues()
      {
        return m_values;
      }

      /**
       * @brief Copies the matrix to the host.
       * @returns Matrix in compressed column storage.
       */
      Math::SparseMatrix<Scalar> toHost() const
      {
        Math::SparseMatrix<Scalar> res(m_rows, m_cols);
        if (m_starts.extent(0) == 0)
          return res;
        const auto starts =
          ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace(), m_starts);
        const auto columns =
          ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace(), m_columns);
        const auto values =
          ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace(), m_values);
        std::vector<Math::SparseTriplet<Scalar>> triplets;
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
      template <class T>
      static Array<T> clone(const Array<T>& a)
      {
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
      Array<IndexType> m_columns;
      Array<Scalar> m_values;
  };
}

namespace Rodin::FormLanguage
{
  /// @brief Type traits for a Kokkos::SparseMatrix, exposing the scalar type.
  template <class Scalar, class Backend>
  struct Traits<::Rodin::Kokkos::SparseMatrix<Scalar, Backend>>
  {
    /// @brief Scalar value type.
    using ScalarType = Scalar;
  };
}

namespace Rodin::Solver::Internal
{
  /**
   * @brief CSR arrays of a Kokkos::SparseMatrix, read in place by the GPU
   * direct solvers when they are in device memory and copied to the device
   * otherwise.
   */
  template <class Backend>
  struct CuDSSOperator<::Rodin::Kokkos::SparseMatrix<Real, Backend>>
  {
    /// @brief Matrix type.
    using MatrixType = ::Rodin::Kokkos::SparseMatrix<Real, Backend>;

    /// @brief Whether the arrays are those of the transpose.
    static constexpr bool Transposed = false;

    /// @brief Whether the arrays are in device memory.
    static constexpr bool Resident = !::Kokkos::SpaceAccessibility<::Kokkos::HostSpace,
      typename MatrixType::MemorySpace>::accessible;

    /// @brief Leaves the matrix unchanged, being already compressed.
    static void prepare(MatrixType&) {}

    /// @brief Number of rows.
    static size_t getRowCount(const MatrixType& matrix)
    {
      return matrix.getRowCount();
    }

    /// @brief Number of columns.
    static size_t getColumnCount(const MatrixType& matrix)
    {
      return matrix.getColumnCount();
    }

    /// @brief Number of stored coefficients.
    static size_t getNonZeroCount(const MatrixType& matrix)
    {
      return matrix.getNonZeroCount();
    }

    /// @brief Start of every row, followed by the coefficient count.
    static const Math::SparseIndex* getStarts(const MatrixType& matrix)
    {
      return matrix.getRowStarts().data();
    }

    /// @brief Column of every coefficient.
    static const Math::SparseIndex* getIndices(const MatrixType& matrix)
    {
      return matrix.getColumns().data();
    }

    /// @brief Value of every coefficient.
    static const Real* getValues(const MatrixType& matrix)
    {
      return matrix.getValues().data();
    }
  };
}

#endif
