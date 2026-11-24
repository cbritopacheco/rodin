/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_MATH_SPATIALMATRIX_H
#define RODIN_KOKKOS_MATH_SPATIALMATRIX_H

/**
 * @file
 * @brief Device-callable bounded spatial matrix.
 */

#include <Kokkos_Core.hpp>
#include "Rodin/Types.h"

namespace Rodin::Kokkos::Math
{
  /**
   * @brief Device-local matrix with at most three rows and columns.
   *
   * Active dimensions include zero for point geometry. Storage never allocates.
   */
  class SpatialMatrix
  {
    public:
      /// @brief Constructs the empty matrix.
      KOKKOS_INLINE_FUNCTION SpatialMatrix()
        : m_rows(0),
          m_cols(0),
          m_data{}
      {}

      /**
       * @brief Constructs a zero matrix.
       * @param rows Active row count.
       * @param cols Active column count.
       */
      KOKKOS_INLINE_FUNCTION SpatialMatrix(Index rows, Index cols)
        : m_rows(rows),
          m_cols(cols),
          m_data{}
      {
        KOKKOS_ASSERT(rows <= 3 && cols <= 3);
      }

      /// @brief Returns the row count. @return Active row count.
      KOKKOS_INLINE_FUNCTION Index rows() const
      {
        return m_rows;
      }

      /// @brief Returns the column count. @return Active column count.
      KOKKOS_INLINE_FUNCTION Index cols() const
      {
        return m_cols;
      }

      /**
       * @brief Accesses a matrix entry.
       * @param row Active row index.
       * @param col Active column index.
       * @return Mutable matrix entry.
       */
      KOKKOS_INLINE_FUNCTION Real& operator()(Index row, Index col)
      {
        KOKKOS_ASSERT(row < m_rows && col < m_cols);
        return m_data[row * 3 + col];
      }

      /**
       * @brief Accesses a matrix entry.
       * @param row Active row index.
       * @param col Active column index.
       * @return SpatialMatrix entry value.
       */
      KOKKOS_INLINE_FUNCTION Real operator()(Index row, Index col) const
      {
        KOKKOS_ASSERT(row < m_rows && col < m_cols);
        return m_data[row * 3 + col];
      }

      /**
       * @brief Computes the determinant of a square matrix.
       * @return Signed determinant, with empty determinant equal to one.
       */
      KOKKOS_INLINE_FUNCTION Real determinant() const
      {
        KOKKOS_ASSERT(m_rows == m_cols);
        if (m_rows == 0)
          return 1;
        if (m_rows == 1)
          return operator()(0, 0);
        if (m_rows == 2)
        {
          return operator()(0, 0) * operator()(1, 1) -
                 operator()(0, 1) * operator()(1, 0);
        }
        return operator()(0, 0) *
          (operator()(1, 1) * operator()(2, 2) - operator()(1, 2) * operator()(2, 1)) -
               operator()(0, 1) *
          (operator()(1, 0) * operator()(2, 2) - operator()(1, 2) * operator()(2, 0)) +
               operator()(0, 2) *
          (operator()(1, 0) * operator()(2, 1) - operator()(1, 1) * operator()(2, 0));
      }

      /**
       * @brief Computes the inverse of a nonsingular square matrix.
       * @return Inverse matrix, including the empty inverse in dimension zero.
       */
      KOKKOS_INLINE_FUNCTION SpatialMatrix inverse() const
      {
        KOKKOS_ASSERT(m_rows == m_cols);
        SpatialMatrix result(m_rows, m_cols);
        const Real det = determinant();
        KOKKOS_ASSERT(det != 0);
        // One division, instead of one per entry.
        const Real inverse = 1 / det;
        if (m_rows == 1)
          result(0, 0) = inverse;
        else if (m_rows == 2)
        {
          result(0, 0) = operator()(1, 1) * inverse;
          result(0, 1) = -operator()(0, 1) * inverse;
          result(1, 0) = -operator()(1, 0) * inverse;
          result(1, 1) = operator()(0, 0) * inverse;
        }
        else if (m_rows == 3)
        {
          for (Index i = 0; i < 3; ++i)
          {
            for (Index j = 0; j < 3; ++j)
            {
              const Index r0 = (j + 1) % 3;
              const Index r1 = (j + 2) % 3;
              const Index c0 = (i + 1) % 3;
              const Index c1 = (i + 2) % 3;
              result(i, j) = (operator()(r0, c0) * operator()(r1, c1) -
                               operator()(r0, c1) * operator()(r1, c0)) *
                inverse;
            }
          }
        }
        return result;
      }

    private:
      Index m_rows;
      Index m_cols;
      ::Kokkos::Array<Real, 9> m_data;
  };
}

#endif
