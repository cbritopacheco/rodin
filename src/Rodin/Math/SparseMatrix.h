/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file SparseMatrix.h
 * @brief Sparse matrix type and operations.
 *
 * This file provides type aliases for sparse matrices and specialized operations
 * for sparse linear algebra. Sparse matrices efficiently store matrices where
 * most elements are zero, which is common in finite element discretizations.
 */
#ifndef RODIN_MATH_SPARSEMATRIX_H
#define RODIN_MATH_SPARSEMATRIX_H

#include <cstdint>

#include <Eigen/Sparse>

#include "Rodin/FormLanguage/Traits.h"

namespace Rodin::Math
{
  /// Signed 32-bit index type used by sparse matrix storage.
  using SparseIndex = std::int32_t;

  /**
   * @brief Sparse matrix entry used during assembly.
   * @tparam ScalarType The entry value type.
   */
  template <class ScalarType>
  using SparseTriplet = Eigen::Triplet<ScalarType, SparseIndex>;

  /**
   * @brief Sparse matrix type.
   *
   * A sparse matrix stores only the non-zero elements, making it memory-efficient
   * for large matrices arising from finite element discretizations where most
   * entries are zero. The matrix uses compressed column storage (CCS) format
   * by default.
   *
   * ## Storage Format
   * Eigen's SparseMatrix uses Compressed Column Storage (CCS), also known as
   * Compressed Sparse Column (CSC) format, which stores:
   * - Non-zero values
   * - Signed 32-bit row indices for each non-zero value
   * - Signed 32-bit column pointers indicating where each column starts
   *
   * ## Typical Applications
   * - Stiffness matrices: @f$ K @f$ in @f$ Ku = f @f$
   * - Mass matrices: @f$ M @f$ in @f$ M\ddot{u} + Ku = f @f$
   * - Discretized differential operators
   *
   * @tparam ScalarType The element type (e.g., Real, Complex)
   */
  template <class ScalarType>
  using SparseMatrix = Eigen::SparseMatrix<ScalarType, Eigen::ColMajor, SparseIndex>;

  /**
   * @brief Performs the AXPY operation on sparse matrices.
   *
   * Computes @f$ y \leftarrow y + \alpha x @f$ where @f$ y @f$ and @f$ x @f$
   * are sparse matrices and @f$ \alpha @f$ is a scalar.
   *
   * This is the sparse matrix version of the Level 1 BLAS AXPY operation.
   *
   * @tparam AScalarType Scalar type for alpha
   * @tparam YScalarType Scalar type for matrix y
   * @tparam XScalarType Scalar type for matrix x
   * @param[in,out] y Matrix to be updated (accumulated result)
   * @param[in] alpha Scalar multiplier for x
   * @param[in] x Matrix to be scaled and added to y
   */
  template <class AScalarType, class YScalarType, class XScalarType>
  void axpy(SparseMatrix<YScalarType>& y, AScalarType alpha, const SparseMatrix<XScalarType>& x)
  {
    y += alpha * x;
  }
}

namespace Rodin::FormLanguage
{
  /// @brief Type traits for a Math::SparseMatrix: exposes the scalar type.
  template <class Number>
  struct Traits<Math::SparseMatrix<Number>>
  {
      /// @brief Scalar value type.
      using ScalarType = Number;
  };
}

#endif
