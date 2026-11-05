/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EIKONAL_METRIC_H
#define RODIN_EIKONAL_METRIC_H

/**
 * @file
 * @brief Travel metric induced by a symmetric positive definite velocity tensor.
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <Eigen/Cholesky>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Math.h"

namespace Rodin::Eikonal
{
  /**
   * @brief Bilinear travel metric induced by a squared-velocity tensor.
   *
   * Construction validates an ambient SPD squared-velocity tensor @f$M@f$.
   * The stored Cholesky factor evaluates @f$a^T M^{-1}b@f$ without explicitly
   * forming the inverse. Vectors have the physical embedding dimension;
   * embedded-cell solvers restrict this ambient metric to physical tangents.
   */
  class Metric final
  {
    using Storage = Math::SpatialMatrix<Real>::Data;
    using Vector = Math::FixedSizeVector<Real, Math::SpatialMatrix<Real>::MaxSize>;

    /// @brief Relative symmetry check allowing a small number of rounding errors.
    static constexpr Real SymmetryRoundoffFactor = Real(64);

    public:
      /**
       * @brief Validates and factors a squared-velocity tensor.
       * @tparam Matrix Matrix value type with scalar entry access.
       * @param velocity Ambient squared-velocity tensor, not its inverse.
       * @param dimension Physical embedding dimension.
       */
      template <class Matrix>
      Metric(const Matrix& velocity, size_t dimension)
        : m_dimension(dimension)
      {
        if (dimension == 0 || dimension > Math::SpatialMatrix<Real>::MaxSize ||
          velocity.rows() != dimension || velocity.cols() != dimension)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Velocity tensor must match the physical embedding dimension."
            << Alert::Raise;
        }
        Storage matrix = Storage::Identity();
        Real scale = 0;
        for (size_t row = 0; row < dimension; ++row)
        {
          for (size_t column = 0; column < dimension; ++column)
          {
            const Real value = velocity(row, column);
            if (!std::isfinite(value))
            {
              Alert::MemberFunctionException(*this, __func__)
                << "Velocity tensor entries must be finite." << Alert::Raise;
            }
            matrix(row, column) = value;
            scale = std::max(scale, std::abs(value));
          }
        }
        const Real tolerance =
          SymmetryRoundoffFactor * std::numeric_limits<Real>::epsilon() * scale;
        for (size_t row = 0; row < dimension; ++row)
        {
          for (size_t column = 0; column < row; ++column)
          {
            if (std::abs(matrix(row, column) - matrix(column, row)) > tolerance)
            {
              Alert::MemberFunctionException(*this, __func__)
                << "Velocity tensor must be symmetric." << Alert::Raise;
            }
            matrix(row, column) = matrix(column, row) =
              (matrix(row, column) + matrix(column, row)) / 2;
          }
        }
        m_factor.compute(matrix);
        if (m_factor.info() != Eigen::Success)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Velocity tensor must be positive definite." << Alert::Raise;
        }
      }

      /**
       * @brief Contracts two physical displacements in the inverse velocity metric.
       * @param lhs First physical displacement.
       * @param rhs Second physical displacement.
       * @returns The real metric product @f$lhs^T M^{-1}rhs@f$.
       */
      Real dot(
        const Math::SpatialVector<Real>& lhs, const Math::SpatialVector<Real>& rhs) const
      {
        assert(lhs.size() == m_dimension && rhs.size() == m_dimension);
        Vector left, right;
        left.setZero();
        right.setZero();
        for (size_t coordinate = 0; coordinate < m_dimension; ++coordinate)
        {
          left(coordinate) = lhs(coordinate);
          right(coordinate) = rhs(coordinate);
        }
        const Vector inverseRight = m_factor.solve(right);
        return Math::dot(left, inverseRight);
      }

    private:
      size_t m_dimension;
      Eigen::LLT<Storage> m_factor;
  };
}

#endif
