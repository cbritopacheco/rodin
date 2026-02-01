/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Matrix.h
 * @brief Dense matrix type aliases and definitions.
 *
 * This file provides type aliases for various dense matrix types built on Eigen.
 * Matrices are stored in column-major order by default and support both dynamic
 * and fixed sizes.
 */
#ifndef RODIN_MATH_MATRIX_H
#define RODIN_MATH_MATRIX_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/src/Core/ArrayBase.h>

#include "Rodin/Configure.h"

#include "Rodin/FormLanguage/Traits.h"

#include "Rodin/Types.h"

#include "ForwardDecls.h"

namespace Rodin::Math
{
  /**
   * @brief Dynamic-size dense matrix type.
   *
   * A dense matrix with both dimensions determined at runtime. The matrix
   * is stored in column-major order (Fortran style) and supports standard
   * linear algebra operations.
   *
   * @tparam ScalarType The element type (e.g., Real, Complex)
   */
  template <class ScalarType>
  using Matrix = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;

  /**
   * @brief Real-valued dense matrix.
   *
   * Convenience alias for Matrix<Real>.
   */
  using RealMatrix = Matrix<Real>;

  /**
   * @brief Complex-valued dense matrix.
   *
   * Convenience alias for Matrix<Complex>, commonly used in eigenvalue
   * problems, frequency-domain analysis, and other applications requiring
   * complex arithmetic.
   */
  using ComplexMatrix = Matrix<Complex>;

  /**
   * @brief Point matrix with bounded row dimension.
   *
   * A dynamic-size matrix where the number of rows is bounded by
   * RODIN_MAXIMAL_SPACE_DIMENSION. Commonly used to store collections of
   * spatial points as columns.
   *
   * Example: storing vertex coordinates of a mesh element where each column
   * represents one vertex's coordinates.
   */
  using PointMatrix =
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, 0, RODIN_MAXIMAL_SPACE_DIMENSION, Eigen::Dynamic>;

  /**
   * @brief Fixed-size matrix type.
   *
   * A compile-time fixed-size matrix. Both dimensions are known at compile time,
   * enabling optimizations and stack allocation.
   *
   * Example uses:
   * - 2×2 rotation matrices
   * - 3×3 transformation matrices
   * - Small dense blocks in larger systems
   *
   * @tparam ScalarType The element type
   * @tparam Rows The number of rows (must be known at compile time)
   * @tparam Cols The number of columns (must be known at compile time)
   */
  template <class ScalarType, size_t Rows, size_t Cols>
  using FixedSizeMatrix = Eigen::Matrix<ScalarType, Rows, Cols>;

  /**
   * @brief Spatial matrix with bounded maximum dimensions.
   *
   * A dynamic-size matrix with maximum dimensions bounded by
   * RODIN_MAXIMAL_SPACE_DIMENSION. Used for geometric transformations,
   * Jacobians, and other spatial operators to optimize memory allocation.
   *
   * Typical uses:
   * - Jacobian matrices: @f$ J = \frac{\partial \mathbf{f}}{\partial \mathbf{x}} @f$
   * - Metric tensors
   * - Local coordinate transformations
   *
   * @tparam ScalarType The element type
   */
  template <class ScalarType>
  class SpatialMatrix
  {
    public:
      using Scalar = ScalarType;
      using Data = Eigen::Matrix3<ScalarType>;
      static constexpr std::uint8_t MaxSize = RODIN_MAXIMAL_SPACE_DIMENSION;
      static_assert(MaxSize == 3, "MaxSize must be equal to 3.");

      constexpr
      SpatialMatrix() noexcept
        : m_rows(0), m_cols(0)
      {}

      constexpr
      SpatialMatrix(std::uint8_t rows, std::uint8_t cols)
        : m_rows(rows), m_cols(cols)
      {
        assert(rows <= MaxSize && cols <= MaxSize);
      }

      template <class EigenDerived>
      constexpr
      SpatialMatrix(const Eigen::MatrixBase<EigenDerived>& other)
        : m_rows(static_cast<std::uint8_t>(other.rows())),
          m_cols(static_cast<std::uint8_t>(other.cols())),
          m_data(other)
      {
        assert(m_rows <= MaxSize && m_cols <= MaxSize);
      }

      template <class EigenDerived>
      constexpr
      SpatialMatrix& operator=(const Eigen::ArrayBase<EigenDerived>& other)
      {
        const std::uint8_t r = static_cast<std::uint8_t>(other.rows());
        const std::uint8_t c = static_cast<std::uint8_t>(other.cols());
        assert(r <= MaxSize && c <= MaxSize);
        m_rows = r;
        m_cols = c;
        this->eigen() = other;
        return *this;
      }

      template <class EigenDerived>
      constexpr
      SpatialMatrix& operator=(const Eigen::MatrixBase<EigenDerived>& other)
      {
        const std::uint8_t r = static_cast<std::uint8_t>(other.rows());
        const std::uint8_t c = static_cast<std::uint8_t>(other.cols());
        assert(r <= MaxSize && c <= MaxSize);
        m_rows = r;
        m_cols = c;
        this->eigen() = other;
        return *this;
      }

      constexpr
      std::uint8_t rows() const noexcept
      {
        return m_rows;
      }

      constexpr
      std::uint8_t cols() const noexcept
      {
        return m_cols;
      }

      constexpr
      void resize(std::uint8_t r, std::uint8_t c)
      {
        assert(r <= MaxSize && c <= MaxSize);
        m_rows = r;
        m_cols = c;
      }

      constexpr
      ScalarType& operator()(std::uint8_t i, std::uint8_t j)
      {
        assert(i < m_rows && j < m_cols);
        return m_data(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
      }

      constexpr
      const ScalarType& operator()(std::uint8_t i, std::uint8_t j) const
      {
        assert(i < m_rows && j < m_cols);
        return m_data(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
      }

      constexpr
      void setZero() noexcept
      {
        m_data.setZero();
      }

      constexpr
      void setConstant(const ScalarType& value) noexcept
      {
        m_data.setConstant(value);
      }

      constexpr
      auto dot(const SpatialMatrix& other) const noexcept
      {
        assert(m_cols == other.m_rows);
        assert(m_rows == other.m_rows);
        return this->eigen().dot(other.eigen());
      }

      template <class EigenDerived>
      constexpr
      auto dot(const Eigen::MatrixBase<EigenDerived>& other) const noexcept
      {
        return this->eigen().dot(other);
      }

      constexpr
      auto transpose() const noexcept
      {
        return this->eigen().transpose();
      }

      constexpr
      auto value() const noexcept
      {
        return this->eigen().value();
      }

      constexpr
      auto squaredNorm() const noexcept
      {
        return this->eigen().squaredNorm();
      }

      constexpr
      auto inverse() const noexcept
      {
        return this->eigen().inverse();
      }

      constexpr
      auto determinant() const noexcept
      {
        return this->eigen().determinant();
      }

      constexpr
      auto row(std::uint8_t i) noexcept
      {
        assert(i < m_rows);
        return this->eigen().row(static_cast<Eigen::Index>(i));
      }

      constexpr
      auto col(std::uint8_t j) noexcept
      {
        assert(j < m_cols);
        return this->eigen().col(static_cast<Eigen::Index>(j));
      }

      constexpr
      auto eigen() noexcept
      {
        return m_data.topLeftCorner(static_cast<Eigen::Index>(m_rows), static_cast<Eigen::Index>(m_cols));
      }

      constexpr
      auto eigen() const noexcept
      {
        return m_data.topLeftCorner(static_cast<Eigen::Index>(m_rows), static_cast<Eigen::Index>(m_cols));
      }

    private:
      std::uint8_t m_rows;
      std::uint8_t m_cols;
      Data m_data;
  };

  template <class LHS, class Scalar>
  auto operator*(
    const LHS& s,
    const SpatialMatrix<Scalar>& m)
  {
    return s * m.eigen();
  }

  template <class Scalar, class RHS>
  auto operator*(
    const SpatialMatrix<Scalar>& m,
    const RHS& s)
  {
    return m.eigen() * s;
  }

  template <class Scalar>
  auto operator*(const SpatialMatrix<Scalar>& a, const SpatialMatrix<Scalar>& b)
  {
    assert(a.cols() == b.rows());
    return a.eigen() * b.eigen();
  }

  template <class Scalar>
  auto operator*(const SpatialMatrix<Scalar>& m, const SpatialVector<Scalar>& v)
  {
    assert(m.cols() == v.size());
    return m.eigen() * v.eigen();
  }

  template <class Scalar>
  auto operator*(const SpatialVector<Scalar>& v, const SpatialMatrix<Scalar>& m)
  {
    assert(v.size() == m.rows());
    return v.eigen() * m.eigen();
  }
}

namespace Rodin::FormLanguage
{
  template <class Number>
  struct Traits<Math::Matrix<Number>>
  {
    using ScalarType = Number;
  };

  template <class Number>
  struct Traits<Math::SpatialMatrix<Number>>
  {
    using ScalarType = Number;
  };
}

#endif

