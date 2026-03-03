#ifndef RODIN_MATH_SPATIALMATRIX_H
#define RODIN_MATH_SPATIALMATRIX_H

#include "Rodin/FormLanguage/Traits.h"

#include "ForwardDecls.h"
#include "Vector.h"
#include "Common.h"
#include "Traits.h"

namespace Rodin::Math
{
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

      constexpr
      SpatialMatrix(const SpatialVector<Scalar>& vec)
        : m_rows(static_cast<std::uint8_t>(vec.size())),
          m_cols(1)
      {
        switch (m_rows)
        {
          case 3:
            m_data(2,0) = vec(2);
            [[fallthrough]];
          case 2:
            m_data(1,0) = vec(1);
            [[fallthrough]];
          case 1:
            m_data(0,0) = vec(0);
            break;
          case 0:
            break;
          default:
            assert(false);
        }
      }

      constexpr
      SpatialMatrix(const SpatialMatrix&) = default;

      constexpr
      SpatialMatrix(SpatialMatrix&&) = default;

      constexpr
      SpatialMatrix& operator=(const SpatialMatrix&) = default;

      constexpr
      SpatialMatrix& operator=(SpatialMatrix&&) = default;

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
        switch (m_rows)
        {
          case 3:
            switch (m_cols)
            {
              case 3:
                m_data(2,2) = other(static_cast<Eigen::Index>(2), static_cast<Eigen::Index>(2));
                [[fallthrough]];
              case 2:
                m_data(2,1) = other(static_cast<Eigen::Index>(2), static_cast<Eigen::Index>(1));
                [[fallthrough]];
              case 1:
                m_data(2,0) = other(static_cast<Eigen::Index>(2), static_cast<Eigen::Index>(0));
                break;
              case 0:
                break;
              default:
                assert(false);
            }
            [[fallthrough]];
          case 2:
            switch (m_cols)
            {
              case 3:
                m_data(1,2) = other(static_cast<Eigen::Index>(1), static_cast<Eigen::Index>(2));
                [[fallthrough]];
              case 2:
                m_data(1,1) = other(static_cast<Eigen::Index>(1), static_cast<Eigen::Index>(1));
                [[fallthrough]];
              case 1:
                m_data(1,0) = other(static_cast<Eigen::Index>(1), static_cast<Eigen::Index>(0));
                break;
              case 0:
                break;
              default:
                assert(false);
            }
            [[fallthrough]];
          case 1:
            switch (m_cols)
            {
              case 3:
                m_data(0,2) = other(static_cast<Eigen::Index>(0), static_cast<Eigen::Index>(2));
                [[fallthrough]];
              case 2:
                m_data(0,1) = other(static_cast<Eigen::Index>(0), static_cast<Eigen::Index>(1));
                [[fallthrough]];
              case 1:
                m_data(0,0) = other(static_cast<Eigen::Index>(0), static_cast<Eigen::Index>(0));
                break;
              case 0:
                break;
              default:
                assert(false);
            }
            break;
          case 0:
            break;
          default:
            assert(false);
        }
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
        switch (m_rows)
        {
          case 3:
            switch (m_cols)
            {
              case 3:
                m_data(2,2) = other(static_cast<Eigen::Index>(2), static_cast<Eigen::Index>(2));
                [[fallthrough]];
              case 2:
                m_data(2,1) = other(static_cast<Eigen::Index>(2), static_cast<Eigen::Index>(1));
                [[fallthrough]];
              case 1:
                m_data(2,0) = other(static_cast<Eigen::Index>(2), static_cast<Eigen::Index>(0));
                break;
              case 0:
                break;
              default:
                assert(false);
            }
            [[fallthrough]];
          case 2:
            switch (m_cols)
            {
              case 3:
                m_data(1,2) = other(static_cast<Eigen::Index>(1), static_cast<Eigen::Index>(2));
                [[fallthrough]];
              case 2:
                m_data(1,1) = other(static_cast<Eigen::Index>(1), static_cast<Eigen::Index>(1));
                [[fallthrough]];
              case 1:
                m_data(1,0) = other(static_cast<Eigen::Index>(1), static_cast<Eigen::Index>(0));
                break;
              case 0:
                break;
              default:
                assert(false);
            }
            [[fallthrough]];
          case 1:
            switch (m_cols)
            {
              case 3:
                m_data(0,2) = other(static_cast<Eigen::Index>(0), static_cast<Eigen::Index>(2));
                [[fallthrough]];
              case 2:
                m_data(0,1) = other(static_cast<Eigen::Index>(0), static_cast<Eigen::Index>(1));
                [[fallthrough]];
              case 1:
                m_data(0,0) = other(static_cast<Eigen::Index>(0), static_cast<Eigen::Index>(0));
                break;
              case 0:
                break;
              default:
                assert(false);
            }
            break;
          case 0:
            break;
          default:
            assert(false);
        }
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
      void setIdentity() noexcept
      {
        m_data.setIdentity();
      }

      [[nodiscard]] constexpr
      ScalarType squaredNorm() const noexcept
      {
        const auto r = m_rows;
        const auto c = m_cols;
        const auto& A = m_data;

        // key in [0..15] for r,c in [0..3]
        switch (static_cast<unsigned>(r) * 4u + static_cast<unsigned>(c))
        {
          case 0u:  // 0x0
          case 1u:  // 0x1
          case 2u:  // 0x2
          case 3u:  // 0x3
          case 4u:  // 1x0
          case 8u:  // 2x0
          case 12u: // 3x0
            return ScalarType(0);

          case 5u: // 1x1
            return Math::pow2(A(0,0));

          case 6u: // 1x2
            return Math::pow2(A(0,0)) + Math::pow2(A(0,1));

          case 7u: // 1x3
            return Math::pow2(A(0,0)) + Math::pow2(A(0,1)) + Math::pow2(A(0,2));

          case 9u: // 2x1
            return Math::pow2(A(0,0)) + Math::pow2(A(1,0));

          case 10u: // 2x2
            return Math::pow2(A(0,0)) + Math::pow2(A(0,1))
                 + Math::pow2(A(1,0)) + Math::pow2(A(1,1));

          case 11u: // 2x3
            return Math::pow2(A(0,0)) + Math::pow2(A(0,1)) + Math::pow2(A(0,2))
                 + Math::pow2(A(1,0)) + Math::pow2(A(1,1)) + Math::pow2(A(1,2));

          case 13u: // 3x1
            return Math::pow2(A(0,0)) + Math::pow2(A(1,0)) + Math::pow2(A(2,0));

          case 14u: // 3x2
            return Math::pow2(A(0,0)) + Math::pow2(A(0,1))
                 + Math::pow2(A(1,0)) + Math::pow2(A(1,1))
                 + Math::pow2(A(2,0)) + Math::pow2(A(2,1));

          case 15u: // 3x3
            return Math::pow2(A(0,0)) + Math::pow2(A(0,1)) + Math::pow2(A(0,2))
                 + Math::pow2(A(1,0)) + Math::pow2(A(1,1)) + Math::pow2(A(1,2))
                 + Math::pow2(A(2,0)) + Math::pow2(A(2,1)) + Math::pow2(A(2,2));

          default:
            assert(false);
            return ScalarType(0);
        }
      }

      [[nodiscard]] constexpr
      ScalarType norm() const noexcept
      {
        return Math::sqrt(this->squaredNorm());
      }

      constexpr
      ScalarType dot(const SpatialMatrix& other) const noexcept
      {
        assert(m_rows == other.m_rows);
        assert(m_cols == other.m_cols);

        const std::uint8_t r = m_rows;
        const std::uint8_t c = m_cols;

        // Access m_data directly (avoid topLeftCorner expressions).
        const auto& A = m_data;
        const auto& B = other.m_data;

        // Flatten (r,c) into one switch value in [0..9]
        switch (static_cast<unsigned>(r) * 4u + static_cast<unsigned>(c))
        {
          case 0u: // 0x0
            return ScalarType(0);

          case 1u: // 0x1
            return ScalarType(0);

          case 2u: // 0x2
            return ScalarType(0);

          case 3u: // 0x3
            return ScalarType(0);

          case 4u: // 1x0
            return ScalarType(0);

          case 5u: // 1x1
            return A(0,0) * B(0,0);

          case 6u: // 1x2
            return A(0,0) * B(0,0)
                 + A(0,1) * B(0,1);

          case 7u: // 1x3
            return A(0,0) * B(0,0)
                 + A(0,1) * B(0,1)
                 + A(0,2) * B(0,2);

          case 8u: // 2x0
            return ScalarType(0);

          case 9u: // 2x1
            return A(0,0) * B(0,0)
                 + A(1,0) * B(1,0);

          case 10u: // 2x2
            return A(0,0) * B(0,0)
                 + A(0,1) * B(0,1)
                 + A(1,0) * B(1,0)
                 + A(1,1) * B(1,1);

          case 11u: // 2x3
            return A(0,0) * B(0,0)
                 + A(0,1) * B(0,1)
                 + A(0,2) * B(0,2)
                 + A(1,0) * B(1,0)
                 + A(1,1) * B(1,1)
                 + A(1,2) * B(1,2);

          case 12u: // 3x0
            return ScalarType(0);

          case 13u: // 3x1
            return A(0,0) * B(0,0)
                 + A(1,0) * B(1,0)
                 + A(2,0) * B(2,0);

          case 14u: // 3x2
            return A(0,0) * B(0,0)
                 + A(0,1) * B(0,1)
                 + A(1,0) * B(1,0)
                 + A(1,1) * B(1,1)
                 + A(2,0) * B(2,0)
                 + A(2,1) * B(2,1);

          case 15u: // 3x3
            return A(0,0) * B(0,0)
                 + A(0,1) * B(0,1)
                 + A(0,2) * B(0,2)
                 + A(1,0) * B(1,0)
                 + A(1,1) * B(1,1)
                 + A(1,2) * B(1,2)
                 + A(2,0) * B(2,0)
                 + A(2,1) * B(2,1)
                 + A(2,2) * B(2,2);

          default:
            // Should be unreachable with MaxSize==3 and validated sizes.
            assert(false);
            return ScalarType(0);
        }
      }

      template <class EigenDerived>
      constexpr
      ScalarType dot(const Eigen::MatrixBase<EigenDerived>& other) const noexcept
      {
        ScalarType s = ScalarType(0);
        switch (m_rows)
        {
          case 3:
            switch (m_cols)
            {
              case 3:
                s += m_data(2,2) * other(static_cast<Eigen::Index>(2), static_cast<Eigen::Index>(2));
                [[fallthrough]];
              case 2:
                s += m_data(2,1) * other(static_cast<Eigen::Index>(2), static_cast<Eigen::Index>(1));
                [[fallthrough]];
              case 1:
                s += m_data(2,0) * other(static_cast<Eigen::Index>(2), static_cast<Eigen::Index>(0));
                [[fallthrough]];
              case 0:
                break;
              default:
                assert(false);
            }
            [[fallthrough]];
          case 2:
            switch (m_cols)
            {
              case 3:
                s += m_data(1,2) * other(static_cast<Eigen::Index>(1), static_cast<Eigen::Index>(2));
                [[fallthrough]];
              case 2:
                s += m_data(1,1) * other(static_cast<Eigen::Index>(1), static_cast<Eigen::Index>(1));
                [[fallthrough]];
              case 1:
                s += m_data(1,0) * other(static_cast<Eigen::Index>(1), static_cast<Eigen::Index>(0));
                [[fallthrough]];
              case 0:
                break;
              default:
                assert(false);
            }
            [[fallthrough]];
          case 1:
            switch (m_cols)
            {
              case 3:
                s += m_data(0,2) * other(static_cast<Eigen::Index>(0), static_cast<Eigen::Index>(2));
                [[fallthrough]];
              case 2:
                s += m_data(0,1) * other(static_cast<Eigen::Index>(0), static_cast<Eigen::Index>(1));
                [[fallthrough]];
              case 1:
                s += m_data(0,0) * other(static_cast<Eigen::Index>(0), static_cast<Eigen::Index>(0));
                [[fallthrough]];
              case 0:
                break;
              default:
                assert(false);
            }
            [[fallthrough]];
          case 0:
            break;
          default:
            assert(false);
        }
        return s;
      }

      constexpr
      SpatialMatrix<ScalarType> transpose() const noexcept
      {
        SpatialMatrix<ScalarType> T(m_cols, m_rows);
        const auto r = m_rows;
        const auto c = m_cols;
        const auto& A = m_data;

        switch (static_cast<unsigned>(r) * 4u + static_cast<unsigned>(c))
        {
          // any 0×c or r×0: nothing to write
          case 0u: case 1u: case 2u: case 3u:
          case 4u: case 8u: case 12u:
            return T;

          case 5u: // 1x1
            T(0, 0) = A(0, 0); return T;

          case 6u: // 1x2 -> 2x1
            T(0, 0) = A(0, 0);
            T(1, 0) = A(0, 1);
            return T;

          case 7u: // 1x3 -> 3x1
            T(0, 0) = A(0, 0);
            T(1, 0) = A(0, 1);
            T(2, 0) = A(0, 2);
            return T;

          case 9u: // 2x1 -> 1x2
            T(0, 0) = A(0, 0);
            T(0, 1) = A(1, 0);
            return T;

          case 10u: // 2x2
            T(0, 0) = A(0, 0); T(0, 1) = A(1, 0);
            T(1, 0) = A(0, 1); T(1, 1) = A(1, 1);
            return T;

          case 11u: // 2x3 -> 3x2
            T(0, 0) = A(0, 0); T(0, 1) = A(1, 0);
            T(1, 0) = A(0, 1); T(1, 1) = A(1, 1);
            T(2, 0) = A(0, 2); T(2, 1) = A(1, 2);
            return T;

          case 13u: // 3x1 -> 1x3
            T(0, 0) = A(0, 0);
            T(0, 1) = A(1, 0);
            T(0, 2) = A(2, 0);
            return T;

          case 14u: // 3x2 -> 2x3
            T(0, 0) = A(0, 0); T(0, 1) = A(1, 0); T(0, 2) = A(2, 0);
            T(1, 0) = A(0, 1); T(1, 1) = A(1, 1); T(1, 2) = A(2, 1);
            return T;

          case 15u: // 3x3
            T(0, 0) = A(0, 0); T(0, 1) = A(1, 0); T(0, 2) = A(2, 0);
            T(1, 0) = A(0, 1); T(1, 1) = A(1, 1); T(1, 2) = A(2, 1);
            T(2, 0) = A(0, 2); T(2, 1) = A(1, 2); T(2, 2) = A(2, 2);
            return T;

          default:
            assert(false);
            return T;
        }
      }

      constexpr
      SpatialMatrix<ScalarType> adjoint() const noexcept
      {
        SpatialMatrix<ScalarType> Aout(m_cols, m_rows);
        const auto r = m_rows;
        const auto c = m_cols;
        const auto& A = m_data;

        switch (static_cast<unsigned>(r) * 4u + static_cast<unsigned>(c))
        {
          case 0u: case 1u: case 2u: case 3u:
          case 4u: case 8u: case 12u:
            return Aout;

          case 5u: // 1x1
            Aout(0, 0) = conj(A(0, 0)); return Aout;

          case 6u: // 1x2 -> 2x1
            Aout(0, 0) = conj(A(0, 0));
            Aout(1, 0) = conj(A(0, 1));
            return Aout;

          case 7u: // 1x3 -> 3x1
            Aout(0, 0) = conj(A(0, 0));
            Aout(1, 0) = conj(A(0, 1));
            Aout(2, 0) = conj(A(0, 2));
            return Aout;

          case 9u: // 2x1 -> 1x2
            Aout(0, 0) = conj(A(0, 0));
            Aout(0, 1) = conj(A(1, 0));
            return Aout;

          case 10u: // 2x2
            Aout(0, 0) = conj(A(0, 0)); Aout(0, 1) = conj(A(1, 0));
            Aout(1, 0) = conj(A(0, 1)); Aout(1, 1) = conj(A(1, 1));
            return Aout;

          case 11u: // 2x3 -> 3x2
            Aout(0, 0) = conj(A(0, 0)); Aout(0, 1) = conj(A(1, 0));
            Aout(1, 0) = conj(A(0, 1)); Aout(1, 1) = conj(A(1, 1));
            Aout(2, 0) = conj(A(0, 2)); Aout(2, 1) = conj(A(1, 2));
            return Aout;

          case 13u: // 3x1 -> 1x3
            Aout(0, 0) = conj(A(0, 0));
            Aout(0, 1) = conj(A(1, 0));
            Aout(0, 2) = conj(A(2, 0));
            return Aout;

          case 14u: // 3x2 -> 2x3
            Aout(0, 0) = conj(A(0, 0)); Aout(0, 1) = conj(A(1, 0)); Aout(0, 2) = conj(A(2, 0));
            Aout(1, 0) = conj(A(0, 1)); Aout(1, 1) = conj(A(1, 1)); Aout(1, 2) = conj(A(2, 1));
            return Aout;

          case 15u: // 3x3
            Aout(0, 0) = conj(A(0, 0)); Aout(0, 1) = conj(A(1, 0)); Aout(0, 2) = conj(A(2, 0));
            Aout(1, 0) = conj(A(0, 1)); Aout(1, 1) = conj(A(1, 1)); Aout(1, 2) = conj(A(2, 1));
            Aout(2, 0) = conj(A(0, 2)); Aout(2, 1) = conj(A(1, 2)); Aout(2, 2) = conj(A(2, 2));
            return Aout;

          default:
            assert(false);
            return Aout;
        }
      }

      constexpr
      ScalarType trace() const noexcept
      {
        assert(m_rows == m_cols);
        switch (m_rows)
        {
          case 0:
            return ScalarType(0);
          case 1:
            return m_data(0, 0);
          case 2:
            return m_data(0, 0) + m_data(1, 1);
          case 3:
            return m_data(0, 0) + m_data(1, 1) + m_data(2, 2);
          default:
            assert(false);
            return ScalarType(0);
        }
      }

      constexpr
      ScalarType value() const noexcept
      {
        return m_data(0, 0);
      }

      [[nodiscard]] inline
      SpatialVector<ScalarType>
      solve(const SpatialVector<ScalarType>& b) const noexcept
      {
        assert(this->rows() == this->cols());
        assert(b.size() == this->rows());

        const auto n = this->rows();

        SpatialVector<ScalarType> x(n);

        if (n == 1)
        {
          const ScalarType det = (*this)(0,0);
          assert(det != ScalarType(0));
          x(0) = b(0) / det;
          return x;
        }

        if (n == 2)
        {
          const ScalarType a00 = (*this)(0,0), a01 = (*this)(0,1);
          const ScalarType a10 = (*this)(1,0), a11 = (*this)(1,1);

          const ScalarType det = a00*a11 - a01*a10;
          assert(det != ScalarType(0));

          const ScalarType b0 = b(0), b1 = b(1);

          // x = A^{-1} b with explicit 2x2 inverse
          x(0) = ( a11*b0 - a01*b1) / det;
          x(1) = (-a10*b0 + a00*b1) / det;
          return x;
        }

        if (n == 3)
        {
          const ScalarType a00 = (*this)(0,0), a01 = (*this)(0,1), a02 = (*this)(0,2);
          const ScalarType a10 = (*this)(1,0), a11 = (*this)(1,1), a12 = (*this)(1,2);
          const ScalarType a20 = (*this)(2,0), a21 = (*this)(2,1), a22 = (*this)(2,2);

          // Cofactors (of A)
          const ScalarType c00 =  (a11*a22 - a12*a21);
          const ScalarType c01 = -(a10*a22 - a12*a20);
          const ScalarType c02 =  (a10*a21 - a11*a20);

          const ScalarType c10 = -(a01*a22 - a02*a21);
          const ScalarType c11 =  (a00*a22 - a02*a20);
          const ScalarType c12 = -(a00*a21 - a01*a20);

          const ScalarType c20 =  (a01*a12 - a02*a11);
          const ScalarType c21 = -(a00*a12 - a02*a10);
          const ScalarType c22 =  (a00*a11 - a01*a10);

          const ScalarType det = a00*c00 + a01*c01 + a02*c02;
          assert(det != ScalarType(0));

          const ScalarType b0 = b(0), b1 = b(1), b2 = b(2);

          // x = (adj(A)/det) b, where adj(A) = C^T
          x(0) = (c00*b0 + c10*b1 + c20*b2) / det;
          x(1) = (c01*b0 + c11*b1 + c21*b2) / det;
          x(2) = (c02*b0 + c12*b1 + c22*b2) / det;

          return x;
        }

        assert(false && "Unsupported SpatialMatrix size (expected 1..3 square)");
        return x;
      }

      // Determinant specialized for 1x1, 2x2, 3x3 (+ generic fallback).
      constexpr
      ScalarType determinant() const noexcept
      {
        const auto r = m_rows;
        const auto c = m_cols;
        assert(r == c);
        (void) c;

        if (r == 0)
          return ScalarType(1); // convention; you can also assert(false)

        if (r == 1)
        {
          return (*this)(0,0);
        }
        else if (r == 2)
        {
          const ScalarType a = (*this)(0,0);
          const ScalarType b = (*this)(0,1);
          const ScalarType c0 = (*this)(1,0);
          const ScalarType d = (*this)(1,1);
          return a * d - b * c0;
        }
        else if (r == 3)
        {
          const ScalarType a = (*this)(0,0);
          const ScalarType b = (*this)(0,1);
          const ScalarType c0 = (*this)(0,2);
          const ScalarType d = (*this)(1,0);
          const ScalarType e = (*this)(1,1);
          const ScalarType f = (*this)(1,2);
          const ScalarType g = (*this)(2,0);
          const ScalarType h = (*this)(2,1);
          const ScalarType i = (*this)(2,2);
          return a * (e * i - f * h)
               - b * (d * i - f * g)
               + c0 * (d * h - e * g);
        }
        else
        {
          assert(false);
          return Math::nan<ScalarType>();
        }
      }

      // Inverse specialized for 1x1, 2x2, 3x3 (+ generic fallback).
      // Returns a SpatialMatrix (same runtime size).
      constexpr
      SpatialMatrix<ScalarType> inverse() const noexcept
      {
        const auto r = m_rows;
        const auto c = m_cols;
        assert(r == c);

        SpatialMatrix<ScalarType> inv(r, c);

        if (r == 0)
        {
          // convention: inverse of empty is empty
          return inv;
        }

        if (r == 1)
        {
          const ScalarType a = (*this)(0,0);
          assert(a != ScalarType(0));
          inv(0,0) = ScalarType(1) / a;
          return inv;
        }
        else if (r == 2)
        {
          const ScalarType a = (*this)(0,0);
          const ScalarType b = (*this)(0,1);
          const ScalarType c0 = (*this)(1,0);
          const ScalarType d = (*this)(1,1);

          const ScalarType det = a * d - b * c0;
          assert(det != ScalarType(0));
          const ScalarType invdet = ScalarType(1) / det;

          inv(0,0) =  d * invdet;
          inv(0,1) = -b * invdet;
          inv(1,0) = -c0 * invdet;
          inv(1,1) =  a * invdet;
          return inv;
        }
        else if (r == 3)
        {
          const ScalarType a = (*this)(0,0);
          const ScalarType b = (*this)(0,1);
          const ScalarType c0 = (*this)(0,2);
          const ScalarType d = (*this)(1,0);
          const ScalarType e = (*this)(1,1);
          const ScalarType f = (*this)(1,2);
          const ScalarType g = (*this)(2,0);
          const ScalarType h = (*this)(2,1);
          const ScalarType i = (*this)(2,2);

          // cofactors (same as in your PointBase code)
          const ScalarType A =  (e * i - f * h);
          const ScalarType B = -(d * i - f * g);
          const ScalarType C =  (d * h - e * g);
          const ScalarType D = -(b * i - c0 * h);
          const ScalarType E =  (a * i - c0 * g);
          const ScalarType F = -(a * h - b * g);
          const ScalarType G =  (b * f - c0 * e);
          const ScalarType H = -(a * f - c0 * d);
          const ScalarType I =  (a * e - b * d);

          const ScalarType det = a * A + b * B + c0 * C;
          assert(det != ScalarType(0));
          const ScalarType invdet = ScalarType(1) / det;

          // adjugate / det (note the transpose of cofactor matrix)
          inv(0,0) = A * invdet;
          inv(0,1) = D * invdet;
          inv(0,2) = G * invdet;

          inv(1,0) = B * invdet;
          inv(1,1) = E * invdet;
          inv(1,2) = H * invdet;

          inv(2,0) = C * invdet;
          inv(2,1) = F * invdet;
          inv(2,2) = I * invdet;

          return inv;
        }
        else
        {
          assert(false);
          return inv;
        }
      }

      constexpr
      SpatialVector<Scalar> row(std::uint8_t i) const noexcept
      {
        assert(i < m_rows);
        SpatialVector<Scalar> v(m_cols);
        switch (m_cols)
        {
          case 3:
            v(2) = (*this)(static_cast<Eigen::Index>(i), 2);
            [[fallthrough]];
          case 2:
            v(1) = (*this)(static_cast<Eigen::Index>(i), 1);
            [[fallthrough]];
          case 1:
            v(0) = (*this)(static_cast<Eigen::Index>(i), 0);
            break;
          case 0:
            break;
          default:
            assert(false);
        }
        return v;
      }

      constexpr
      SpatialVector<Scalar> col(std::uint8_t j) const noexcept
      {
        assert(j < m_cols);
        SpatialVector<Scalar> v(m_rows);
        switch (m_rows)
        {
          case 3:
            v(2) = (*this)(2, static_cast<Eigen::Index>(j));
            [[fallthrough]];
          case 2:
            v(1) = (*this)(1, static_cast<Eigen::Index>(j));
            [[fallthrough]];
          case 1:
            v(0) = (*this)(0, static_cast<Eigen::Index>(j));
            break;
          case 0:
            break;
          default:
            assert(false);
        }
        return v;
      }

      constexpr
      SpatialMatrix pseudoInverse() const noexcept
      {
        // Moore-Penrose pseudoinverse via (A^T A)^{-1} A^T for full column rank
        // and A^T (A A^T)^{-1} for full row rank.
        SpatialMatrix<Scalar> A_T = this->transpose();
        if (m_rows >= m_cols)
        {
          // Full column rank assumed
          SpatialMatrix<Scalar> AtA = A_T * (*this);
          SpatialMatrix<Scalar> AtA_inv = AtA.inverse();
          return AtA_inv * A_T;
        }
        else
        {
          // Full row rank assumed
          SpatialMatrix<Scalar> AAt = (*this) * A_T;
          SpatialMatrix<Scalar> AAt_inv = AAt.inverse();
          return A_T * AAt_inv;
        }
      }

      constexpr
      auto& getData() noexcept
      {
        return m_data;
      }

      constexpr
      const auto& getData() const noexcept
      {
        return m_data;
      }

    private:
      std::uint8_t m_rows;
      std::uint8_t m_cols;
      Data m_data;
  };

  template <class Scalar>
  [[nodiscard]] inline
  SpatialMatrix<Scalar>
  operator*(const Scalar& s, const SpatialMatrix<Scalar>& A) noexcept
  {
    SpatialMatrix<Scalar> C(A.rows(), A.cols());
    const auto r = A.rows();
    const auto c = A.cols();

    // any 0-dimension => nothing to write
    if (r == 0 || c == 0) return C;

    switch (static_cast<unsigned>(r)  *  4u + static_cast<unsigned>(c))
    {
      case 5u:  C(0, 0) = s * A(0, 0); return C;

      case 6u:  C(0, 0) = s * A(0, 0); C(0, 1) = s * A(0, 1); return C;
      case 7u:  C(0, 0) = s * A(0, 0); C(0, 1) = s * A(0, 1); C(0, 2) = s * A(0, 2); return C;

      case 9u:  C(0, 0) = s * A(0, 0); C(1, 0) = s * A(1, 0); return C;

      case 10u:
        C(0, 0) = s * A(0, 0); C(0, 1) = s * A(0, 1);
        C(1, 0) = s * A(1, 0); C(1, 1) = s * A(1, 1);
        return C;

      case 11u:
        C(0, 0) = s * A(0, 0); C(0, 1) = s * A(0, 1); C(0, 2) = s * A(0, 2);
        C(1, 0) = s * A(1, 0); C(1, 1) = s * A(1, 1); C(1, 2) = s * A(1, 2);
        return C;

      case 13u:
        C(0, 0) = s * A(0, 0); C(1, 0) = s * A(1, 0); C(2, 0) = s * A(2, 0);
        return C;

      case 14u:
        C(0, 0) = s * A(0, 0); C(0, 1) = s * A(0, 1);
        C(1, 0) = s * A(1, 0); C(1, 1) = s * A(1, 1);
        C(2, 0) = s * A(2, 0); C(2, 1) = s * A(2, 1);
        return C;

      case 15u:
        C(0, 0) = s * A(0, 0); C(0, 1) = s * A(0, 1); C(0, 2) = s * A(0, 2);
        C(1, 0) = s * A(1, 0); C(1, 1) = s * A(1, 1); C(1, 2) = s * A(1, 2);
        C(2, 0) = s * A(2, 0); C(2, 1) = s * A(2, 1); C(2, 2) = s * A(2, 2);
        return C;

      default:
        // remaining keys are 0×c / r×0 already returned above
        assert(false);
        return C;
    }
  }

  template <class LHSScalar, class RHSScalar>
  [[nodiscard]] inline
  SpatialMatrix<typename FormLanguage::Mult<LHSScalar, RHSScalar>::Type>
  operator*(const SpatialMatrix<LHSScalar>& A, const Math::Vector<RHSScalar>& s) noexcept
  {
    using OutScalar = typename FormLanguage::Mult<LHSScalar, RHSScalar>::Type;

    assert(A.cols() == static_cast<std::uint8_t>(s.size()));
    assert(A.rows() <= SpatialMatrix<OutScalar>::MaxSize);
    assert(A.cols() <= SpatialMatrix<OutScalar>::MaxSize);
    assert(s.size() <= SpatialMatrix<OutScalar>::MaxSize);

    const std::uint8_t r = A.rows();
    const std::uint8_t k = A.cols();

    SpatialMatrix<OutScalar> C(r, 1);
    C.setZero();

    // Fast unrolled paths for k in {1,2,3} and r in {1,2,3}
    // (k==0 is legal in your type, returns zero column)
    switch (k)
    {
      case 0:
        return C;

      case 1:
      {
        const OutScalar x0 = static_cast<OutScalar>(s(0));
        switch (r)
        {
          case 0: return C;
          case 1: C(0,0) = static_cast<OutScalar>(A(0,0)) * x0; return C;
          case 2:
            C(0,0) = static_cast<OutScalar>(A(0,0)) * x0;
            C(1,0) = static_cast<OutScalar>(A(1,0)) * x0;
            return C;
          case 3:
            C(0,0) = static_cast<OutScalar>(A(0,0)) * x0;
            C(1,0) = static_cast<OutScalar>(A(1,0)) * x0;
            C(2,0) = static_cast<OutScalar>(A(2,0)) * x0;
            return C;
          default:
            assert(false);
            return C;
        }
      }

      case 2:
      {
        const OutScalar x0 = static_cast<OutScalar>(s(0));
        const OutScalar x1 = static_cast<OutScalar>(s(1));
        switch (r)
        {
          case 0: return C;
          case 1:
            C(0,0) = static_cast<OutScalar>(A(0,0)) * x0
                   + static_cast<OutScalar>(A(0,1)) * x1;
            return C;
          case 2:
            C(0,0) = static_cast<OutScalar>(A(0,0)) * x0
                   + static_cast<OutScalar>(A(0,1)) * x1;
            C(1,0) = static_cast<OutScalar>(A(1,0)) * x0
                   + static_cast<OutScalar>(A(1,1)) * x1;
            return C;
          case 3:
            C(0,0) = static_cast<OutScalar>(A(0,0)) * x0
                   + static_cast<OutScalar>(A(0,1)) * x1;
            C(1,0) = static_cast<OutScalar>(A(1,0)) * x0
                   + static_cast<OutScalar>(A(1,1)) * x1;
            C(2,0) = static_cast<OutScalar>(A(2,0)) * x0
                   + static_cast<OutScalar>(A(2,1)) * x1;
            return C;
          default:
            assert(false);
            return C;
        }
      }

      case 3:
      {
        const OutScalar x0 = static_cast<OutScalar>(s(0));
        const OutScalar x1 = static_cast<OutScalar>(s(1));
        const OutScalar x2 = static_cast<OutScalar>(s(2));
        switch (r)
        {
          case 0: return C;
          case 1:
            C(0,0) = static_cast<OutScalar>(A(0,0)) * x0
                   + static_cast<OutScalar>(A(0,1)) * x1
                   + static_cast<OutScalar>(A(0,2)) * x2;
            return C;
          case 2:
            C(0,0) = static_cast<OutScalar>(A(0,0)) * x0
                   + static_cast<OutScalar>(A(0,1)) * x1
                   + static_cast<OutScalar>(A(0,2)) * x2;
            C(1,0) = static_cast<OutScalar>(A(1,0)) * x0
                   + static_cast<OutScalar>(A(1,1)) * x1
                   + static_cast<OutScalar>(A(1,2)) * x2;
            return C;
          case 3:
            C(0,0) = static_cast<OutScalar>(A(0,0)) * x0
                   + static_cast<OutScalar>(A(0,1)) * x1
                   + static_cast<OutScalar>(A(0,2)) * x2;
            C(1,0) = static_cast<OutScalar>(A(1,0)) * x0
                   + static_cast<OutScalar>(A(1,1)) * x1
                   + static_cast<OutScalar>(A(1,2)) * x2;
            C(2,0) = static_cast<OutScalar>(A(2,0)) * x0
                   + static_cast<OutScalar>(A(2,1)) * x1
                   + static_cast<OutScalar>(A(2,2)) * x2;
            return C;
          default:
            assert(false);
            return C;
        }
      }

      default:
        assert(false); // MaxSize==3, so k cannot exceed 3 here
        return C;
    }
  }

  template <class Scalar>
  [[nodiscard]] inline
  SpatialMatrix<Scalar>
  operator*(const SpatialMatrix<Scalar>& A, const Scalar& s) noexcept
  {
    SpatialMatrix<Scalar> C(A.rows(), A.cols());

    const std::uint8_t r = A.rows();
    const std::uint8_t c = A.cols();

    // Any 0-dimension => nothing to write (result is correctly sized)
    if (r == 0 || c == 0)
      return C;

    // key in [0..15] for r, c in [0..3]
    switch (static_cast<unsigned>(r) * 4u + static_cast<unsigned>(c))
    {
      case 5u: // 1x1
        C(0, 0) = A(0, 0) * s;
        return C;

      case 6u: // 1x2
        C(0, 0) = A(0, 0) * s;
        C(0, 1) = A(0, 1) * s;
        return C;

      case 7u: // 1x3
        C(0, 0) = A(0, 0) * s;
        C(0, 1) = A(0, 1) * s;
        C(0, 2) = A(0, 2) * s;
        return C;

      case 9u: // 2x1
        C(0, 0) = A(0, 0) * s;
        C(1, 0) = A(1, 0) * s;
        return C;

      case 10u: // 2x2
        C(0, 0) = A(0, 0) * s; C(0, 1) = A(0, 1) * s;
        C(1, 0) = A(1, 0) * s; C(1, 1) = A(1, 1) * s;
        return C;

      case 11u: // 2x3
        C(0, 0) = A(0, 0) * s; C(0, 1) = A(0, 1) * s; C(0, 2) = A(0, 2) * s;
        C(1, 0) = A(1, 0) * s; C(1, 1) = A(1, 1) * s; C(1, 2) = A(1, 2) * s;
        return C;

      case 13u: // 3x1
        C(0, 0) = A(0, 0) * s;
        C(1, 0) = A(1, 0) * s;
        C(2, 0) = A(2, 0) * s;
        return C;

      case 14u: // 3x2
        C(0, 0) = A(0, 0) * s; C(0, 1) = A(0, 1) * s;
        C(1, 0) = A(1, 0) * s; C(1, 1) = A(1, 1) * s;
        C(2, 0) = A(2, 0) * s; C(2, 1) = A(2, 1) * s;
        return C;

      case 15u: // 3x3
        C(0, 0) = A(0, 0) * s; C(0, 1) = A(0, 1) * s; C(0, 2) = A(0, 2) * s;
        C(1, 0) = A(1, 0) * s; C(1, 1) = A(1, 1) * s; C(1, 2) = A(1, 2) * s;
        C(2, 0) = A(2, 0) * s; C(2, 1) = A(2, 1) * s; C(2, 2) = A(2, 2) * s;
        return C;

      default:
        // Remaining keys are 0×c / r×0,  already returned above.
        assert(false);
        return C;
    }
  }

  template <class Scalar>
  [[nodiscard]] inline
  SpatialMatrix<Scalar>
  operator*(const SpatialMatrix<Scalar>& A, const SpatialMatrix<Scalar>& B)
  {
    assert(A.cols() == B.rows());

    const std::uint8_t r = A.rows();
    const std::uint8_t k = A.cols();
    const std::uint8_t c = B.cols();

    SpatialMatrix<Scalar> C(r, c);

    // Any zero dimension => empty product (correctly sized)
    if (r == 0 || k == 0 || c == 0)
      return C;

    // key in [0..63] for (r,k,c) in [0..3]^3
    const unsigned key = static_cast<unsigned>(r) * 16u
                       + static_cast<unsigned>(k) *  4u
                       + static_cast<unsigned>(c);

    switch (key)
    {
      // -------------------- k = 1 --------------------
      // r x 1  times  1 x c  => outer product

      case 1u * 16u + 1u * 4u + 1u: // 1x1 * 1x1 => 1x1
        C(0, 0) = A(0, 0) * B(0, 0);
        return C;

      case 1u * 16u + 1u * 4u + 2u: // 1x1 * 1x2 => 1x2
        C(0, 0) = A(0, 0) * B(0, 0);
        C(0, 1) = A(0, 0) * B(0, 1);
        return C;

      case 1u * 16u + 1u * 4u + 3u: // 1x1 * 1x3 => 1x3
        C(0, 0) = A(0, 0) * B(0, 0);
        C(0, 1) = A(0, 0) * B(0, 1);
        C(0, 2) = A(0, 0) * B(0, 2);
        return C;

      case 2u * 16u + 1u * 4u + 1u: // 2x1 * 1x1 => 2x1
        C(0, 0) = A(0, 0) * B(0, 0);
        C(1, 0) = A(1, 0) * B(0, 0);
        return C;

      case 2u * 16u + 1u * 4u + 2u: // 2x1 * 1x2 => 2x2
        C(0, 0) = A(0, 0) * B(0, 0);  C(0, 1) = A(0, 0) * B(0, 1);
        C(1, 0) = A(1, 0) * B(0, 0);  C(1, 1) = A(1, 0) * B(0, 1);
        return C;

      case 2u * 16u + 1u * 4u + 3u: // 2x1 * 1x3 => 2x3
        C(0, 0) = A(0, 0) * B(0, 0);  C(0, 1) = A(0, 0) * B(0, 1);  C(0, 2) = A(0, 0) * B(0, 2);
        C(1, 0) = A(1, 0) * B(0, 0);  C(1, 1) = A(1, 0) * B(0, 1);  C(1, 2) = A(1, 0) * B(0, 2);
        return C;

      case 3u * 16u + 1u * 4u + 1u: // 3x1 * 1x1 => 3x1
        C(0, 0) = A(0, 0) * B(0, 0);
        C(1, 0) = A(1, 0) * B(0, 0);
        C(2, 0) = A(2, 0) * B(0, 0);
        return C;

      case 3u * 16u + 1u * 4u + 2u: // 3x1 * 1x2 => 3x2
        C(0, 0) = A(0, 0) * B(0, 0);  C(0, 1) = A(0, 0) * B(0, 1);
        C(1, 0) = A(1, 0) * B(0, 0);  C(1, 1) = A(1, 0) * B(0, 1);
        C(2, 0) = A(2, 0) * B(0, 0);  C(2, 1) = A(2, 0) * B(0, 1);
        return C;

      case 3u * 16u + 1u * 4u + 3u: // 3x1 * 1x3 => 3x3
        C(0, 0) = A(0, 0) * B(0, 0);  C(0, 1) = A(0, 0) * B(0, 1);  C(0, 2) = A(0, 0) * B(0, 2);
        C(1, 0) = A(1, 0) * B(0, 0);  C(1, 1) = A(1, 0) * B(0, 1);  C(1, 2) = A(1, 0) * B(0, 2);
        C(2, 0) = A(2, 0) * B(0, 0);  C(2, 1) = A(2, 0) * B(0, 1);  C(2, 2) = A(2, 0) * B(0, 2);
        return C;

      // -------------------- k = 2 --------------------

      case 1u  *  16u + 2u  *  4u + 1u: // 1x2  *  2x1 => 1x1
        C(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0);
        return C;

      case 1u  *  16u + 2u  *  4u + 2u: // 1x2  *  2x2 => 1x2
        C(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0);
        C(0, 1) = A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1);
        return C;

      case 1u  *  16u + 2u  *  4u + 3u: // 1x2  *  2x3 => 1x3
        C(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0);
        C(0, 1) = A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1);
        C(0, 2) = A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2);
        return C;

      case 2u  *  16u + 2u  *  4u + 1u: // 2x2  *  2x1 => 2x1
        C(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0);
        C(1, 0) = A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0);
        return C;

      case 2u  *  16u + 2u  *  4u + 2u: // 2x2  *  2x2 => 2x2
        C(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0);
        C(0, 1) = A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1);
        C(1, 0) = A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0);
        C(1, 1) = A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1);
        return C;

      case 2u  *  16u + 2u  *  4u + 3u: // 2x2  *  2x3 => 2x3
        C(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0);
        C(0, 1) = A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1);
        C(0, 2) = A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2);

        C(1, 0) = A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0);
        C(1, 1) = A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1);
        C(1, 2) = A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2);
        return C;

      case 3u  *  16u + 2u  *  4u + 1u: // 3x2  *  2x1 => 3x1
        C(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0);
        C(1, 0) = A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0);
        C(2, 0) = A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0);
        return C;

      case 3u  *  16u + 2u  *  4u + 2u: // 3x2  *  2x2 => 3x2
        C(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0);
        C(0, 1) = A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1);

        C(1, 0) = A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0);
        C(1, 1) = A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1);

        C(2, 0) = A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0);
        C(2, 1) = A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1);
        return C;

      case 3u  *  16u + 2u  *  4u + 3u: // 3x2  *  2x3 => 3x3
        C(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0);
        C(0, 1) = A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1);
        C(0, 2) = A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2);

        C(1, 0) = A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0);
        C(1, 1) = A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1);
        C(1, 2) = A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2);

        C(2, 0) = A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0);
        C(2, 1) = A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1);
        C(2, 2) = A(2, 0) * B(0, 2) + A(2, 1) * B(1, 2);
        return C;

      // -------------------- k = 3 --------------------

      case 1u  *  16u + 3u  *  4u + 1u: // 1x3  *  3x1 => 1x1
        C(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0);
        return C;

      case 1u  *  16u + 3u  *  4u + 2u: // 1x3  *  3x2 => 1x2
        C(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0);
        C(0, 1) = A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1);
        return C;

      case 1u  *  16u + 3u  *  4u + 3u: // 1x3  *  3x3 => 1x3
        C(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0);
        C(0, 1) = A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1);
        C(0, 2) = A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2);
        return C;

      case 2u  *  16u + 3u  *  4u + 1u: // 2x3  *  3x1 => 2x1
        C(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0);
        C(1, 0) = A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0);
        return C;

      case 2u  *  16u + 3u  *  4u + 2u: // 2x3  *  3x2 => 2x2
        C(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0);
        C(0, 1) = A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1);

        C(1, 0) = A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0);
        C(1, 1) = A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1);
        return C;

      case 2u  *  16u + 3u  *  4u + 3u: // 2x3  *  3x3 => 2x3
        C(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0);
        C(0, 1) = A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1);
        C(0, 2) = A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2);

        C(1, 0) = A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0);
        C(1, 1) = A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1);
        C(1, 2) = A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2);
        return C;

      case 3u  *  16u + 3u  *  4u + 1u: // 3x3  *  3x1 => 3x1
        C(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0);
        C(1, 0) = A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0);
        C(2, 0) = A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0);
        return C;

      case 3u  *  16u + 3u  *  4u + 2u: // 3x3  *  3x2 => 3x2
        C(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0);
        C(0, 1) = A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1);

        C(1, 0) = A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0);
        C(1, 1) = A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1);

        C(2, 0) = A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0);
        C(2, 1) = A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1);
        return C;

      case 3u  *  16u + 3u  *  4u + 3u: // 3x3  *  3x3 => 3x3
        C(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0);
        C(0, 1) = A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1);
        C(0, 2) = A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2);

        C(1, 0) = A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0);
        C(1, 1) = A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1);
        C(1, 2) = A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2);

        C(2, 0) = A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0);
        C(2, 1) = A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1);
        C(2, 2) = A(2, 0) * B(0, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2);
        return C;

      default:
        assert(false);
        return C;
    }
  }

  template <class Scalar>
  [[nodiscard]] inline
  SpatialMatrix<Scalar>
  operator+(const SpatialMatrix<Scalar>& A, const SpatialMatrix<Scalar>& B)
  {
    assert(A.rows() == B.rows());
    assert(A.cols() == B.cols());

    SpatialMatrix<Scalar> C(A.rows(), A.cols());
    const auto r = A.rows();
    const auto c = A.cols();

    if (r == 0 || c == 0) return C;

    switch (static_cast<unsigned>(r) * 4u + static_cast<unsigned>(c))
    {
      case 5u:
        C(0, 0) = A(0, 0) + B(0, 0); return C;

      case 6u:
        C(0, 0) = A(0, 0) + B(0, 0); C(0, 1) = A(0, 1) + B(0, 1); return C;

      case 7u:
        C(0, 0) = A(0, 0) + B(0, 0); C(0, 1) = A(0, 1) + B(0, 1); C(0, 2) = A(0, 2) + B(0, 2); return C;

      case 9u:
        C(0, 0) = A(0, 0) + B(0, 0); C(1, 0) = A(1, 0) + B(1, 0); return C;

      case 10u:
        C(0, 0) = A(0, 0) + B(0, 0); C(0, 1) = A(0, 1) + B(0, 1);
        C(1, 0) = A(1, 0) + B(1, 0); C(1, 1) = A(1, 1) + B(1, 1);
        return C;

      case 11u:
        C(0, 0) = A(0, 0) + B(0, 0); C(0, 1) = A(0, 1) + B(0, 1); C(0, 2) = A(0, 2) + B(0, 2);
        C(1, 0) = A(1, 0) + B(1, 0); C(1, 1) = A(1, 1) + B(1, 1); C(1, 2) = A(1, 2) + B(1, 2);
        return C;

      case 13u:
        C(0, 0) = A(0, 0) + B(0, 0); C(1, 0) = A(1, 0) + B(1, 0); C(2, 0) = A(2, 0) + B(2, 0);
        return C;

      case 14u:
        C(0, 0) = A(0, 0) + B(0, 0); C(0, 1) = A(0, 1) + B(0, 1);
        C(1, 0) = A(1, 0) + B(1, 0); C(1, 1) = A(1, 1) + B(1, 1);
        C(2, 0) = A(2, 0) + B(2, 0); C(2, 1) = A(2, 1) + B(2, 1);
        return C;

      case 15u:
        C(0, 0) = A(0, 0) + B(0, 0); C(0, 1) = A(0, 1) + B(0, 1); C(0, 2) = A(0, 2) + B(0, 2);
        C(1, 0) = A(1, 0) + B(1, 0); C(1, 1) = A(1, 1) + B(1, 1); C(1, 2) = A(1, 2) + B(1, 2);
        C(2, 0) = A(2, 0) + B(2, 0); C(2, 1) = A(2, 1) + B(2, 1); C(2, 2) = A(2, 2) + B(2, 2);
        return C;

      default:
        assert(false);
        return C;
    }
  }

  template <class LHSScalar, class RHSScalar>
  [[nodiscard]] inline
  SpatialVector<typename FormLanguage::Mult<LHSScalar, RHSScalar>::Type>
  operator*(const SpatialMatrix<LHSScalar>& A, const SpatialVector<RHSScalar>& x) noexcept
  {
    using Out = typename FormLanguage::Mult<LHSScalar, RHSScalar>::Type;

    assert(A.cols() == x.size());

    const std::uint8_t r = A.rows();
    const std::uint8_t c = A.cols();

    SpatialVector<Out> y(r);
    if (r == 0) return y;

    // default init (only needed if your SpatialVector doesn't default-initialize)
    for (std::uint8_t i = 0; i < r; ++i) y[i] = Out(0);

    switch (c)
    {
      case 0:
        return y;

      case 1:
      {
        const Out x0 = static_cast<Out>(x[0]);
        switch (r)
        {
          case 1: y[0] = static_cast<Out>(A(0, 0)) * x0; return y;
          case 2:
            y[0] = static_cast<Out>(A(0, 0)) * x0;
            y[1] = static_cast<Out>(A(1, 0)) * x0;
            return y;
          case 3:
            y[0] = static_cast<Out>(A(0, 0)) * x0;
            y[1] = static_cast<Out>(A(1, 0)) * x0;
            y[2] = static_cast<Out>(A(2, 0)) * x0;
            return y;
          case 0: return y;
          default: assert(false); return y;
        }
      }

      case 2:
      {
        const Out x0 = static_cast<Out>(x[0]);
        const Out x1 = static_cast<Out>(x[1]);
        switch (r)
        {
          case 1:
            y[0] = static_cast<Out>(A(0, 0)) * x0 + static_cast<Out>(A(0, 1)) * x1;
            return y;
          case 2:
            y[0] = static_cast<Out>(A(0, 0)) * x0 + static_cast<Out>(A(0, 1)) * x1;
            y[1] = static_cast<Out>(A(1, 0)) * x0 + static_cast<Out>(A(1, 1)) * x1;
            return y;
          case 3:
            y[0] = static_cast<Out>(A(0, 0)) * x0 + static_cast<Out>(A(0, 1)) * x1;
            y[1] = static_cast<Out>(A(1, 0)) * x0 + static_cast<Out>(A(1, 1)) * x1;
            y[2] = static_cast<Out>(A(2, 0)) * x0 + static_cast<Out>(A(2, 1)) * x1;
            return y;
          case 0: return y;
          default: assert(false); return y;
        }
      }

      case 3:
      {
        const Out x0 = static_cast<Out>(x[0]);
        const Out x1 = static_cast<Out>(x[1]);
        const Out x2 = static_cast<Out>(x[2]);
        switch (r)
        {
          case 1:
            y[0] = static_cast<Out>(A(0, 0)) * x0
                 + static_cast<Out>(A(0, 1)) * x1
                 + static_cast<Out>(A(0, 2)) * x2;
            return y;
          case 2:
            y[0] = static_cast<Out>(A(0, 0)) * x0
                 + static_cast<Out>(A(0, 1)) * x1
                 + static_cast<Out>(A(0, 2)) * x2;
            y[1] = static_cast<Out>(A(1, 0)) * x0
                 + static_cast<Out>(A(1, 1)) * x1
                 + static_cast<Out>(A(1, 2)) * x2;
            return y;
          case 3:
            y[0] = static_cast<Out>(A(0, 0)) * x0
                 + static_cast<Out>(A(0, 1)) * x1
                 + static_cast<Out>(A(0, 2)) * x2;
            y[1] = static_cast<Out>(A(1, 0)) * x0
                 + static_cast<Out>(A(1, 1)) * x1
                 + static_cast<Out>(A(1, 2)) * x2;
            y[2] = static_cast<Out>(A(2, 0)) * x0
                 + static_cast<Out>(A(2, 1)) * x1
                 + static_cast<Out>(A(2, 2)) * x2;
            return y;
          case 0: return y;
          default: assert(false); return y;
        }
      }

      default:
        assert(false);
        return y;
    }
  }

  template <class Scalar>
  [[nodiscard]] inline
  SpatialMatrix<Scalar>
  operator*(const SpatialVector<Scalar>& v, const SpatialMatrix<Scalar>& m) noexcept
  {
    assert(v.size() == m.rows());

    const std::uint8_t r = v.size();
    const std::uint8_t c = m.cols();

    SpatialMatrix<Scalar> result(1, c);

    switch (c)
    {
      case 3:
      {
        Scalar s2 = Scalar(0);
        switch (r)
        {
          case 3: s2 += v[2] * m(2,2); [[fallthrough]];
          case 2: s2 += v[1] * m(1,2); [[fallthrough]];
          case 1: s2 += v[0] * m(0,2); [[fallthrough]];
          case 0: break;
          default: assert(false);
        }
        result(0,2) = s2;
        [[fallthrough]];
      }
      case 2:
      {
        Scalar s1 = Scalar(0);
        switch (r)
        {
          case 3: s1 += v[2] * m(2,1); [[fallthrough]];
          case 2: s1 += v[1] * m(1,1); [[fallthrough]];
          case 1: s1 += v[0] * m(0,1); [[fallthrough]];
          case 0: break;
          default: assert(false);
        }
        result(0,1) = s1;
        [[fallthrough]];
      }
      case 1:
      {
        Scalar s0 = Scalar(0);
        switch (r)
        {
          case 3: s0 += v[2] * m(2,0); [[fallthrough]];
          case 2: s0 += v[1] * m(1,0); [[fallthrough]];
          case 1: s0 += v[0] * m(0,0); [[fallthrough]];
          case 0: break;
          default: assert(false);
        }
        result(0,0) = s0;
        break;
      }
      case 0:
        // 1x0 matrix: nothing to do
        break;
      default:
        assert(false);
    }

    return result;
  }


  template <class EigenDerived, class Scalar>
  [[nodiscard]] inline
  auto operator*(
    const Eigen::MatrixBase<EigenDerived>& s,
    const SpatialMatrix<Scalar>& m)
  {
    return s * m.getData().topLeftCorner(
      static_cast<Eigen::Index>(m.rows()),
      static_cast<Eigen::Index>(m.cols()));
  }

  template <class Scalar, class EigenDerived>
  [[nodiscard]] inline
  auto operator*(
    const SpatialMatrix<Scalar>& m,
    const Eigen::MatrixBase<EigenDerived>& s)
  {
    return m.getData().topLeftCorner(
      static_cast<Eigen::Index>(m.rows()),
      static_cast<Eigen::Index>(m.cols())) * s;
  }

  template <class Scalar>
  std::ostream& operator<<(std::ostream& os, const SpatialMatrix<Scalar>& m)
  {
    os << m.getData().topLeftCorner(
      static_cast<Eigen::Index>(m.rows()),
      static_cast<Eigen::Index>(m.cols()));
    return os;
  }
}

namespace Rodin::FormLanguage
{
  template <class Number>
  struct Traits<Math::SpatialMatrix<Number>>
  {
    using ScalarType = Number;
  };
}
#endif
