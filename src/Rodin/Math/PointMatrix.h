#ifndef RODIN_MATH_POINTMATRIX_H
#define RODIN_MATH_POINTMATRIX_H

#include <cstdint>
#include <cassert>
#include <initializer_list>
#include <ostream>
#include <utility>

#include <Eigen/Core>

#include "Rodin/FormLanguage/Traits.h"

#include "Common.h"
#include "ForwardDecls.h"
#include "Rodin/Geometry/Point.h"
#include "SpatialVector.h"

namespace Rodin::Math
{
  /**
   * @brief Point matrix with bounded maximum number of rows (space dimension).
   *
   * This type encodes the invariant: rows <= RODIN_MAXIMAL_SPACE_DIMENSION (==3),
   * while allowing an arbitrary number of columns (dynamic).
   *
   * Typical uses:
   * - Point sets stored as columns (3 x N, 2 x N, 1 x N)
   * - Jacobians where the codomain is embedded in R^3 (3 x k, k dynamic)
   *
   * @tparam ScalarType Scalar type.
   */
  class PointMatrix
  {
    public:
      using Scalar = Real;
      static constexpr std::uint8_t MaxRows = RODIN_MAXIMAL_SPACE_DIMENSION;
      static_assert(MaxRows == 3, "RODIN_MAXIMAL_SPACE_DIMENSION must be equal to 3.");

      // Storage is always 3 x cols (rows are a runtime prefix 0..3)
      using Data = Eigen::Matrix<Scalar, MaxRows, Eigen::Dynamic>;

      PointMatrix() noexcept
        : m_rows(0),
          m_cols(0),
          m_data(MaxRows, 0)
      {}

      explicit
      PointMatrix(std::uint8_t rows, Eigen::Index cols)
        : m_rows(rows),
          m_cols(cols),
          m_data(MaxRows, cols)
      {
        assert(rows <= MaxRows);
        assert(cols >= 0);
      }

      PointMatrix(const PointMatrix&) = default;

      PointMatrix(PointMatrix&&) = default;

      PointMatrix& operator=(const PointMatrix&) = default;

      PointMatrix& operator=(PointMatrix&&) = default;

      template <class EigenDerived>
      explicit
      PointMatrix(const Eigen::MatrixBase<EigenDerived>& other)
        : m_rows(static_cast<std::uint8_t>(other.rows())),
          m_cols(static_cast<Eigen::Index>(other.cols())),
          m_data(MaxRows, static_cast<Eigen::Index>(other.cols()))
      {
        assert(m_rows <= MaxRows);
        this->getData().topRows(static_cast<Eigen::Index>(m_rows)).leftCols(m_cols) = other;
      }

      template <class EigenDerived>
      PointMatrix& operator=(const Eigen::MatrixBase<EigenDerived>& other)
      {
        const auto r = static_cast<std::uint8_t>(other.rows());
        const auto c = static_cast<Eigen::Index>(other.cols());
        assert(r <= MaxRows);
        this->resize(r, c);
        this->getData().topRows(static_cast<Eigen::Index>(m_rows)).leftCols(m_cols) = other;
        return *this;
      }

      template <class EigenDerived>
      PointMatrix& operator=(const Eigen::ArrayBase<EigenDerived>& other)
      {
        const auto r = static_cast<std::uint8_t>(other.rows());
        const auto c = static_cast<Eigen::Index>(other.cols());
        assert(r <= MaxRows);
        this->resize(r, c);
        this->getData().topRows(static_cast<Eigen::Index>(m_rows)).leftCols(m_cols) = other.matrix();
        return *this;
      }

      PointMatrix& operator*=(const Scalar& s) noexcept
      {
        m_data.topRows(static_cast<Eigen::Index>(m_rows)).leftCols(m_cols) *= s;
        return *this;
      }

      [[nodiscard]] constexpr
      std::uint8_t rows() const noexcept
      {
        return m_rows;
      }

      [[nodiscard]] constexpr
      Eigen::Index cols() const noexcept
      {
        return m_cols;
      }

      void resize(std::uint8_t r, Eigen::Index c)
      {
        assert(r <= MaxRows);
        assert(c >= 0);
        m_rows = r;
        m_cols = c;
        // Keep 3 rows allocated; resize only columns.
        m_data.conservativeResize(MaxRows, c);
      }

      [[nodiscard]] inline
      Scalar& operator()(std::uint8_t i, Eigen::Index j) noexcept
      {
        assert(i < m_rows);
        assert(j >= 0 && j < m_cols);
        return m_data(static_cast<Eigen::Index>(i), j);
      }

      [[nodiscard]] inline
      const Scalar& operator()(std::uint8_t i, Eigen::Index j) const noexcept
      {
        assert(i < m_rows);
        assert(j >= 0 && j < m_cols);
        return m_data(static_cast<Eigen::Index>(i), j);
      }

      void setZero() noexcept
      {
        m_data.setZero();
      }

      void setConstant(const Scalar& v) noexcept
      {
        m_data.setConstant(v);
      }

      // Identity only makes sense for square with rows==cols<=3; columns dynamic so we guard.
      void setIdentity() noexcept
      {
        assert(m_rows == static_cast<std::uint8_t>(m_cols));
        assert(m_rows <= MaxRows);
        m_data.setZero();
        for (std::uint8_t i = 0; i < m_rows; ++i)
          (*this)(i, static_cast<Eigen::Index>(i)) = Scalar(1);
      }

      [[nodiscard]] inline
      Data& getData() noexcept
      {
        return m_data;
      }

      [[nodiscard]] inline
      const Data& getData() const noexcept
      {
        return m_data;
      }

      /**
       * @brief Frobenius dot product with another PointMatrix of same runtime size.
       */
      [[nodiscard]] inline
      Scalar dot(const PointMatrix& other) const noexcept
      {
        assert(m_rows == other.m_rows);
        assert(m_cols == other.m_cols);

        if (m_rows == 0 || m_cols == 0)
          return Scalar(0);

        const auto c = m_cols;
        Scalar s = Scalar(0);

        // Unroll rows (<=3); loop columns (possibly large).
        switch (m_rows)
        {
          case 3:
            for (Eigen::Index j = 0; j < c; ++j)
              s += m_data(2, j) * other.m_data(2, j);
            [[fallthrough]];
          case 2:
            for (Eigen::Index j = 0; j < c; ++j)
              s += m_data(1, j) * other.m_data(1, j);
            [[fallthrough]];
          case 1:
            for (Eigen::Index j = 0; j < c; ++j)
              s += m_data(0, j) * other.m_data(0, j);
            break;
          case 0:
            break;
          default:
            assert(false);
        }

        return s;
      }

      template <class EigenDerived>
      [[nodiscard]] inline
      Scalar dot(const Eigen::MatrixBase<EigenDerived>& other) const noexcept
      {
        assert(static_cast<std::uint8_t>(other.rows()) == m_rows);
        assert(static_cast<Eigen::Index>(other.cols()) == m_cols);

        if (m_rows == 0 || m_cols == 0)
          return Scalar(0);

        const auto c = m_cols;
        Scalar s = Scalar(0);

        switch (m_rows)
        {
          case 3:
            for (Eigen::Index j = 0; j < c; ++j)
              s += m_data(2, j) * other(static_cast<Eigen::Index>(2), j);
            [[fallthrough]];
          case 2:
            for (Eigen::Index j = 0; j < c; ++j)
              s += m_data(1, j) * other(static_cast<Eigen::Index>(1), j);
            [[fallthrough]];
          case 1:
            for (Eigen::Index j = 0; j < c; ++j)
              s += m_data(0, j) * other(static_cast<Eigen::Index>(0), j);
            break;
          case 0:
            break;
          default:
            assert(false);
        }

        return s;
      }

      [[nodiscard]] inline
      Scalar squaredNorm() const noexcept
      {
        if (m_rows == 0 || m_cols == 0)
          return Scalar(0);

        const auto c = m_cols;
        Scalar s = Scalar(0);

        switch (m_rows)
        {
          case 3:
            for (Eigen::Index j = 0; j < c; ++j) s += Math::pow2(m_data(2, j));
            [[fallthrough]];
          case 2:
            for (Eigen::Index j = 0; j < c; ++j) s += Math::pow2(m_data(1, j));
            [[fallthrough]];
          case 1:
            for (Eigen::Index j = 0; j < c; ++j) s += Math::pow2(m_data(0, j));
            break;
          case 0:
            break;
          default:
            assert(false);
        }

        return s;
      }

      /**
       * @brief Transpose (cols x rows). Since rows<=3, result has at most 3 columns.
       */
      [[nodiscard]] inline
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, 0, Eigen::Dynamic, MaxRows>
      transpose() const
      {
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, 0, Eigen::Dynamic, MaxRows> T;
        T.resize(m_cols, static_cast<Eigen::Index>(m_rows));

        if (m_rows == 0 || m_cols == 0)
          return T;

        for (Eigen::Index j = 0; j < m_cols; ++j)
        {
          switch (m_rows)
          {
            case 3: T(j, 2) = m_data(2, j); [[fallthrough]];
            case 2: T(j, 1) = m_data(1, j); [[fallthrough]];
            case 1: T(j, 0) = m_data(0, j); break;
            case 0: break;
            default: assert(false);
          }
        }

        return T;
      }

      /**
       * @brief Adjoint (conjugate transpose), consistent with Eigen semantics.
       */
      [[nodiscard]] inline
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, 0, Eigen::Dynamic, MaxRows>
      adjoint() const
      {
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, 0, Eigen::Dynamic, MaxRows> A;
        A.resize(m_cols, static_cast<Eigen::Index>(m_rows));

        if (m_rows == 0 || m_cols == 0)
          return A;

        for (Eigen::Index j = 0; j < m_cols; ++j)
        {
          switch (m_rows)
          {
            case 3: A(j, 2) = conj(m_data(2, j)); [[fallthrough]];
            case 2: A(j, 1) = conj(m_data(1, j)); [[fallthrough]];
            case 1: A(j, 0) = conj(m_data(0, j)); break;
            case 0: break;
            default: assert(false);
          }
        }

        return A;
      }

      [[nodiscard]] inline
      Math::SpatialVector<Scalar> col(Eigen::Index j) const noexcept
      {
        assert(j >= 0 && j < m_cols);
        Math::SpatialVector<Scalar> v(m_rows);
        switch (m_rows)
        {
          case 3:
            v(2) = m_data(2, j);
            [[fallthrough]];
          case 2:
            v(1) = m_data(1, j);
            [[fallthrough]];
          case 1:
            v(0) = m_data(0, j);
            break;
          case 0:
            break;
          default:
            assert(false);
        }
        return v;
      }

      template<class Archive>
      void serialize(Archive& ar, const unsigned int)
      {
        ar & m_rows;
        ar & m_cols;
        ar & m_data;
      }

    private:
      std::uint8_t m_rows;
      Eigen::Index m_cols;
      Data m_data;
  };

  // -------- free operators (minimal but complete for typical usage) --------

  [[nodiscard]] inline
  PointMatrix operator*(const Real& s, const PointMatrix& A)
  {
    PointMatrix C(A.rows(), A.cols());
    if (A.rows() == 0 || A.cols() == 0)
      return C;

    // Multiply only active block; keep full storage consistent.
    C.getData().leftCols(A.cols()) = s * A.getData().leftCols(A.cols());
    return C;
  }

  [[nodiscard]] inline
  PointMatrix operator*(const PointMatrix& A, const Real& s)
  {
    return s * A;
  }

  template <class Scalar>
  [[nodiscard]] inline
  PointMatrix operator+(const PointMatrix& A, const PointMatrix& B)
  {
    assert(A.rows() == B.rows());
    assert(A.cols() == B.cols());

    PointMatrix C(A.rows(), A.cols());
    if (A.rows() == 0 || A.cols() == 0)
      return C;

    C.getData().leftCols(A.cols()) = A.getData().leftCols(A.cols()) + B.getData().leftCols(A.cols());
    return C;
  }

  template <class Scalar>
  [[nodiscard]] inline
  PointMatrix operator-(const PointMatrix& A, const PointMatrix& B)
  {
    assert(A.rows() == B.rows());
    assert(A.cols() == B.cols());

    PointMatrix C(A.rows(), A.cols());
    if (A.rows() == 0 || A.cols() == 0)
      return C;

    C.getData().leftCols(A.cols()) = A.getData().leftCols(A.cols()) - B.getData().leftCols(A.cols());
    return C;
  }

  // PointMatrix (r x k) times Eigen matrix (k x n) -> Eigen result (r x n)
  template <class Scalar, class EigenDerived>
  [[nodiscard]] inline
  auto operator*(const PointMatrix& A, const Eigen::MatrixBase<EigenDerived>& B)
  {
    assert(static_cast<Eigen::Index>(A.cols()) == B.rows());
    return A.getData().topRows(static_cast<Eigen::Index>(A.rows())).leftCols(A.cols()) * B;
  }

  // Eigen matrix (m x r) times PointMatrix (r x k) -> Eigen result (m x k)
  template <class EigenDerived, class Scalar>
  [[nodiscard]] inline
  auto operator*(const Eigen::MatrixBase<EigenDerived>& A, const PointMatrix& B)
  {
    assert(A.cols() == static_cast<Eigen::Index>(B.rows()));
    return A * B.getData()
      .topRows(static_cast<Eigen::Index>(B.rows()))
      .leftCols(B.cols());
  }

  template <class Scalar>
  inline
  std::ostream& operator<<(std::ostream& os, const PointMatrix& M)
  {
    os << M.getData().topRows(static_cast<Eigen::Index>(M.rows())).leftCols(M.cols());
    return os;
  }
}

namespace Rodin::FormLanguage
{
  template <>
  struct Traits<Math::PointMatrix>
  {
    using ScalarType = Real;
  };
}

#endif
