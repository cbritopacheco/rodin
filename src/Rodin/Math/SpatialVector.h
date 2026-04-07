/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file SpatialVector.h
 * @brief Fixed-capacity spatial vector with bounded maximum dimension.
 *
 * This file provides a spatial vector class with maximum dimensions bounded by
 * RODIN_MAXIMAL_SPACE_DIMENSION. Used for geometric points, normals, and other
 * spatial vectors to optimize memory allocation.
 */
#ifndef RODIN_MATH_SPATIALVECTOR_H
#define RODIN_MATH_SPATIALVECTOR_H

#include <iostream>

#include "Rodin/FormLanguage/Traits.h"

#include "ForwardDecls.h"
#include "Common.h"

#include "Traits.h"

namespace Rodin::Math
{
  /**
   * @brief Spatial vector with bounded maximum size.
   *
   * A dynamic-size vector with maximum size bounded by RODIN_MAXIMAL_SPACE_DIMENSION.
   * Used for geometric quantities in 2D or 3D space to optimize memory allocation.
   *
   * @note Internally stores a fixed Eigen::Vector3 regardless of the logical size.
   * For size=1 or size=2 vectors, only the first 1 or 2 elements are active; the
   * remaining elements in the underlying storage are unused. This design avoids
   * dynamic allocation for small spatial vectors at the cost of a few extra bytes.
   *
   * @tparam ScalarType The element type
   */
  template <class ScalarType>
  class SpatialVector
  {
    public:
      using Scalar = ScalarType;
      static constexpr std::uint8_t MaxSize = RODIN_MAXIMAL_SPACE_DIMENSION;

      using Data = Eigen::Vector3<ScalarType>;

      static_assert(MaxSize == 3, "MaxSize must be equal to 3.");

      constexpr
      SpatialVector() noexcept
        : m_size(0)
      {}

      constexpr
      explicit SpatialVector(std::uint8_t size)
        : m_size(size)
      {
        assert(size <= MaxSize);
      }

      constexpr
      SpatialVector(std::initializer_list<ScalarType> init)
        : m_size(init.size())
      {
        assert(init.size() <= MaxSize);
        switch (m_size)
        {
          case 3:
            m_data[2] = *(init.begin() + 2);
            [[fallthrough]];
          case 2:
            m_data[1] = *(init.begin() + 1);
            [[fallthrough]];
          case 1:
            m_data[0] = *(init.begin());
            break;
          case 0:
            break;
          default:
            assert(false);
        }
      }

      template <class EigenDerived>
      constexpr
      SpatialVector(const Eigen::MatrixBase<EigenDerived>& other)
        : m_size(static_cast<std::uint8_t>(other.size()))
      {
        assert(m_size <= MaxSize);
        switch (m_size)
        {
          case 3:
            m_data[2] = other(static_cast<Eigen::Index>(2));
            [[fallthrough]];
          case 2:
            m_data[1] = other(static_cast<Eigen::Index>(1));
            [[fallthrough]];
          case 1:
            m_data[0] = other(static_cast<Eigen::Index>(0));
            break;
          case 0:
            break;
          default:
            assert(false);
        }
      }

      constexpr
      SpatialVector(const SpatialVector& other) noexcept
        : m_size(other.m_size),
          m_data(other.m_data)
      {}

      constexpr
      SpatialVector(SpatialVector&& other) noexcept
        : m_size(std::move(other.m_size)),
          m_data(std::move(other.m_data))
      {}

      constexpr
      SpatialVector& operator=(const SpatialVector& other) noexcept
      {
        if (this != &other)
        {
          m_size = other.m_size;
          m_data = other.m_data;
        }
        return *this;
      }

      constexpr
      SpatialVector& operator=(SpatialVector&& other) noexcept
      {
        if (this != &other)
        {
          m_size = std::move(other.m_size);
          m_data = std::move(other.m_data);
        }
        return *this;
      }

      SpatialVector& operator+=(const SpatialVector& other) noexcept
      {
        assert(m_size == other.m_size);

        switch (m_size)
        {
          case 3:
            m_data[2] += other.m_data[2];
            [[fallthrough]];
          case 2:
            m_data[1] += other.m_data[1];
            [[fallthrough]];
          case 1:
            m_data[0] += other.m_data[0];
            break;
          case 0:
            break;
          default:
            assert(false);
        }

        return *this;
      }

      SpatialVector& operator-=(const SpatialVector& other) noexcept
      {
        assert(m_size == other.m_size);

        switch (m_size)
        {
          case 3:
            m_data[2] -= other.m_data[2];
            [[fallthrough]];
          case 2:
            m_data[1] -= other.m_data[1];
            [[fallthrough]];
          case 1:
            m_data[0] -= other.m_data[0];
            break;
          case 0:
            break;
          default:
            assert(false);
        }

        return *this;
      }

      SpatialVector& operator*=(const ScalarType& s) noexcept
      {
        switch (m_size)
        {
          case 3:
            m_data[2] *= s;
            [[fallthrough]];
          case 2:
            m_data[1] *= s;
            [[fallthrough]];
          case 1:
            m_data[0] *= s;
            break;
          case 0:
            break;
          default:
            assert(false);
        }

        return *this;
      }

      SpatialVector& operator/=(const ScalarType& s) noexcept
      {
        switch (m_size)
        {
          case 3:
            m_data[2] /= s;
            [[fallthrough]];
          case 2:
            m_data[1] /= s;
            [[fallthrough]];
          case 1:
            m_data[0] /= s;
            break;
          case 0:
            break;
          default:
            assert(false);
        }

        return *this;
      }

      SpatialVector operator-() const noexcept
      {
        SpatialVector result(*this);
        result *= ScalarType(-1);
        return result;
      }

      template <class EigenDerived>
      constexpr
      SpatialVector& operator=(const Eigen::ArrayBase<EigenDerived>& other)
      {
        const std::uint8_t n = static_cast<std::uint8_t>(other.size());
        assert(n <= MaxSize);
        m_size = n;
        switch (m_size)
        {
          case 3:
            m_data[2] = other(static_cast<Eigen::Index>(2));
            [[fallthrough]];
          case 2:
            m_data[1] = other(static_cast<Eigen::Index>(1));
            [[fallthrough]];
          case 1:
            m_data[0] = other(static_cast<Eigen::Index>(0));
            break;
          case 0:
            break;
          default:
            assert(false);
        }
        return *this;
      }

      template <class EigenDerived>
      SpatialVector& operator=(const Eigen::MatrixBase<EigenDerived>& v)
      {
        const std::uint8_t n = static_cast<std::uint8_t>(v.size());
        assert(n <= MaxSize);
        m_size = n;
        switch (m_size)
        {
          case 3:
            m_data[2] = v(static_cast<Eigen::Index>(2));
            [[fallthrough]];
          case 2:
            m_data[1] = v(static_cast<Eigen::Index>(1));
            [[fallthrough]];
          case 1:
            m_data[0] = v(static_cast<Eigen::Index>(0));
            break;
          case 0:
            break;
          default:
            assert(false);
        }
        return *this;
      }

      template <class EigenDerived>
      SpatialVector& operator+=(const Eigen::MatrixBase<EigenDerived>& v)
      {
        assert(static_cast<std::uint8_t>(v.size()) == m_size);
        switch (m_size)
        {
          case 3:
            m_data[2] += v(static_cast<Eigen::Index>(2));
            [[fallthrough]];
          case 2:
            m_data[1] += v(static_cast<Eigen::Index>(1));
            [[fallthrough]];
          case 1:
            m_data[0] += v(static_cast<Eigen::Index>(0));
            break;
          case 0:
            break;
          default:
            assert(false);
        }
        return *this;
      }

      template <class EigenDerived>
      SpatialVector& operator-=(const Eigen::MatrixBase<EigenDerived>& v)
      {
        assert(static_cast<std::uint8_t>(v.size()) == m_size);
        switch (m_size)
        {
          case 3:
            m_data[2] -= v(static_cast<Eigen::Index>(2));
            [[fallthrough]];
          case 2:
            m_data[1] -= v(static_cast<Eigen::Index>(1));
            [[fallthrough]];
          case 1:
            m_data[0] -= v(static_cast<Eigen::Index>(0));
            break;
          case 0:
            break;
          default:
            assert(false);
        }
        return *this;
      }

      constexpr
      std::uint8_t size() const noexcept
      {
        return m_size;
      }

      constexpr
      void resize(std::uint8_t n)
      {
        assert(n <= MaxSize);
        m_size = n;
      }

      constexpr
      ScalarType& operator()(std::uint8_t i)
      {
        assert(i < m_size);
        return m_data[i];
      }

      constexpr
      const ScalarType& operator()(std::uint8_t i) const
      {
        assert(i < m_size);
        return m_data[i];
      }

      constexpr
      ScalarType& operator[](std::uint8_t i)
      {
        assert(i < m_size);
        return m_data[i];
      }

      constexpr
      const ScalarType& operator[](std::uint8_t i) const
      {
        assert(i < m_size);
        return m_data[i];
      }

      constexpr
      ScalarType& x()
      {
        assert(m_size >= 1);
        return m_data[0];
      }

      constexpr
      const ScalarType& x() const
      {
        assert(m_size >= 1);
        return m_data[0];
      }

      constexpr
      ScalarType& y()
      {
        assert(m_size >= 2);
        return m_data[1];
      }

      constexpr
      const ScalarType& y() const
      {
        assert(m_size >= 2);
        return m_data[1];
      }

      constexpr
      ScalarType& z()
      {
        assert(m_size >= 3);
        return m_data[2];
      }

      constexpr
      const ScalarType& z() const
      {
        assert(m_size >= 3);
        return m_data[2];
      }

      void setZero() noexcept
      {
        m_data.setZero();
      }

      void setConstant(const ScalarType& value) noexcept
      {
        m_data.setConstant(value);
      }

      [[nodiscard]] constexpr
      SpatialVector cross(const SpatialVector& other) const noexcept
      {
        assert(m_size == 3 && other.m_size == 3);

        SpatialVector r(3);
        r[0] = m_data[1] * other.m_data[2] - m_data[2] * other.m_data[1];
        r[1] = m_data[2] * other.m_data[0] - m_data[0] * other.m_data[2];
        r[2] = m_data[0] * other.m_data[1] - m_data[1] * other.m_data[0];
        return r;
      }

      template <class EigenDerived>
      [[nodiscard]] constexpr
      SpatialVector cross(const Eigen::MatrixBase<EigenDerived>& other) const noexcept
      {
        static_assert(EigenDerived::ColsAtCompileTime == 1 || EigenDerived::RowsAtCompileTime == 1,
                      "cross expects a vector expression");
        assert(m_size == 3);
        assert(other.size() == 3);

        const ScalarType bx = other(static_cast<Eigen::Index>(0));
        const ScalarType by = other(static_cast<Eigen::Index>(1));
        const ScalarType bz = other(static_cast<Eigen::Index>(2));

        SpatialVector r(3);
        r[0] = m_data[1] * bz - m_data[2] * by;
        r[1] = m_data[2] * bx - m_data[0] * bz;
        r[2] = m_data[0] * by - m_data[1] * bx;
        return r;
      }

      inline
      constexpr
      ScalarType dot(const SpatialVector& other) const noexcept
      {
        assert(m_size == other.m_size);

        const auto& a = m_data;
        const auto& b = other.m_data;

        ScalarType s = ScalarType(0);
        switch (m_size)
        {
          case 3:
            s += a[2] * b[2];
            [[fallthrough]];
          case 2:
            s += a[1] * b[1];
            [[fallthrough]];
          case 1:
            s += a[0] * b[0];
            [[fallthrough]];
          case 0:
            break;
          default:
            assert(false);
        }
        return s;
      }

      template <class EigenDerived>
      constexpr
      ScalarType dot(const Eigen::MatrixBase<EigenDerived>& other) const noexcept
      {
        assert(static_cast<std::uint8_t>(other.size()) == m_size);
        ScalarType s = ScalarType(0);
        switch (m_size)
        {
          case 3:
            s += m_data[2] * other(static_cast<Eigen::Index>(2));
            [[fallthrough]];
          case 2:
            s += m_data[1] * other(static_cast<Eigen::Index>(1));
            [[fallthrough]];
          case 1:
            s += m_data[0] * other(static_cast<Eigen::Index>(0));
            [[fallthrough]];
          case 0:
            break;
          default:
            assert(false);
        }
        return s;
      }

      SpatialMatrix<Scalar> transpose() const noexcept
      {
        SpatialMatrix<Scalar> m(1, m_size);
        switch (m_size)
        {
          case 3:
            m(0,2) = m_data[2];
            [[fallthrough]];
          case 2:
            m(0,1) = m_data[1];
            [[fallthrough]];
          case 1:
            m(0,0) = m_data[0];
            break;
          case 0:
            break;
          default:
            assert(false);
        }
        return m;
      }

      ScalarType value() const noexcept
      {
        assert(m_size >= 1);
        return m_data[0];
      }

      constexpr
      void normalize() noexcept
      {
        ScalarType n = this->norm();
        switch (m_size)
        {
          case 3:
            m_data[2] /= n;
            [[fallthrough]];
          case 2:
            m_data[1] /= n;
            [[fallthrough]];
          case 1:
            m_data[0] /= n;
            break;
          case 0:
            break;
          default:
            assert(false);
        }
      }

      constexpr
      auto squaredNorm() const noexcept
      {
        ScalarType s = ScalarType(0);
        switch (m_size)
        {
          case 3:
            s += Math::pow2(m_data[2]);
            [[fallthrough]];
          case 2:
            s += Math::pow2(m_data[1]);
            [[fallthrough]];
          case 1:
            s += Math::pow2(m_data[0]);
            [[fallthrough]];
          case 0:
            break;
          default:
            assert(false);
        }
        return s;
      }

      constexpr
      ScalarType stableNorm() const noexcept
      {
        return m_data.stableNorm();
      }

      constexpr
      ScalarType blueNorm() const noexcept
      {
        return m_data.blueNorm();
      }

      template <size_t P>
      constexpr
      ScalarType lpNorm() const noexcept
      {
        ScalarType s = ScalarType(0);
        std::integral_constant<size_t, P> p;
        switch (m_size)
        {
          case 3:
            s += Math::pow(Math::abs(m_data[2]), p);
            [[fallthrough]];
          case 2:
            s += Math::pow(Math::abs(m_data[1]), p);
            [[fallthrough]];
          case 1:
            s += Math::pow(Math::abs(m_data[0]), p);
            [[fallthrough]];
          case 0:
            break;
          default:
            assert(false);
        }
        return Math::pow(s, ScalarType(1) / ScalarType(P));
      }

      constexpr
      SpatialVector normalized() const noexcept
      {
        SpatialVector v(*this);
        v.normalize();
        return v;
      }

      constexpr
      ScalarType norm() const noexcept
      {
        ScalarType s = 0;
        switch (m_size)
        {
          case 3:
            s += Math::pow2(m_data[2]);
            [[fallthrough]];
          case 2:
            s += Math::pow2(m_data[1]);
            [[fallthrough]];
          case 1:
            s += Math::pow2(m_data[0]);
            [[fallthrough]];
          case 0:
            break;
          default:
            assert(false);
        }
        return Math::sqrt(s);
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

      SpatialVector conjugate() const noexcept
      {
        SpatialVector r(*this);
        if constexpr (std::is_same_v<Scalar, std::complex<float>>
                   || std::is_same_v<Scalar, std::complex<double>>
                   || std::is_same_v<Scalar, std::complex<long double>>)
        {
          switch (m_size)
          {
            case 3:
              r.m_data[2] = std::conj(r.m_data[2]);
              [[fallthrough]];
            case 2:
              r.m_data[1] = std::conj(r.m_data[1]);
              [[fallthrough]];
            case 1:
              r.m_data[0] = std::conj(r.m_data[0]);
              break;
            case 0:
              break;
            default:
              assert(false);
          }
        }
        return r;
      }

      template<class Archive>
      void serialize(Archive& ar, const unsigned int)
      {
        ar & m_size;
        for (std::uint8_t i = 0; i < m_size; i++)
          ar & m_data[i];
      }

    private:
      std::uint8_t m_size;
      Data m_data;
  };

  template <class Scalar>
  [[nodiscard]] inline
  SpatialVector<Scalar>
  operator+(const SpatialVector<Scalar>& a, const SpatialVector<Scalar>& b) noexcept
  {
    assert(a.size() == b.size());
    SpatialVector<Scalar> r(a.size());
    switch (a.size())
    {
      case 3:
        r[2] = a[2] + b[2];
        [[fallthrough]];
      case 2:
        r[1] = a[1] + b[1];
        [[fallthrough]];
      case 1:
        r[0] = a[0] + b[0];
        break;
      case 0:
        break;
      default:
        assert(false);
    }
    return r;
  }

  template <class Scalar>
  [[nodiscard]] inline
  SpatialVector<Scalar>
  operator-(const SpatialVector<Scalar>& a, const SpatialVector<Scalar>& b) noexcept
  {
    assert(a.size() == b.size());
    SpatialVector<Scalar> r(a.size());
    switch (a.size())
    {
      case 3:
        r[2] = a[2] - b[2];
        [[fallthrough]];
      case 2:
        r[1] = a[1] - b[1];
        [[fallthrough]];
      case 1:
        r[0] = a[0] - b[0];
        break;
      case 0:
        break;
      default:
        assert(false);
    }
    return r;
  }

  template <class LHS, class Scalar>
  [[nodiscard]] inline
  SpatialVector<Scalar>
  operator*(const LHS& s, const SpatialVector<Scalar>& v) noexcept
  {
    SpatialVector<Scalar> r(v.size());
    switch (v.size())
    {
      case 3:
        r[2] = s * v[2];
        [[fallthrough]];
      case 2:
        r[1] = s * v[1];
        [[fallthrough]];
      case 1:
        r[0] = s * v[0];
        break;
      case 0:
        break;
      default:
        assert(false);
    }
    return r;
  }

  template <class Scalar, class RHS>
  auto operator*(
      const SpatialVector<Scalar>& v,
      const RHS& s)
  {
    SpatialVector<Scalar> r(v.size());
    switch (v.size())
    {
      case 3:
        r[2] = v[2] * s;
        [[fallthrough]];
      case 2:
        r[1] = v[1] * s;
        [[fallthrough]];
      case 1:
        r[0] = v[0] * s;
        break;
      case 0:
        break;
      default:
        assert(false);
    }
    return r;
  }

  template <class Scalar, class RHS>
  [[nodiscard]] inline
  SpatialVector<Scalar>
  operator/(const SpatialVector<Scalar>& v, const RHS& s) noexcept
  {
    SpatialVector<Scalar> r(v.size());
    switch (v.size())
    {
      case 3:
        r[2] = v[2] / s;
        [[fallthrough]];
      case 2:
        r[1] = v[1] / s;
        [[fallthrough]];
      case 1:
        r[0] = v[0] / s;
        break;
      case 0:
        break;
      default:
        assert(false);
    }
    return r;
  }

  template <class EigenDerived, class Scalar>
  SpatialVector<Scalar> operator+(
    const Eigen::MatrixBase<EigenDerived>& a,
    const SpatialVector<Scalar>& b)
  {
    assert(static_cast<std::uint8_t>(a.size()) == b.size());
    SpatialVector<Scalar> r(b.size());
    switch (b.size())
    {
      case 3:
        r[2] = a(static_cast<Eigen::Index>(2)) + b[2];
        [[fallthrough]];
      case 2:
        r[1] = a(static_cast<Eigen::Index>(1)) + b[1];
        [[fallthrough]];
      case 1:
        r[0] = a(static_cast<Eigen::Index>(0)) + b[0];
        break;
      case 0:
        break;
      default:
        assert(false);
    }
    return r;
  }

  template <class Scalar, class EigenDerived>
  SpatialVector<Scalar> operator+(
    const SpatialVector<Scalar>& a,
    const Eigen::MatrixBase<EigenDerived>& b)
  {
    assert(static_cast<std::uint8_t>(b.size()) == a.size());
    SpatialVector<Scalar> r(a.size());
    switch (a.size())
    {
      case 3:
        r[2] = a[2] + b(static_cast<Eigen::Index>(2));
        [[fallthrough]];
      case 2:
        r[1] = a[1] + b(static_cast<Eigen::Index>(1));
        [[fallthrough]];
      case 1:
        r[0] = a[0] + b(static_cast<Eigen::Index>(0));
        break;
      case 0:
        break;
      default:
        assert(false);
    }
    return r;
  }

  template <class EigenDerived, class Scalar>
  auto operator-(
    const Eigen::MatrixBase<EigenDerived>& a,
    const SpatialVector<Scalar>& b)
  {
    assert(static_cast<std::uint8_t>(a.size()) == b.size());
    return a - b.getData().head(static_cast<Eigen::Index>(b.size()));
  }

  template <class Scalar, class EigenDerived>
  auto operator-(
    const SpatialVector<Scalar>& a,
    const Eigen::MatrixBase<EigenDerived>& b)
  {
    assert(static_cast<std::uint8_t>(b.size()) == a.size());
    return a.getData().head(static_cast<Eigen::Index>(a.size())) - b;
  }

  template <class Scalar, class EigenDerived>
  auto operator*(
      const SpatialVector<Scalar>& v,
      const Eigen::MatrixBase<EigenDerived>& m)
  {
    assert(static_cast<std::uint8_t>(m.rows()) == v.size());
    return v.getData().head(static_cast<Eigen::Index>(v.size())) * m;
  }

  template <class EigenDerived, class Scalar>
  auto operator*(
      const Eigen::MatrixBase<EigenDerived>& m,
      const SpatialVector<Scalar>& v)
  {
    assert(static_cast<std::uint8_t>(m.cols()) == v.size());
    return m * v.getData().head(static_cast<Eigen::Index>(v.size()));
  }

  /**
   * @brief Real-valued spatial vector for point coordinates.
   *
   * Convenience alias for SpatialVector<Real>, commonly used to represent
   * points in 2D or 3D space.
   */
  using SpatialPoint = SpatialVector<Real>;

  template <class Scalar>
  std::ostream& operator<<(std::ostream& os, const SpatialVector<Scalar>& v)
  {
    os << v.getData().head(static_cast<Eigen::Index>(v.size()));
    return os;
  }
}



namespace Rodin::FormLanguage
{
  template <class Number>
  struct Traits<Math::SpatialVector<Number>>
  {
    using ScalarType = Number;
  };
}

#endif
