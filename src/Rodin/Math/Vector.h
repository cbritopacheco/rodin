/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Vector.h
 * @brief Dense vector type aliases and definitions.
 *
 * This file provides type aliases for various vector types built on Eigen.
 * Vectors are column-oriented by default and support both dynamic and fixed sizes.
 */
#ifndef RODIN_MATH_VECTOR_H
#define RODIN_MATH_VECTOR_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <cstdint>

#include "Rodin/Math/Common.h"
#include "Rodin/Types.h"

#include "Rodin/FormLanguage/Traits.h"

namespace Rodin::Math
{
  /**
   * @brief Dynamic-size dense vector type.
   *
   * A column vector with dynamic size determined at runtime. The vector
   * supports standard linear algebra operations.
   *
   * @tparam ScalarType The element type (e.g., Real, Complex)
   */
  template <class ScalarType>
  using Vector = Eigen::VectorX<ScalarType>;

  /**
   * @brief Dynamic-size complex-valued vector.
   *
   * Convenience alias for Vector<Complex>.
   */
  using ComplexVector = Vector<Complex>;

  /**
   * @brief Fixed-size vector type.
   *
   * A compile-time fixed-size column vector. The size is known at compile time,
   * enabling optimizations and stack allocation.
   *
   * @tparam ScalarType The element type
   * @tparam Size The number of elements (must be known at compile time)
   */
  template <class ScalarType, size_t Size>
  using FixedSizeVector = Eigen::Vector<ScalarType, Size>;

  /**
   * @brief 2D fixed-size vector.
   *
   * A vector with exactly 2 elements.
   *
   * @tparam ScalarType The element type
   */
  template <class ScalarType>
  using Vector2 = FixedSizeVector<ScalarType, 2>;

  /**
   * @brief 3D fixed-size vector.
   *
   * A vector with exactly 3 elements.
   *
   * @tparam ScalarType The element type
   */
  template <class ScalarType>
  using Vector3 = FixedSizeVector<ScalarType, 3>;

  /**
   * @brief 4D fixed-size vector.
   *
   * A vector with exactly 4 elements.
   *
   * @tparam ScalarType The element type
   */
  template <class ScalarType>
  using Vector4 = FixedSizeVector<ScalarType, 4>;

  /**
   * @brief 8D fixed-size vector.
   *
   * A vector with exactly 8 elements.
   *
   * @tparam ScalarType The element type
   */
  template <class ScalarType>
  using Vector8 = FixedSizeVector<ScalarType, 8>;

  /**
   * @brief 16D fixed-size vector.
   *
   * A vector with exactly 16 elements.
   *
   * @tparam ScalarType The element type
   */
  template <class ScalarType>
  using Vector16 = FixedSizeVector<ScalarType, 16>;

  /**
   * @brief 32D fixed-size vector.
   *
   * A vector with exactly 32 elements.
   *
   * @tparam ScalarType The element type
   */
  template <class ScalarType>
  using Vector32 = FixedSizeVector<ScalarType, 32>;

  /**
   * @brief 64D fixed-size vector.
   *
   * A vector with exactly 64 elements.
   *
   * @tparam ScalarType The element type
   */
  template <class ScalarType>
  using Vector64 = FixedSizeVector<ScalarType, 64>;

  /**
   * @brief 128D fixed-size vector.
   *
   * A vector with exactly 128 elements.
   *
   * @tparam ScalarType The element type
   */
  template <class ScalarType>
  using Vector128 = FixedSizeVector<ScalarType, 128>;

  /**
   * @brief Spatial vector with bounded maximum size.
   *
   * A dynamic-size vector with maximum size bounded by RODIN_MAXIMAL_SPACE_DIMENSION.
   * Used for geometric quantities in 2D or 3D space to optimize memory allocation.
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
        : m_size(other.m_size),
          m_data(std::move(other.m_data))
      {}

      constexpr
      SpatialVector& operator=(const SpatialVector& other) noexcept
      {
        m_size = other.m_size;
        m_data = other.m_data;
        return *this;
      }

      constexpr
      SpatialVector& operator=(SpatialVector&& other) noexcept
      {
        m_size = other.m_size;
        m_data = std::move(other.m_data);
        return *this;
      }

      SpatialVector& operator+=(const SpatialVector& other) noexcept
      {
        assert(m_size == other.m_size);
        this->eigen() += other.eigen();
        return *this;
      }

      SpatialVector& operator-=(const SpatialVector& other) noexcept
      {
        assert(m_size == other.m_size);
        this->eigen() -= other.eigen();
        return *this;
      }

      SpatialVector& operator*=(const ScalarType& s) noexcept
      {
        this->eigen() *= s;
        return *this;
      }

      SpatialVector& operator/=(const ScalarType& s) noexcept
      {
        this->eigen() /= s;
        return *this;
      }

      SpatialVector operator-() const noexcept
      {
        SpatialVector result(*this);
        result *= ScalarType(-1);
        return result;
      }

      SpatialVector& operator*=(const SpatialVector& other)
      {
        assert(m_size == other.m_size);
        this->eigen() *= other.eigen();
        return *this;
      }

      template <class EigenDerived>
      constexpr
      SpatialVector& operator=(const Eigen::ArrayBase<EigenDerived>& other)
      {
        const std::uint8_t n = static_cast<std::uint8_t>(other.size());
        assert(n <= MaxSize);
        m_size = n;
        this->eigen() = other;
        return *this;
      }

      template <class EigenDerived>
      SpatialVector& operator=(const Eigen::MatrixBase<EigenDerived>& v)
      {
        const std::uint8_t n = static_cast<std::uint8_t>(v.size());
        assert(n <= MaxSize);
        m_size = n;
        this->eigen() = v;
        return *this;
      }

      template <class EigenDerived>
      SpatialVector& operator+=(const Eigen::MatrixBase<EigenDerived>& v)
      {
        assert(static_cast<std::uint8_t>(v.size()) == m_size);
        this->eigen() += v;
        return *this;
      }

      template <class EigenDerived>
      SpatialVector& operator-=(const Eigen::MatrixBase<EigenDerived>& v)
      {
        assert(static_cast<std::uint8_t>(v.size()) == m_size);
        this->eigen() -= v;
        return *this;
      }

      template <class EigenDerived>
      SpatialVector& operator*=(const Eigen::MatrixBase<EigenDerived>& s)
      {
        this->eigen() *= s;
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

      auto dot(const SpatialVector& other) const noexcept
      {
        assert(m_size == other.m_size);
        return this->eigen().dot(other.eigen());
      }

      template <class EigenDerived>
      auto dot(const Eigen::MatrixBase<EigenDerived>& other) const noexcept
      {
        assert(m_size == static_cast<std::uint8_t>(other.size()));
        return this->eigen().dot(other);
      }

      auto transpose() const noexcept
      {
        return this->eigen().transpose();
      }

      auto value() const noexcept
      {
        return this->eigen().value();
      }

      constexpr
      void normalize() noexcept
      {
        this->eigen().normalize();
      }

      constexpr
      auto squaredNorm() const noexcept
      {
        return this->eigen().squaredNorm();
      }

      constexpr
      auto stableNorm() const noexcept
      {
        return this->eigen().stableNorm();
      }

      constexpr
      auto blueNorm() const noexcept
      {
        return this->eigen().blueNorm();
      }

      template <size_t P>
      constexpr
      auto lpNorm() const noexcept
      {
        return this->eigen().template lpNorm<P>();
      }

      constexpr
      auto normalized() const noexcept
      {
        return this->eigen().normalized();
      }

      constexpr
      auto norm() const noexcept
      {
        return this->eigen().norm();
      }

      constexpr
      auto data() noexcept
      {
        return m_data.data();
      }

      constexpr
      auto data() const noexcept
      {
        return m_data.data();
      }

      constexpr
      auto eigen() noexcept
      {
        return m_data.head(static_cast<Eigen::Index>(m_size));
      }

      constexpr
      auto eigen() const noexcept
      {
        return m_data.head(static_cast<Eigen::Index>(m_size));
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
  auto operator+(
    const SpatialVector<Scalar>& a,
    const SpatialVector<Scalar>& b)
  {
    assert(a.size() == b.size());
    return a.eigen() + b.eigen();
  }

  template <class Scalar>
  auto operator-(
    const SpatialVector<Scalar>& a,
    const SpatialVector<Scalar>& b)
  {
    assert(a.size() == b.size());
    return a.eigen() - b.eigen();
  }

  template <class LHS, class Scalar>
  auto operator*(
      const LHS& s,
      const SpatialVector<Scalar>& v)
  {
    return s * v.eigen();
  }

  template <class Scalar, class RHS>
  auto operator*(
      const SpatialVector<Scalar>& v,
      const RHS& s)
  {
    return v.eigen() * s;
  }

  template <class EigenDerived, class Scalar>
  auto operator+(
    const Eigen::MatrixBase<EigenDerived>& a,
    const SpatialVector<Scalar>& b)
  {
    assert(static_cast<std::uint8_t>(a.size()) == b.size());
    return a + b.eigen();
  }

  template <class Scalar, class EigenDerived>
  auto operator+(
    const SpatialVector<Scalar>& a,
    const Eigen::MatrixBase<EigenDerived>& b)
  {
    assert(static_cast<std::uint8_t>(b.size()) == a.size());
    return a.eigen() + b;
  }

  template <class EigenDerived, class Scalar>
  auto operator-(
    const Eigen::MatrixBase<EigenDerived>& a,
    const SpatialVector<Scalar>& b)
  {
    assert(static_cast<std::uint8_t>(a.size()) == b.size());
    return a - b.eigen();
  }

  template <class Scalar, class EigenDerived>
  auto operator-(
    const SpatialVector<Scalar>& a,
    const Eigen::MatrixBase<EigenDerived>& b)
  {
    assert(static_cast<std::uint8_t>(b.size()) == a.size());
    return a.eigen() - b;
  }

  template <class Scalar, class EigenDerived>
  auto operator*(
      const SpatialVector<Scalar>& v,
      const Eigen::MatrixBase<EigenDerived>& m)
  {
    assert(static_cast<std::uint8_t>(m.rows()) == v.size());
    return v.eigen() * m;
  }

  template <class EigenDerived, class Scalar>
  auto operator*(
      const Eigen::MatrixBase<EigenDerived>& m,
      const SpatialVector<Scalar>& v)
  {
    assert(static_cast<std::uint8_t>(m.cols()) == v.size());
    return m * v.eigen();
  }

  template <class Scalar>
  auto operator/(const SpatialVector<Scalar>& v, const Scalar& s)
  {
    return v.eigen() / s;
  }

  /**
   * @brief Real-valued spatial vector for point coordinates.
   *
   * Convenience alias for SpatialVector<Real>, commonly used to represent
   * points in 2D or 3D space.
   */
  using SpatialPoint = SpatialVector<Real>;
}

namespace Rodin::FormLanguage
{
  template <class Number>
  struct Traits<Math::Vector<Number>>
  {
    using ScalarType = Number;
  };

  template <class Number>
  struct Traits<Math::SpatialVector<Number>>
  {
    using ScalarType = Number;
  };

  template <class Number, size_t S>
  struct Traits<Math::FixedSizeVector<Number, S>>
  {
    using ScalarType = Number;
    static constexpr size_t Size = S;
  };
}

#endif
