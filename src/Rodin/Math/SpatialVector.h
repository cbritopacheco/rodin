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
#include <type_traits>

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
      /// @brief The element scalar type.
      using Scalar = ScalarType;
      /// @brief Maximum storable dimension (RODIN_MAXIMAL_SPACE_DIMENSION, fixed at 3).
      static constexpr std::uint8_t MaxSize = RODIN_MAXIMAL_SPACE_DIMENSION;

      /// @brief Underlying fixed-capacity Eigen storage type.
      using Data = Eigen::Vector3<ScalarType>;

      static_assert(MaxSize == 3, "MaxSize must be equal to 3.");

      /// @brief Constructs an empty (size 0) spatial vector.
      constexpr
      SpatialVector() noexcept
        : m_size(0)
      {
        zeroStorage();
      }

      /// @brief Constructs a zero-initialized spatial vector of the given size.
      /// @param size Number of entries.
      constexpr
      explicit SpatialVector(std::uint8_t size)
        : m_size(size)
      {
        assert(size <= MaxSize);
        zeroStorage();
      }

      /// @brief Constructs a spatial vector from an initializer list of components.
      /// @param init Initial component values.
      constexpr
      SpatialVector(std::initializer_list<ScalarType> init)
        : m_size(init.size())
      {
        assert(init.size() <= MaxSize);
        zeroStorage();
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

      /// @brief Constructs a spatial vector from an Eigen vector expression.
      /// @param other Object to copy from.
      template <class EigenDerived>
      constexpr
      SpatialVector(const Eigen::MatrixBase<EigenDerived>& other)
        : m_size(static_cast<std::uint8_t>(other.size()))
      {
        assert(m_size <= MaxSize);
        zeroStorage();
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

      /// @brief Copy constructor.
      /// @param other Object to copy from.
      constexpr
      SpatialVector(const SpatialVector& other) noexcept
        : m_size(other.m_size),
          m_data(other.m_data)
      {}

      /// @brief Move constructor.
      /// @param other Object to move from.
      constexpr
      SpatialVector(SpatialVector&& other) noexcept
        : m_size(std::move(other.m_size)),
          m_data(std::move(other.m_data))
      {}

      /// @brief Returns a zero spatial vector of the given size.
      /// @param size Number of entries.
      /// @returns A zero spatial vector of the given size.
      static constexpr SpatialVector Zero(std::uint8_t size)
      {
        assert(size <= MaxSize);
        SpatialVector result(size);
        switch (size)
        {
          case 3:
            result.m_data[2] = ScalarType(0);
            [[fallthrough]];
          case 2:
            result.m_data[1] = ScalarType(0);
            [[fallthrough]];
          case 1:
            result.m_data[0] = ScalarType(0);
            break;
          case 0:
            break;
          default:
            assert(false);
        }
        return result;
      }

      /// @brief Copy assignment operator.
      /// @param other Object to copy from.
      /// @returns Reference to this object after the operation.
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

      /// @brief Move assignment operator.
      /// @param other Object to move from.
      /// @returns Reference to this object after the operation.
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

      /// @brief Adds another spatial vector componentwise in place.
      /// @param other Other operand.
      /// @returns Reference to this object after the operation.
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

      /// @brief Subtracts another spatial vector componentwise in place.
      /// @param other Other operand.
      /// @returns Reference to this object after the operation.
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

      /// @brief Scales this vector by a scalar in place.
      /// @param s Scalar factor.
      /// @returns Reference to this object after the operation.
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

      /// @brief Divides this vector by a scalar in place.
      /// @param s Scalar factor.
      /// @returns Reference to this object after the operation.
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

      /// @brief Returns the negation of this vector.
      /// @returns Difference of the operands, or the negated operand for the unary overload.
      SpatialVector operator-() const noexcept
      {
        SpatialVector result(*this);
        result *= ScalarType(-1);
        return result;
      }

      /// @brief Assigns from an Eigen array expression, resizing to match.
      /// @param other Object to copy from.
      /// @returns Reference to this object after the operation.
      template <class EigenDerived>
      constexpr
      SpatialVector& operator=(const Eigen::ArrayBase<EigenDerived>& other)
      {
        const std::uint8_t n = static_cast<std::uint8_t>(other.size());
        assert(n <= MaxSize);
        m_size = n;
        zeroStorage();
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

      /// @brief Assigns from an Eigen vector expression, resizing to match.
      /// @returns Reference to this object after the operation.
      /// @param v Value to assign.
      template <class EigenDerived>
      SpatialVector& operator=(const Eigen::MatrixBase<EigenDerived>& v)
      {
        const std::uint8_t n = static_cast<std::uint8_t>(v.size());
        assert(n <= MaxSize);
        m_size = n;
        zeroStorage();
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

      /// @brief Adds an Eigen vector expression componentwise in place.
      /// @returns Reference to this object after the operation.
      /// @param v Value to assign.
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

      /// @brief Subtracts an Eigen vector expression componentwise in place.
      /// @returns Reference to this object after the operation.
      /// @param v Value to assign.
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

      /// @brief Returns the number of active components.
      /// @returns The number of active components.
      constexpr
      std::uint8_t size() const noexcept
      {
        return m_size;
      }

      /// @brief Sets the logical size (must not exceed MaxSize).
      /// @param n Number of entries.
      constexpr
      void resize(std::uint8_t n)
      {
        assert(n <= MaxSize);
        m_size = n;
      }

      /// @brief Returns a reference to component @p i.
      /// @param i Index of the requested entry.
      /// @returns Reference to the entry at the supplied indices.
      constexpr
      ScalarType& operator()(std::uint8_t i)
      {
        assert(i < m_size);
        return m_data[i];
      }

      /// @brief Returns a const reference to component @p i.
      /// @param i Index of the requested entry.
      /// @returns Reference to the entry at the supplied indices.
      constexpr
      const ScalarType& operator()(std::uint8_t i) const
      {
        assert(i < m_size);
        return m_data[i];
      }

      /// @brief Returns a reference to component @p i.
      /// @param i Index of the requested entry.
      /// @returns Entry at the supplied index.
      constexpr
      ScalarType& operator[](std::uint8_t i)
      {
        assert(i < m_size);
        return m_data[i];
      }

      /// @brief Returns a const reference to component @p i.
      /// @param i Index of the requested entry.
      /// @returns Entry at the supplied index.
      constexpr
      const ScalarType& operator[](std::uint8_t i) const
      {
        assert(i < m_size);
        return m_data[i];
      }

      /// @brief Eigen-compatible element access (for drop-in compatibility with Math::Vector)
      /// @param i Index of the requested entry.
      /// @returns Reference to the entry at the supplied index.
      constexpr
      ScalarType& coeffRef(std::size_t i)
      {
        assert(i < m_size);
        return m_data[static_cast<std::uint8_t>(i)];
      }

      /// @brief Eigen-compatible element access (const, for drop-in compatibility with Math::Vector)
      /// @param i Index of the requested entry.
      /// @returns Reference to the entry at the supplied index.
      constexpr
      const ScalarType& coeffRef(std::size_t i) const
      {
        assert(i < m_size);
        return m_data[static_cast<std::uint8_t>(i)];
      }

      /// @brief Returns a reference to the first (x) component.
      /// @returns A reference to the first (x) component.
      constexpr
      ScalarType& x()
      {
        assert(m_size >= 1);
        return m_data[0];
      }

      /// @brief Returns a const reference to the first (x) component.
      /// @returns A const reference to the first (x) component.
      constexpr
      const ScalarType& x() const
      {
        assert(m_size >= 1);
        return m_data[0];
      }

      /// @brief Returns a reference to the second (y) component.
      /// @returns A reference to the second (y) component.
      constexpr
      ScalarType& y()
      {
        assert(m_size >= 2);
        return m_data[1];
      }

      /// @brief Returns a const reference to the second (y) component.
      /// @returns A const reference to the second (y) component.
      constexpr
      const ScalarType& y() const
      {
        assert(m_size >= 2);
        return m_data[1];
      }

      /// @brief Returns a reference to the third (z) component.
      /// @returns A reference to the third (z) component.
      constexpr
      ScalarType& z()
      {
        assert(m_size >= 3);
        return m_data[2];
      }

      /// @brief Returns a const reference to the third (z) component.
      /// @returns A const reference to the third (z) component.
      constexpr
      const ScalarType& z() const
      {
        assert(m_size >= 3);
        return m_data[2];
      }

      /// @brief Sets all components to zero.
      void setZero() noexcept
      {
        m_data.setZero();
      }

      /// @brief Sets all components to the given value.
      /// @param value Value to store or assign.
      void setConstant(const ScalarType& value) noexcept
      {
        m_data.setConstant(value);
      }

      /// @brief Returns the cross product with another 3D spatial vector.
      /// @param other Other operand.
      /// @returns The cross product with another 3D spatial vector.
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

      /// @brief Returns the cross product with a 3D Eigen vector expression.
      /// @param other Other operand.
      /// @returns The cross product with a 3D Eigen vector expression.
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

      /// @brief Returns the Euclidean dot product with another spatial vector.
      /// @param other Other operand.
      /// @returns The Euclidean dot product with another spatial vector.
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

      /// @brief Returns the Euclidean dot product with an Eigen vector expression.
      /// @param other Other operand.
      /// @returns The Euclidean dot product with an Eigen vector expression.
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

      /// @brief Returns this vector as a 1-by-size row matrix.
      /// @returns This vector as a 1-by-size row matrix.
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

      /// @brief Returns the first component (for size-1 vectors used as scalars).
      /// @returns The first component (for size-1 vectors used as scalars).
      ScalarType value() const noexcept
      {
        assert(m_size >= 1);
        return m_data[0];
      }

      /// @brief Normalizes this vector to unit Euclidean norm in place.
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

      /// @brief Returns the squared Euclidean norm.
      /// @returns The squared Euclidean norm.
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

      /// @brief Returns the Euclidean norm computed with Eigen's stable algorithm.
      /// @returns The Euclidean norm computed with Eigen's stable algorithm.
      constexpr
      ScalarType stableNorm() const noexcept
      {
        return m_data.stableNorm();
      }

      /// @brief Returns the Euclidean norm computed with Eigen's Blue algorithm.
      /// @returns The Euclidean norm computed with Eigen's Blue algorithm.
      constexpr
      ScalarType blueNorm() const noexcept
      {
        return m_data.blueNorm();
      }

      /// @brief Returns the @f$ \ell^P @f$ norm of this vector.
      /// @returns The @f$ \ell^P @f$ norm of this vector.
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

      /// @brief Returns a unit-norm copy of this vector.
      /// @returns A unit-norm copy of this vector.
      constexpr
      SpatialVector normalized() const noexcept
      {
        SpatialVector v(*this);
        v.normalize();
        return v;
      }

      /// @brief Returns the Euclidean norm.
      /// @returns The Euclidean norm.
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

      /// @brief Returns a reference to the underlying Eigen storage.
      /// @returns A reference to the underlying Eigen storage.
      constexpr
      auto& getData() noexcept
      {
        return m_data;
      }

      /// @brief Returns a const reference to the underlying Eigen storage.
      /// @returns A const reference to the underlying Eigen storage.
      constexpr
      const auto& getData() const noexcept
      {
        return m_data;
      }

      /// @brief Returns the complex conjugate (identity for real scalar types).
      /// @returns The complex conjugate (identity for real scalar types).
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

      /// @brief Serializes the vector (for boost::serialization).
      /// @param ar Serialization archive.
      template<class Archive>
      void serialize(Archive& ar, const unsigned int)
      {
        ar & m_size;
        for (std::uint8_t i = 0; i < m_size; i++)
          ar & m_data[i];
      }

    private:
      constexpr void zeroStorage() noexcept
      {
        m_data[0] = ScalarType(0);
        m_data[1] = ScalarType(0);
        m_data[2] = ScalarType(0);
      }

      std::uint8_t m_size;
      Data m_data;
  };

  /// @brief Componentwise sum of two spatial vectors.
  /// @returns Sum of the operands.
  /// @param a Left operand.
  /// @param b Right operand.
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

  /// @brief Componentwise difference of two spatial vectors.
  /// @returns Difference of the operands, or the negated operand for the unary overload.
  /// @param a Left operand.
  /// @param b Right operand.
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

  /// @brief Scalar-times-vector product with real/complex promotion.
  /// @param value Value to store or assign.
  /// @param v Vector operand.
  /// @returns Product of the operands.
  template <class LHS, class Scalar>
    requires(std::is_arithmetic_v<LHS> || std::is_same_v<LHS, Complex>)
  [[nodiscard]] inline auto operator*(
    const LHS& value, const SpatialVector<Scalar>& v) noexcept
  {
    using Coefficient =
      std::conditional_t<std::is_same_v<Scalar, Complex> && std::is_arithmetic_v<LHS>,
        Real, LHS>;
    const Coefficient s = value;
    SpatialVector<typename FormLanguage::Mult<Coefficient, Scalar>::Type> r(v.size());
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

  /// @brief Vector-times-scalar product with real/complex promotion.
  /// @param v Vector operand.
  /// @param value Value to store or assign.
  /// @returns Product of the operands.
  template <class Scalar, class RHS>
    requires(std::is_arithmetic_v<RHS> || std::is_same_v<RHS, Complex>)
  [[nodiscard]] inline auto operator*(const SpatialVector<Scalar>& v, const RHS& value)
  {
    using Coefficient =
      std::conditional_t<std::is_same_v<Scalar, Complex> && std::is_arithmetic_v<RHS>,
        Real, RHS>;
    const Coefficient s = value;
    SpatialVector<typename FormLanguage::Mult<Scalar, Coefficient>::Type> r(v.size());
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

  /// @brief Vector-divided-by-scalar product.
  /// @param v Vector operand.
  /// @param s Scalar factor.
  /// @returns Quotient of the operands.
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

  /// @brief Componentwise sum of an Eigen vector expression and a spatial vector.
  /// @returns Sum of the operands.
  /// @param a Left operand.
  /// @param b Right operand.
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

  /// @brief Componentwise sum of a spatial vector and an Eigen vector expression.
  /// @returns Sum of the operands.
  /// @param a Left operand.
  /// @param b Right operand.
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

  /// @brief Componentwise difference of an Eigen vector expression and a spatial vector.
  /// @returns Difference of the operands, or the negated operand for the unary overload.
  /// @param a Left operand.
  /// @param b Right operand.
  template <class EigenDerived, class Scalar>
  auto operator-(
    const Eigen::MatrixBase<EigenDerived>& a,
    const SpatialVector<Scalar>& b)
  {
    using OutScalar =
      typename FormLanguage::Minus<typename EigenDerived::Scalar, Scalar>::Type;
    assert(static_cast<std::uint8_t>(a.size()) == b.size());
    SpatialVector<OutScalar> result(static_cast<std::uint8_t>(a.size()));
    result = a - b.getData().head(static_cast<Eigen::Index>(b.size()));
    return result;
  }

  /// @brief Componentwise difference of a spatial vector and an Eigen vector expression.
  /// @returns Difference of the operands, or the negated operand for the unary overload.
  /// @param a Left operand.
  /// @param b Right operand.
  template <class Scalar, class EigenDerived>
  auto operator-(
    const SpatialVector<Scalar>& a,
    const Eigen::MatrixBase<EigenDerived>& b)
  {
    using OutScalar =
      typename FormLanguage::Minus<Scalar, typename EigenDerived::Scalar>::Type;
    assert(static_cast<std::uint8_t>(b.size()) == a.size());
    SpatialVector<OutScalar> result(static_cast<std::uint8_t>(b.size()));
    result = a.getData().head(static_cast<Eigen::Index>(a.size())) - b;
    return result;
  }

  /// @brief Row-vector-times-matrix product with an Eigen matrix expression.
  /// @param v Vector operand.
  /// @returns Product of the operands.
  /// @param m Matrix operand.
  template <class Scalar, class EigenDerived>
  auto operator*(
      const SpatialVector<Scalar>& v,
      const Eigen::MatrixBase<EigenDerived>& m)
  {
    using OutScalar =
      typename FormLanguage::Mult<Scalar, typename EigenDerived::Scalar>::Type;
    assert(static_cast<std::uint8_t>(m.rows()) == v.size());
    assert(static_cast<std::uint8_t>(m.cols()) <= SpatialMatrix<OutScalar>::MaxSize);
    SpatialMatrix<OutScalar> result(1, static_cast<std::uint8_t>(m.cols()));
    result = v.getData().head(static_cast<Eigen::Index>(v.size())).transpose() * m;
    return result;
  }

  /// @brief Matrix-times-vector product with an Eigen matrix expression.
  /// @param v Vector operand.
  /// @returns Product of the operands.
  /// @param m Matrix operand.
  template <class EigenDerived, class Scalar>
  auto operator*(
      const Eigen::MatrixBase<EigenDerived>& m,
      const SpatialVector<Scalar>& v)
  {
    using OutScalar =
      typename FormLanguage::Mult<typename EigenDerived::Scalar, Scalar>::Type;
    assert(static_cast<std::uint8_t>(m.cols()) == v.size());
    assert(static_cast<std::uint8_t>(m.rows()) <= SpatialVector<OutScalar>::MaxSize);
    SpatialVector<OutScalar> result(static_cast<std::uint8_t>(m.rows()));
    result = m * v.getData().head(static_cast<Eigen::Index>(v.size()));
    return result;
  }

  /// @brief Componentwise difference of a spatial matrix and an Eigen matrix expression.
  /// @param A System matrix.
  /// @param B Second matrix operand.
  /// @returns Difference of the operands, or the negated operand for the unary overload.
  template <class Scalar, class EigenDerived>
  [[nodiscard]] inline
  SpatialMatrix<Scalar>
  operator-(
    const SpatialMatrix<Scalar>& A,
    const Eigen::MatrixBase<EigenDerived>& B)
  {
    assert(A.rows() == static_cast<std::uint8_t>(B.rows()));
    assert(A.cols() == static_cast<std::uint8_t>(B.cols()));

    SpatialMatrix<Scalar> C(A.rows(), A.cols());

    for (std::uint8_t i = 0; i < A.rows(); ++i)
      for (std::uint8_t j = 0; j < A.cols(); ++j)
        C(i, j) = A(i, j) - B(i, j);

    return C;
  }

  /// @brief Componentwise difference of an Eigen matrix expression and a spatial matrix.
  /// @param A System matrix.
  /// @param B Second matrix operand.
  /// @returns Difference of the operands, or the negated operand for the unary overload.
  template <class EigenDerived, class Scalar>
  [[nodiscard]] inline
  SpatialMatrix<Scalar>
  operator-(
    const Eigen::MatrixBase<EigenDerived>& A,
    const SpatialMatrix<Scalar>& B)
  {
    assert(static_cast<std::uint8_t>(A.rows()) == B.rows());
    assert(static_cast<std::uint8_t>(A.cols()) == B.cols());

    SpatialMatrix<Scalar> C(B.rows(), B.cols());

    for (std::uint8_t i = 0; i < B.rows(); ++i)
      for (std::uint8_t j = 0; j < B.cols(); ++j)
        C(i, j) = A(i, j) - B(i, j);

    return C;
  }

  /**
   * @brief Real-valued spatial vector for point coordinates.
   *
   * Convenience alias for SpatialVector<Real>, commonly used to represent
   * points in 2D or 3D space.
   */
  using SpatialPoint = SpatialVector<Real>;

  /// @brief Streams the vector's active components to an output stream.
  /// @param os Output stream.
  /// @param v Vector operand.
  /// @returns Output stream after writing the object.
  template <class Scalar>
  std::ostream& operator<<(std::ostream& os, const SpatialVector<Scalar>& v)
  {
    os << v.getData().head(static_cast<Eigen::Index>(v.size()));
    return os;
  }
}

namespace Rodin::FormLanguage
{
  /// @brief Type traits for a Math::SpatialVector: exposes the scalar type.
  template <class Number>
  struct Traits<Math::SpatialVector<Number>>
  {
      /// @brief Scalar value type.
      using ScalarType = Number;
  };
}

#endif
