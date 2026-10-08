/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file SpatialTensor.h
 * @brief Fixed-capacity spatial tensors with explicit runtime axis extents.
 */
#ifndef RODIN_MATH_SPATIALTENSOR_H
#define RODIN_MATH_SPATIALTENSOR_H

#include <array>
#include <cassert>
#include <cmath>
#include <type_traits>

#include "Rodin/Alert/Exception.h"
#include "SpatialMatrix.h"

namespace Rodin::Math
{
  /**
   * @brief A rank-R spatial tensor with stack storage and runtime extents.
   *
   * Each axis has at most RODIN_MAXIMAL_SPACE_DIMENSION entries. Active
   * entries are stored lexicographically, with the last index varying fastest.
   * Rank three represents @f$ G_{ijk}=\partial_k A_{ij} @f$; contracting
   * its last axis with a vector gives @f$ (Gv)_{ij}=\sum_k G_{ijk}v_k @f$.
   * Inner products conjugate the second operand, as in @ref Rodin::Math.
   * Tensor actions on vectors and matrices use ordinary linear contractions.
   *
   * @par Architecture
   * Runtime extents determine the active contiguous prefix of a fixed-size
   * array. Indexing checks every axis; arithmetic checks matching extents.
   * No heap allocation or implicit flattening into a spatial vector occurs.
   */
  template <class ScalarType, size_t Rank>
  class SpatialTensor
  {
    public:
      /// @brief Scalar entry type.
      using Scalar = ScalarType;
      /// @brief Runtime extent of each tensor axis.
      using Extents = std::array<size_t, Rank>;
      /// @brief Maximum extent of a spatial tensor axis.
      static constexpr size_t MaxSize = RODIN_MAXIMAL_SPACE_DIMENSION;
      /// @brief Maximum number of entries stored without allocation.
      static constexpr size_t Capacity = [] {
        size_t count = 1;
        for (size_t i = 0; i < Rank; ++i)
          count *= MaxSize;
        return count;
      }();
      static_assert(Rank >= 3);

      SpatialTensor() = default;
      /**
       * @brief Constructs tensor storage with explicit runtime extents.
       * @param extents Extent of each tensor axis.
       */
      explicit SpatialTensor(const Extents& extents)
      {
        resize(extents);
      }
      /**
       * @brief Constructs tensor storage with explicit runtime extents.
       * @param sizes Extents of the tensor axes.
       */
      template <class... Sizes>
        requires(sizeof...(Sizes) == Rank)
      explicit SpatialTensor(Sizes... sizes)
        : SpatialTensor(Extents{static_cast<size_t>(sizes)...})
      {}

      /**
       * @brief Sets active extents, rejecting axes larger than spatial capacity.
       * @returns Reference to this object after the operation.
       * @param extents Extent of each tensor axis.
       */
      SpatialTensor& resize(const Extents& extents)
      {
        for (auto extent : extents)
          if (extent > MaxSize)
            Alert::Exception() << "SpatialTensor extent exceeds spatial capacity."
                               << Alert::Raise;
        m_extents = extents;
        m_size = 1;
        for (auto extent : extents)
          m_size *= extent;
        return *this;
      }
      /**
       * @brief Returns all active tensor axis extents.
       * @returns All active tensor axis extents.
       */
      const Extents& getExtents() const
      {
        return m_extents;
      }
      /**
       * @brief Returns the active extent of the selected tensor axis.
       * @returns The active extent of the selected tensor axis.
       * @param axis Tensor axis whose extent is requested.
       */
      size_t getDimension(size_t axis) const
      {
        return m_extents.at(axis);
      }
      /**
       * @brief Returns the number of active tensor entries.
       * @returns The number of active tensor entries.
       */
      size_t size() const
      {
        return m_size;
      }
      /**
       * @brief Accesses an active tensor entry with bounds assertions.
       * @param i Index of the requested entry.
       * @returns Entry at the supplied index.
       */
      Scalar& operator[](size_t i)
      {
        assert(i < m_size);
        return m_data[i];
      }
      /**
       * @brief Accesses an active tensor entry with bounds assertions.
       * @param i Index of the requested entry.
       * @returns Entry at the supplied index.
       */
      const Scalar& operator[](size_t i) const
      {
        assert(i < m_size);
        return m_data[i];
      }
      /**
       * @brief Accesses an active tensor entry with bounds assertions.
       * @returns Reference to the entry at the supplied indices.
       * @param indices Index along each tensor axis.
       */
      template <class... Indices>
        requires(sizeof...(Indices) == Rank)
      Scalar& operator()(Indices... indices)
      {
        return m_data[getIndex({static_cast<size_t>(indices)...})];
      }
      /**
       * @brief Accesses an active tensor entry with bounds assertions.
       * @returns Reference to the entry at the supplied indices.
       * @param indices Index along each tensor axis.
       */
      template <class... Indices>
        requires(sizeof...(Indices) == Rank)
      const Scalar& operator()(Indices... indices) const
      {
        return m_data[getIndex({static_cast<size_t>(indices)...})];
      }
      /**
       * @brief Sets every active entry to zero.
       * @returns Reference to this object after the operation.
       */
      SpatialTensor& setZero()
      {
        for (size_t i = 0; i < m_size; ++i)
          m_data[i] = Scalar(0);
        return *this;
      }
      /**
       * @brief Sets every active entry to the supplied constant.
       * @param value Value to store or assign.
       * @returns Reference to this object after the operation.
       */
      SpatialTensor& setConstant(const Scalar& value)
      {
        for (size_t i = 0; i < m_size; ++i)
          m_data[i] = value;
        return *this;
      }
      /**
       * @brief Adds matching active entries after checking extents.
       * @param other Other operand.
       * @returns Reference to this object after the operation.
       */
      SpatialTensor& operator+=(const SpatialTensor& other)
      {
        assert(m_extents == other.m_extents);
        for (size_t i = 0; i < m_size; ++i)
          m_data[i] += other[i];
        return *this;
      }
      /**
       * @brief Subtracts matching entries, or negates active entries.
       * @param other Other operand.
       * @returns Reference to this object after the operation.
       */
      SpatialTensor& operator-=(const SpatialTensor& other)
      {
        assert(m_extents == other.m_extents);
        for (size_t i = 0; i < m_size; ++i)
          m_data[i] -= other[i];
        return *this;
      }
      /**
       * @brief Contracts matching entries and conjugates the second operand.
       * @returns Scalar Frobenius contraction with the supplied tensor.
       * @param other Other operand.
       */
      template <class OtherScalar>
      auto dot(const SpatialTensor<OtherScalar, Rank>& other) const
      {
        assert(m_extents == other.getExtents());
        using Result = decltype(Scalar{} * OtherScalar{});
        Result value = 0;
        for (size_t i = 0; i < m_size; ++i)
          value += m_data[i] * Math::conj(other[i]);
        return value;
      }
      /**
       * @brief Returns the squared Frobenius norm of active entries.
       * @returns The squared Frobenius norm of active entries.
       */
      auto squaredNorm() const
      {
        return std::real(dot(*this));
      }
      /**
       * @brief Returns the Frobenius norm of active entries.
       * @returns The Frobenius norm of active entries.
       */
      auto norm() const
      {
        return std::sqrt(squaredNorm());
      }
      /**
       * @brief Returns entrywise complex conjugation.
       * @returns Entrywise complex conjugation.
       */
      SpatialTensor conjugate() const
      {
        SpatialTensor value(m_extents);
        for (size_t i = 0; i < m_size; ++i)
          value[i] = Math::conj(m_data[i]);
        return value;
      }
      /**
       * @brief Adds matching active entries after checking extents.
       * @param other Other operand.
       * @returns Sum of the operands.
       */
      template <class OtherScalar>
      auto operator+(const SpatialTensor<OtherScalar, Rank>& other) const
      {
        assert(m_extents == other.getExtents());
        SpatialTensor<decltype(Scalar{} + OtherScalar{}), Rank> value(m_extents);
        for (size_t i = 0; i < m_size; ++i)
          value[i] = m_data[i] + other[i];
        return value;
      }
      /**
       * @brief Subtracts matching entries, or negates active entries.
       * @param other Other operand.
       * @returns Difference of the operands, or the negated operand for the unary overload.
       */
      template <class OtherScalar>
      auto operator-(const SpatialTensor<OtherScalar, Rank>& other) const
      {
        assert(m_extents == other.getExtents());
        SpatialTensor<decltype(Scalar{} - OtherScalar{}), Rank> value(m_extents);
        for (size_t i = 0; i < m_size; ++i)
          value[i] = m_data[i] - other[i];
        return value;
      }
      /**
       * @brief Subtracts matching entries, or negates active entries.
       * @returns Difference of the operands, or the negated operand for the unary overload.
       */
      auto operator-() const
      {
        return (*this) * Scalar(-1);
      }
      /**
       * @brief Scales every active tensor entry with scalar type promotion.
       * @returns Product of the operands.
       * @param factor Scalar multiplier.
       */
      template <class Value>
        requires(std::is_arithmetic_v<Value> || std::is_same_v<Value, Complex>)
      auto operator*(const Value& factor) const
      {
        SpatialTensor<std::common_type_t<Scalar, Value>, Rank> value(m_extents);
        for (size_t i = 0; i < m_size; ++i)
          value[i] = static_cast<std::common_type_t<Scalar, Value>>(m_data[i]) *
            static_cast<std::common_type_t<Scalar, Value>>(factor);
        return value;
      }
      /**
       * @brief Scales every active tensor entry with scalar type promotion.
       * @returns Quotient of the operands.
       * @param divisor Scalar divisor.
       */
      template <class Value>
        requires(std::is_arithmetic_v<Value> || std::is_same_v<Value, Complex>)
      auto operator/(const Value& divisor) const
      {
        SpatialTensor<std::common_type_t<Scalar, Value>, Rank> value(m_extents);
        for (size_t i = 0; i < m_size; ++i)
          value[i] = static_cast<std::common_type_t<Scalar, Value>>(m_data[i]) /
            static_cast<std::common_type_t<Scalar, Value>>(divisor);
        return value;
      }
      /**
       * @brief Contracts the last tensor axis with a vector.
       * @returns Product of the operands.
       * @param vector Vector operand.
       */
      template <class OtherScalar>
        requires(Rank == 3)
      auto operator*(const SpatialVector<OtherScalar>& vector) const
      {
        assert(m_extents[2] == vector.size());
        SpatialMatrix<decltype(Scalar{} * OtherScalar{})> value(
          m_extents[0], m_extents[1]);
        value.setZero();
        for (size_t i = 0; i < m_extents[0]; ++i)
          for (size_t j = 0; j < m_extents[1]; ++j)
            for (size_t k = 0; k < m_extents[2]; ++k)
              value(i, j) += (*this)(i, j, k) * vector(k);
        return value;
      }
      /**
       * @brief Rank-four contraction @f$ (CA)_{ij}=\sum_{kl}C_{ijkl}A_{kl} @f$.
       * @returns Product of the operands.
       * @param matrix Matrix operand.
       */
      template <class OtherScalar>
        requires(Rank == 4)
      auto operator*(const SpatialMatrix<OtherScalar>& matrix) const
      {
        assert(m_extents[2] == matrix.rows() && m_extents[3] == matrix.cols());
        SpatialMatrix<std::common_type_t<Scalar, OtherScalar>> value(
          m_extents[0], m_extents[1]);
        value.setZero();
        for (size_t i = 0; i < m_extents[0]; ++i)
          for (size_t j = 0; j < m_extents[1]; ++j)
            for (size_t k = 0; k < m_extents[2]; ++k)
              for (size_t l = 0; l < m_extents[3]; ++l)
                value(i, j) += (*this)(i, j, k, l) * matrix(k, l);
        return value;
      }

    private:
      /**
       * @brief Flattens a tensor multi-index.
       * @param indices Tensor multi-index.
       * @returns Linear storage index corresponding to the supplied tensor coordinates.
       */
      size_t getIndex(const Extents& indices) const
      {
        size_t index = 0;
        for (size_t axis = 0; axis < Rank; ++axis)
        {
          assert(indices[axis] < m_extents[axis]);
          index = index * m_extents[axis] + indices[axis];
        }
        return index;
      }
      Extents m_extents{};
      size_t m_size = 0;
      std::array<Scalar, Capacity> m_data{};
  };
  /**
   * @brief Scales every active tensor entry with scalar type promotion.
   * @returns Product of the operands.
   * @param factor Scalar multiplier.
   * @param tensor Tensor operand.
   */
  template <class Value, class Scalar, size_t Rank>
    requires(std::is_arithmetic_v<Value> || std::is_same_v<Value, Complex>)
  auto operator*(const Value& factor, const SpatialTensor<Scalar, Rank>& tensor)
  {
    return tensor * factor;
  }
}

namespace Rodin::FormLanguage
{
  /// @brief Type traits for the matrix or tensor expression specialization.
  template <class Scalar, size_t Rank>
  struct Traits<Math::SpatialTensor<Scalar, Rank>>
  {
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = Scalar;
  };
}
#endif
