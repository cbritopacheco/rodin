/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Vector.h
 * @brief Dense vector classes and dimension aliases.
 *
 * This file provides Rodin vector types with Eigen storage and expression support.
 * Vectors are column-oriented by default and support both dynamic and fixed sizes.
 */
#ifndef RODIN_MATH_VECTOR_H
#define RODIN_MATH_VECTOR_H

#include <utility>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "Rodin/Types.h"

#include "Rodin/FormLanguage/Traits.h"

#include "Common.h"

namespace Eigen::internal
{
  /// @brief Retains Eigen map and reference support for dynamic Rodin vectors.
  template <class Scalar>
  struct traits<Rodin::Math::Vector<Scalar>> : traits<Eigen::VectorX<Scalar>>
  {};

  /// @brief Retains Eigen map and reference support for fixed-size Rodin vectors.
  template <class Scalar, size_t Size>
  struct traits<Rodin::Math::FixedSizeVector<Scalar, Size>>
    : traits<Eigen::Vector<Scalar, Size>>
  {};

  /// @brief Evaluates dynamic Rodin vectors through their Eigen storage.
  template <class Scalar>
  struct evaluator<Rodin::Math::Vector<Scalar>> : evaluator<Eigen::VectorX<Scalar>>
  {
    /// @brief Native Eigen evaluator.
    using Parent = evaluator<Eigen::VectorX<Scalar>>;

    /// @brief Inherits evaluation of the Eigen storage.
    using Parent::Parent;
  };

  /// @brief Evaluates fixed-size Rodin vectors through their Eigen storage.
  template <class Scalar, size_t Size>
  struct evaluator<Rodin::Math::FixedSizeVector<Scalar, Size>>
    : evaluator<Eigen::Vector<Scalar, Size>>
  {
    /// @brief Native Eigen evaluator.
    using Parent = evaluator<Eigen::Vector<Scalar, Size>>;

    /// @brief Inherits evaluation of the Eigen storage.
    using Parent::Parent;
  };
}

namespace Rodin::Math
{
  /**
   * @brief Dynamic-size dense vector type.
   *
   * A column vector with dynamic size determined at runtime. The vector
   * supports standard linear algebra operations with no additional storage.
   * Dot products conjugate the right operand.
   *
   * @tparam ScalarType The element type (e.g., Real, Complex)
   */
  template <class ScalarType>
  class Vector final : public Eigen::VectorX<ScalarType>
  {
    public:
      /// @brief Eigen storage and expression base.
      using Parent = Eigen::VectorX<ScalarType>;

      /// @brief Inherits Eigen's size, coefficient and expression constructors.
      using Parent::Parent;

      /// @brief Assigns Eigen expressions into the owned storage.
      using Parent::operator=;

      /// @brief Constructs a vector with Eigen's default initialization.
      Vector() = default;

      /**
       * @brief Copies the vector and its owned coefficients.
       * @param other Vector to copy.
       */
      Vector(const Vector& other) = default;

      /**
       * @brief Moves the vector and its owned coefficients.
       * @param other Vector to move.
       */
      Vector(Vector&& other) noexcept = default;

      /**
       * @brief Copies Eigen storage into a Rodin vector.
       * @param other Eigen vector to copy.
       */
      Vector(const Parent& other)
        : Parent(other)
      {}

      /**
       * @brief Takes ownership of an Eigen vector's storage.
       * @param other Eigen vector to move.
       */
      Vector(Parent&& other) noexcept
        : Parent(std::move(other))
      {}

      /**
       * @brief Copies the vector coefficients.
       * @param other Vector to copy.
       * @returns This vector.
       */
      Vector& operator=(const Vector& other) = default;

      /**
       * @brief Moves the vector coefficients.
       * @param other Vector to move.
       * @returns This vector.
       */
      Vector& operator=(Vector&& other) noexcept = default;

      /**
       * @brief Contracts with the conjugate of the right operand.
       *
       * Computes @f$\sum_i u_i\overline{v_i}@f$, consistently with
       * @ref Rodin::Math::dot. Eigen expressions and views retain their Eigen
       * types; use @ref Rodin::Math::dot to contract those with this convention.
       *
       * @tparam Derived Eigen type of the right operand.
       * @param other Vector of the same size.
       * @returns The right-conjugating inner product.
       */
      template <class Derived>
      auto dot(const Eigen::MatrixBase<Derived>& other) const
      {
        return Math::dot(static_cast<const Parent&>(*this), other);
      }
  };

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
   * enabling optimizations and stack allocation without additional storage.
   * Dot products conjugate the right operand.
   *
   * @tparam ScalarType The element type
   * @tparam Size The number of elements (must be known at compile time)
   */
  template <class ScalarType, size_t Size>
  class FixedSizeVector final : public Eigen::Vector<ScalarType, Size>
  {
    public:
      /// @brief Eigen storage and expression base.
      using Parent = Eigen::Vector<ScalarType, Size>;

      /// @brief Inherits Eigen's size, coefficient and expression constructors.
      using Parent::Parent;

      /// @brief Assigns Eigen expressions into the owned storage.
      using Parent::operator=;

      /// @brief Constructs a vector with Eigen's default initialization.
      FixedSizeVector() = default;

      /**
       * @brief Copies the vector and its owned coefficients.
       * @param other Vector to copy.
       */
      FixedSizeVector(const FixedSizeVector& other) = default;

      /**
       * @brief Moves the vector and its owned coefficients.
       * @param other Vector to move.
       */
      FixedSizeVector(FixedSizeVector&& other) noexcept = default;

      /**
       * @brief Copies Eigen storage into a Rodin vector.
       * @param other Eigen vector to copy.
       */
      FixedSizeVector(const Parent& other)
        : Parent(other)
      {}

      /**
       * @brief Takes ownership of an Eigen vector's storage.
       * @param other Eigen vector to move.
       */
      FixedSizeVector(Parent&& other) noexcept
        : Parent(std::move(other))
      {}

      /**
       * @brief Copies the vector coefficients.
       * @param other Vector to copy.
       * @returns This vector.
       */
      FixedSizeVector& operator=(const FixedSizeVector& other) = default;

      /**
       * @brief Moves the vector coefficients.
       * @param other Vector to move.
       * @returns This vector.
       */
      FixedSizeVector& operator=(FixedSizeVector&& other) noexcept = default;

      /**
       * @brief Contracts with the conjugate of the right operand.
       *
       * Computes @f$\sum_i u_i\overline{v_i}@f$, consistently with
       * @ref Rodin::Math::dot. Eigen expressions and views retain their Eigen
       * types; use @ref Rodin::Math::dot to contract those with this convention.
       *
       * @tparam Derived Eigen type of the right operand.
       * @param other Vector of the same size.
       * @returns The right-conjugating inner product.
       */
      template <class Derived>
      auto dot(const Eigen::MatrixBase<Derived>& other) const
      {
        return Math::dot(static_cast<const Parent&>(*this), other);
      }
  };

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
}

namespace Rodin::FormLanguage
{
  /// @brief Scalar type of an Eigen column vector, including bounded vectors.
  template <class Number, int Rows, int Options, int MaxRows>
  struct Traits<Eigen::Matrix<Number, Rows, 1, Options, MaxRows, 1>>
  {
    /// @brief Scalar value type.
    using ScalarType = Number;

    /// @brief Compile-time number of elements, or Eigen::Dynamic.
    static constexpr int Size = Rows;
  };

  /// @brief Type traits for a dynamic Math::Vector: exposes the scalar type.
  template <class Number>
  struct Traits<Math::Vector<Number>>
  {
      /// @brief Scalar value type.
      using ScalarType = Number;
  };

  /**
   * @brief Type traits for a Math::FixedSizeVector: exposes the scalar type
   * and the compile-time size.
   */
  template <class Number, size_t S>
  struct Traits<Math::FixedSizeVector<Number, S>>
  {
      /// @brief Scalar value type.
      using ScalarType = Number;

      /// @brief Compile-time number of elements.
      static constexpr size_t Size = S;
  };
}

#endif
