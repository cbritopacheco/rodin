/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Unit.h
 * @brief Base class for type-safe unit representations.
 *
 * This file provides a CRTP base class for creating type-safe unit classes
 * with full arithmetic operations. Units can represent physical quantities
 * like angles, lengths, or other measured values.
 */
#ifndef RODIN_MATH_UNIT_H
#define RODIN_MATH_UNIT_H

namespace Rodin::Math
{
  /**
   * @brief Base class for units using CRTP.
   *
   * This class provides a generic base for creating type-safe units with
   * value semantics and arithmetic operations. It uses the Curiously Recurring
   * Template Pattern (CRTP) to enable derived classes to return their own
   * type from operations.
   *
   * ## Type Safety
   * Units provide compile-time type safety, preventing accidental mixing of
   * incompatible units (e.g., adding radians to meters).
   *
   * ## Arithmetic Operations
   * All standard arithmetic operations are supported:
   * - Addition/subtraction: @f$ u_1 \pm u_2 @f$
   * - Multiplication/division: @f$ u_1 \times u_2 @f$, @f$ u_1 / u_2 @f$
   * - Unary operations: @f$ +u @f$, @f$ -u @f$
   * - Compound assignments: +=, -=, *=, /=
   * - Comparison: ==, !=, <, >, <=, >=
   *
   * ## Example Usage
   * ```cpp
   * // Creating a derived unit class
   * class Meter : public Unit<Meter, Real>
   * {
   *   public:
   *     using Parent = Unit<Meter, Real>;
   *     using Parent::Parent;
   * };
   *
   * // Using the unit
   * Meter length1(5.0);
   * Meter length2(3.0);
   * Meter sum = length1 + length2;  // 8.0 meters
   * ```
   *
   * @tparam Derived The derived class type (CRTP pattern)
   * @tparam T The underlying value type (typically Real)
   *
   * @see Rad
   */
  template <class Derived, class T>
  class Unit
  {
    public:
      using Type = T;

      /**
       * @brief Creates a unit with value 1.
       * @return Unit with value @f$ 1 @f$
       */
      static Derived One()
      {
        return Derived(T(1));
      }

      /**
       * @brief Creates a unit with value 0.
       * @return Unit with value @f$ 0 @f$
       */
      static Derived Zero()
      {
        return Derived(T(0));
      }

      /**
       * @brief Constructs a unit from a value.
       * @param[in] v The value to wrap in the unit
       */
      constexpr
      Unit(T v)
        : m_v(v)
      {}

      constexpr
      Unit(const Unit&) = default;

      constexpr
      Unit(Unit&&) = default;

      constexpr
      Unit& operator=(const Unit&) = default;

      constexpr
      Unit& operator=(Unit&&) = default;

      /**
       * @brief Explicit conversion to underlying type.
       * @return The underlying value
       */
      explicit
      operator T() const
      {
        return m_v;
      }

      /**
       * @brief Equality comparison.
       * @param[in] other Unit to compare with
       * @return True if values are equal
       */
      constexpr
      bool operator==(const Unit& other) const
      {
        return m_v == other.m_v;
      }

      /**
       * @brief Inequality comparison.
       * @param[in] other Unit to compare with
       * @return True if values are not equal
       */
      constexpr
      bool operator!=(const Unit& other) const
      {
        return !operator==(other);
      }

      /**
       * @brief Less-than comparison.
       * @param[in] other Unit to compare with
       * @return True if this value is less than other
       */
      constexpr
      bool operator<(const Unit& other) const
      {
        return m_v < other.m_v;
      }

      /**
       * @brief Greater-than comparison.
       * @param[in] other Unit to compare with
       * @return True if this value is greater than other
       */
      constexpr
      bool operator>(const Unit& other) const
      {
        return m_v > other.m_v;
      }

      /**
       * @brief Less-than-or-equal comparison.
       * @param[in] other Unit to compare with
       * @return True if this value is less than or equal to other
       */
      constexpr
      bool operator<=(const Unit& other) const
      {
        return !operator>(other);
      }

      /**
       * @brief Greater-than-or-equal comparison.
       * @param[in] other Unit to compare with
       * @return True if this value is greater than or equal to other
       */
      constexpr
      bool operator>=(const Unit& other) const
      {
        return !operator<(other);
      }

      /**
       * @brief Addition operator.
       * @param[in] other Unit to add
       * @return Sum of the two units
       */
      constexpr
      auto operator+(const Unit& other) const
      {
        return Unit(m_v + other.m_v);
      }

      /**
       * @brief Subtraction operator.
       * @param[in] other Unit to subtract
       * @return Difference of the two units
       */
      constexpr
      auto operator-(const Unit& other) const
      {
        return Unit(m_v - other.m_v);
      }

      /**
       * @brief Multiplication operator.
       * @param[in] other Unit to multiply by
       * @return Product of the two units
       */
      constexpr
      auto operator*(const Unit& other) const
      {
        return Unit(m_v * other.m_v);
      }


      /**
       * @brief Division operator.
       * @param[in] other Unit to divide by
       * @return Quotient of the two units
       */
      constexpr
      auto operator/(const Unit& other) const
      {
        return Unit(m_v / other.m_v);
      }

      /**
       * @brief Unary plus operator.
       * @return Copy of this unit
       */
      constexpr
      Unit operator+() const
      {
        return *this;
      }

      /**
       * @brief Unary minus operator.
       * @return Negated unit
       */
      constexpr
      Unit operator-() const
      {
        return Unit(-m_v);
      }

      /**
       * @brief Compound addition assignment.
       * @param[in] other Unit to add
       * @return Reference to this unit
       */
      constexpr
      Unit& operator+=(const Unit& other)
      {
        m_v += other.m_v;
        return *this;
      }

      /**
       * @brief Compound subtraction assignment.
       * @param[in] other Unit to subtract
       * @return Reference to this unit
       */
      constexpr
      Unit& operator-=(const Unit& other)
      {
        m_v -= other.m_v;
        return *this;
      }

      /**
       * @brief Compound multiplication assignment.
       * @param[in] other Unit to multiply by
       * @return Reference to this unit
       */
      constexpr
      Unit& operator*=(const Unit& other)
      {
        m_v *= other.m_v;
        return *this;
      }

      /**
       * @brief Compound division assignment.
       * @param[in] other Unit to divide by
       * @return Reference to this unit
       */
      constexpr
      Unit& operator/=(const Unit& other)
      {
        m_v /= other.m_v;
        return *this;
      }

    private:
      T m_v;  ///< The underlying value of the unit
  };
}

#endif

