/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Traits.h
 * @brief Type traits for the Math module and form language.
 *
 * This file defines type traits used by the form language to determine
 * properties of mathematical types, including scalar types and result types
 * of operations. These traits enable compile-time type deduction for
 * variational formulations.
 */
#ifndef RODIN_MATH_TRAITS_H
#define RODIN_MATH_TRAITS_H

#include <type_traits>

#include <Eigen/Core>

#include "Rodin/FormLanguage/Traits.h"
#include "Rodin/Types.h"

#include "Common.h"

namespace Rodin::FormLanguage
{
  /**
   * @brief Type trait to check if a type is an Eigen object.
   *
   * This trait determines whether a type derives from Eigen::EigenBase,
   * which includes matrices, vectors, and array expressions.
   *
   * @tparam T Type to check
   */
  template <class T>
  struct IsEigenObject
  {
    /**
     * @brief True if T is an Eigen object, false otherwise.
     */
    static constexpr bool Value =
      std::is_base_of_v<Eigen::EigenBase<typename std::decay<T>::type>, typename std::decay<T>::type>;
  };

  /**
   * @brief Traits specialization for Boolean type.
   */
  template <>
  struct Traits<Boolean>
  {
    using ScalarType = Boolean;  ///< Scalar type is Boolean itself
  };

  /**
   * @brief Traits specialization for Integer type.
   */
  template <>
  struct Traits<Integer>
  {
    using ScalarType = Integer;  ///< Scalar type is Integer itself
  };

  /**
   * @brief Traits specialization for Real type.
   */
  template <>
  struct Traits<Real>
  {
    using ScalarType = Real;  ///< Scalar type is Real itself
  };

  /**
   * @brief Traits specialization for Complex type.
   */
  template <>
  struct Traits<Complex>
  {
    using ScalarType = Complex;  ///< Scalar type is Complex itself
  };

  /**
   * @brief Type trait for deducing the result type of addition.
   *
   * Determines the type of @f$ \text{LHS} + \text{RHS} @f$ at compile time.
   *
   * @tparam LHS Type of left-hand side operand
   * @tparam RHS Type of right-hand side operand
   */
  template <class LHS, class RHS>
  struct Sum
  {
    /**
     * @brief Result type of @f$ \text{LHS} + \text{RHS} @f$
     */
    using Type = decltype(Math::sum(std::declval<LHS>(), std::declval<RHS>()));
  };

  /**
   * @brief Type trait for deducing the result type of subtraction.
   *
   * Determines the type of @f$ \text{LHS} - \text{RHS} @f$ at compile time.
   *
   * @tparam LHS Type of left-hand side operand
   * @tparam RHS Type of right-hand side operand
   */
  template <class LHS, class RHS>
  struct Minus
  {
    /**
     * @brief Result type of @f$ \text{LHS} - \text{RHS} @f$
     */
    using Type = decltype(Math::minus(std::declval<LHS>(), std::declval<RHS>()));
  };

  /**
   * @brief Type trait for deducing the result type of unary negation.
   *
   * Determines the type of @f$ -\text{Operand} @f$ at compile time.
   *
   * @tparam Operand Type of operand
   */
  template <class Operand>
  struct UnaryMinus
  {
    /**
     * @brief Result type of @f$ -\text{Operand} @f$
     */
    using Type = decltype(Math::minus(std::declval<Operand>()));
  };

  /**
   * @brief Type trait for deducing the result type of multiplication.
   *
   * Determines the type of @f$ \text{LHS} \times \text{RHS} @f$ at compile time.
   *
   * @tparam LHS Type of left-hand side operand
   * @tparam RHS Type of right-hand side operand
   */
  template <class LHS, class RHS>
  struct Mult
  {
    /**
     * @brief Result type of @f$ \text{LHS} \times \text{RHS} @f$
     */
    using Type = decltype(Math::mult(std::declval<LHS>(), std::declval<RHS>()));
  };

  /**
   * @brief Type trait for deducing the result type of division.
   *
   * Determines the type of @f$ \text{LHS} / \text{RHS} @f$ at compile time.
   *
   * @tparam LHS Type of left-hand side operand
   * @tparam RHS Type of right-hand side operand
   */
  template <class LHS, class RHS>
  struct Division
  {
    /**
     * @brief Result type of @f$ \text{LHS} / \text{RHS} @f$
     */
    using Type = decltype(Math::division(std::declval<LHS>(), std::declval<RHS>()));
  };

  /**
   * @brief Type trait for deducing the result type of dot product.
   *
   * Determines the type of @f$ \text{LHS} \cdot \text{RHS} @f$ at compile time.
   *
   * @tparam LHS Type of left-hand side operand
   * @tparam RHS Type of right-hand side operand
   */
  template <class LHS, class RHS>
  struct Dot
  {
    /**
     * @brief Result type of @f$ \text{LHS} \cdot \text{RHS} @f$
     */
    using Type = decltype(Math::dot(std::declval<LHS>(), std::declval<RHS>()));
  };
}

#endif

