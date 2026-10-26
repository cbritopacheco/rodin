/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Compile-time regression for abstract function-base type deduction.
 */
#include <gtest/gtest.h>
#include "Rodin/Variational/RealFunction.h"
#include "Rodin/Variational/ComplexFunction.h"

namespace Rodin::Tests::Unit
{
  /**
   * @brief Deduces abstract CRTP parameters without constructing their values.
   *
   * Scalar expressions inherit FunctionBase through an abstract intermediate
   * class. The trait must identify that class, not the concrete expression,
   * and remain valid when the expression is supplied through a const reference.
   */
  TEST(FormLanguage_Traits, AbstractFunctionParameter)
  {
    using Expression = Variational::RealFunction<Real>;
    using Expected =
      Variational::ScalarFunctionBase<Real, Variational::RealFunctionBase<Expression>>;
    static_assert(std::is_abstract_v<Expected>);
    static_assert(
      std::is_same_v<FormLanguage::FunctionDerived<Expression>::Type, Expected>);
    static_assert(
      std::is_same_v<FormLanguage::FunctionDerived<const Expression&>::Type, Expected>);
    EXPECT_TRUE(
      (std::is_same_v<FormLanguage::FunctionDerived<Expected>::Type, Expected>));
  }

  /// @brief Retains the exact abstract base parameter of a complex expression.
  TEST(FormLanguage_Traits, ComplexAbstractFunctionParameter)
  {
    using Expression = Variational::ComplexFunction<Complex>;
    using Expected = Variational::ScalarFunctionBase<Complex,
      Variational::ComplexFunctionBase<Expression>>;
    static_assert(std::is_abstract_v<Expected>);
    static_assert(
      std::is_same_v<FormLanguage::FunctionDerived<Expression>::Type, Expected>);
    EXPECT_TRUE(
      (std::is_same_v<FormLanguage::FunctionDerived<const Expression&>::Type, Expected>));
  }
}
