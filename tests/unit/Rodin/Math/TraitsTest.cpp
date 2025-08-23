/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <type_traits>

#include "Rodin/Math/Traits.h"
#include "Rodin/Types.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/Matrix.h"

#include <Eigen/Core>

using namespace Rodin;
using namespace Rodin::FormLanguage;

class TraitsTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// Test IsEigenObject trait
TEST_F(TraitsTest, IsEigenObjectTrait)
{
  // Test with Eigen types
  static_assert(IsEigenObject<Eigen::VectorXd>::Value);
  static_assert(IsEigenObject<Eigen::MatrixXd>::Value);
  static_assert(IsEigenObject<Eigen::Matrix3d>::Value);
  static_assert(IsEigenObject<Eigen::Vector3d>::Value);

  // Test with Rodin vector and matrix types
  static_assert(IsEigenObject<Math::Vector<Real>>::Value);
  static_assert(IsEigenObject<Math::Matrix<Real>>::Value);
  static_assert(IsEigenObject<Rodin::Math::RealMatrix>::Value);
  static_assert(IsEigenObject<Rodin::Math::ComplexMatrix>::Value);

  // Test with non-Eigen types
  static_assert(!IsEigenObject<int>::Value);
  static_assert(!IsEigenObject<Real>::Value);
  static_assert(!IsEigenObject<Complex>::Value);
  static_assert(!IsEigenObject<std::vector<Real>>::Value);
}

// Test IsEigenObject with const and reference types
TEST_F(TraitsTest, IsEigenObjectWithCVRef)
{
  // Test with const types
  static_assert(IsEigenObject<const Eigen::VectorXd>::Value);
  static_assert(IsEigenObject<const Math::Vector<Real>>::Value);

  // Test with reference types
  static_assert(IsEigenObject<Eigen::VectorXd&>::Value);
  static_assert(IsEigenObject<const Eigen::VectorXd&>::Value);
  static_assert(IsEigenObject<Math::Vector<Real>&>::Value);

  // Test with rvalue reference types
  static_assert(IsEigenObject<Eigen::VectorXd&&>::Value);
  static_assert(IsEigenObject<Math::Vector<Real>&&>::Value);
}

// Test Traits specializations for basic types
TEST_F(TraitsTest, BasicTypeTraits)
{
  // Test Boolean traits
  using BooleanTraits = Traits<Boolean>;
  static_assert(std::is_same_v<BooleanTraits::ScalarType, Boolean>);

  // Test Integer traits
  using IntegerTraits = Traits<Integer>;
  static_assert(std::is_same_v<IntegerTraits::ScalarType, Integer>);

  // Test Real traits
  using RealTraits = Traits<Real>;
  static_assert(std::is_same_v<RealTraits::ScalarType, Real>);

  // Test Complex traits
  using ComplexTraits = Traits<Complex>;
  static_assert(std::is_same_v<ComplexTraits::ScalarType, Complex>);
}

// Test Sum trait
TEST_F(TraitsTest, SumTrait)
{
  // Note: These tests check that the trait compiles and has the expected behavior
  // The actual Math::sum function would need to be implemented for full testing

  // Test that Sum trait exists and can be instantiated for basic types
  using RealSum = Sum<Real, Real>;
  using IntSum = Sum<Integer, Integer>;
  using ComplexSum = Sum<Complex, Complex>;

  // Test mixed type sums
  using MixedSum = Sum<Real, Integer>;

  // These tests verify the trait exists and compiles
  static_assert(std::is_class_v<RealSum>);
  static_assert(std::is_class_v<IntSum>);
  static_assert(std::is_class_v<ComplexSum>);
  static_assert(std::is_class_v<MixedSum>);
}

// Test Minus trait
TEST_F(TraitsTest, MinusTrait)
{
  // Test binary minus trait
  using RealMinus = Minus<Real, Real>;
  using IntMinus = Minus<Integer, Integer>;
  using ComplexMinus = Minus<Complex, Complex>;

  // Test mixed type minus
  using MixedMinus = Minus<Real, Integer>;

  static_assert(std::is_class_v<RealMinus>);
  static_assert(std::is_class_v<IntMinus>);
  static_assert(std::is_class_v<ComplexMinus>);
  static_assert(std::is_class_v<MixedMinus>);

  // Test unary minus trait
  using RealUnaryMinus = UnaryMinus<Real>;
  using IntUnaryMinus = UnaryMinus<Integer>;
  using ComplexUnaryMinus = UnaryMinus<Complex>;

  static_assert(std::is_class_v<RealUnaryMinus>);
  static_assert(std::is_class_v<IntUnaryMinus>);
  static_assert(std::is_class_v<ComplexUnaryMinus>);
}

// Test Mult trait
TEST_F(TraitsTest, MultTrait)
{
  using RealMult = Mult<Real, Real>;
  using IntMult = Mult<Integer, Integer>;
  using ComplexMult = Mult<Complex, Complex>;

  // Test mixed type multiplication
  using MixedMult = Mult<Real, Integer>;
  using ComplexRealMult = Mult<Complex, Real>;

  static_assert(std::is_class_v<RealMult>);
  static_assert(std::is_class_v<IntMult>);
  static_assert(std::is_class_v<ComplexMult>);
  static_assert(std::is_class_v<MixedMult>);
  static_assert(std::is_class_v<ComplexRealMult>);
}

// Test Division trait
TEST_F(TraitsTest, DivisionTrait)
{
  using RealDiv = Division<Real, Real>;
  using IntDiv = Division<Integer, Integer>;
  using ComplexDiv = Division<Complex, Complex>;

  // Test mixed type division
  using MixedDiv = Division<Real, Integer>;
  using ComplexRealDiv = Division<Complex, Real>;

  static_assert(std::is_class_v<RealDiv>);
  static_assert(std::is_class_v<IntDiv>);
  static_assert(std::is_class_v<ComplexDiv>);
  static_assert(std::is_class_v<MixedDiv>);
  static_assert(std::is_class_v<ComplexRealDiv>);
}

// Test Dot trait
TEST_F(TraitsTest, DotTrait)
{
  using RealDot = Dot<Real, Real>;
  using IntDot = Dot<Integer, Integer>;
  using ComplexDot = Dot<Complex, Complex>;

  // Test mixed type dot product
  using MixedDot = Dot<Real, Integer>;
  using ComplexRealDot = Dot<Complex, Real>;

  static_assert(std::is_class_v<RealDot>);
  static_assert(std::is_class_v<IntDot>);
  static_assert(std::is_class_v<ComplexDot>);
  static_assert(std::is_class_v<MixedDot>);
  static_assert(std::is_class_v<ComplexRealDot>);
}

// Test trait composition and nesting
TEST_F(TraitsTest, TraitComposition)
{
  // Test that traits can be composed
  using NestedSum = Sum<Sum<Real, Real>, Real>;
  using NestedMult = Mult<Mult<Real, Real>, Real>;

  static_assert(std::is_class_v<NestedSum>);
  static_assert(std::is_class_v<NestedMult>);

  // Test complex expressions
  using ComplexExpr = Sum<Mult<Real, Real>, Division<Real, Real>>;
  static_assert(std::is_class_v<ComplexExpr>);
}

// Test traits with vector and matrix types
TEST_F(TraitsTest, VectorMatrixTraits)
{
  // Test Sum with vector types
  using VectorSum = Sum<Math::Vector<Real>, Math::Vector<Real>>;
  static_assert(std::is_class_v<VectorSum>);

  // Test Mult with matrix and vector
  using MatrixVectorMult = Mult<Math::Matrix<Real>, Math::Vector<Real>>;
  static_assert(std::is_class_v<MatrixVectorMult>);

  // Test Dot with vectors
  using VectorDot = Dot<Math::Vector<Real>, Math::Vector<Real>>;
  static_assert(std::is_class_v<VectorDot>);
}

// Test trait edge cases
TEST_F(TraitsTest, EdgeCases)
{
  // Test traits with the same type
  using SelfSum = Sum<Real, Real>;
  using SelfMult = Mult<Complex, Complex>;

  static_assert(std::is_class_v<SelfSum>);
  static_assert(std::is_class_v<SelfMult>);

  // Test traits with cv-qualified types
  using ConstSum = Sum<const Real, Real>;
  using VolatileSum = Sum<volatile Real, Real>;

  static_assert(std::is_class_v<ConstSum>);
  static_assert(std::is_class_v<VolatileSum>);
}

// Test that all traits are properly defined as structs
TEST_F(TraitsTest, TraitStructureValidation)
{
  // Verify that all traits are properly structured
  static_assert(std::is_class_v<IsEigenObject<Real>>);
  static_assert(std::is_class_v<Traits<Real>>);
  static_assert(std::is_class_v<Sum<Real, Real>>);
  static_assert(std::is_class_v<Minus<Real, Real>>);
  static_assert(std::is_class_v<UnaryMinus<Real>>);
  static_assert(std::is_class_v<Mult<Real, Real>>);
  static_assert(std::is_class_v<Division<Real, Real>>);
  static_assert(std::is_class_v<Dot<Real, Real>>);

  // Test that they have the required static members
  static_assert(IsEigenObject<Eigen::VectorXd>::Value == true);
  static_assert(IsEigenObject<Real>::Value == false);
}
