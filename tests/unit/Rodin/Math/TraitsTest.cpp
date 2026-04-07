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
  static_assert(std::is_same_v<typename Sum<Real, Real>::Type, Real>);
  static_assert(std::is_same_v<typename Sum<Integer, Integer>::Type, Integer>);
  static_assert(std::is_same_v<typename Sum<Complex, Complex>::Type, Complex>);
  static_assert(std::is_same_v<typename Sum<Real, Integer>::Type, Real>);
}

// Test Minus trait
TEST_F(TraitsTest, MinusTrait)
{
  static_assert(std::is_same_v<typename Minus<Real, Real>::Type, Real>);
  static_assert(std::is_same_v<typename Minus<Integer, Integer>::Type, Integer>);
  static_assert(std::is_same_v<typename Minus<Complex, Complex>::Type, Complex>);
  static_assert(std::is_same_v<typename Minus<Real, Integer>::Type, Real>);

  static_assert(std::is_same_v<typename UnaryMinus<Real>::Type, Real>);
  static_assert(std::is_same_v<typename UnaryMinus<Integer>::Type, Integer>);
  static_assert(std::is_same_v<typename UnaryMinus<Complex>::Type, Complex>);
}

// Test Mult trait
TEST_F(TraitsTest, MultTrait)
{
  static_assert(std::is_same_v<typename Mult<Real, Real>::Type, Real>);
  static_assert(std::is_same_v<typename Mult<Integer, Integer>::Type, Integer>);
  static_assert(std::is_same_v<typename Mult<Complex, Complex>::Type, Complex>);
  static_assert(std::is_same_v<typename Mult<Real, Integer>::Type, Real>);
  static_assert(std::is_same_v<typename Mult<Complex, Real>::Type, Complex>);
}

// Test Division trait
TEST_F(TraitsTest, DivisionTrait)
{
  static_assert(std::is_same_v<typename Division<Real, Real>::Type, Real>);
  static_assert(std::is_same_v<typename Division<Integer, Integer>::Type, Integer>);
  static_assert(std::is_same_v<typename Division<Complex, Complex>::Type, Complex>);
  static_assert(std::is_same_v<typename Division<Real, Integer>::Type, Real>);
  static_assert(std::is_same_v<typename Division<Complex, Real>::Type, Complex>);
}

// Test Dot trait
TEST_F(TraitsTest, DotTrait)
{
  static_assert(std::is_same_v<typename Dot<Real, Real>::Type, Real>);
  static_assert(std::is_same_v<typename Dot<Integer, Integer>::Type, Real>);
  static_assert(std::is_same_v<typename Dot<Complex, Complex>::Type, Complex>);
  static_assert(std::is_same_v<typename Dot<Real, Integer>::Type, Real>);
}

// Test trait composition and nesting
TEST_F(TraitsTest, TraitComposition)
{
  static_assert(std::is_same_v<typename Sum<typename Sum<Real, Real>::Type, Real>::Type, Real>);
  static_assert(std::is_same_v<typename Mult<typename Mult<Real, Real>::Type, Real>::Type, Real>);
  static_assert(std::is_same_v<
    typename Sum<typename Mult<Real, Real>::Type, typename Division<Real, Real>::Type>::Type,
    Real>);
}

// Test traits with vector and matrix types
TEST_F(TraitsTest, VectorMatrixTraits)
{
  using VectorSumType = typename Sum<Math::Vector<Real>, Math::Vector<Real>>::Type;
  using MatrixVectorMultType = typename Mult<Math::Matrix<Real>, Math::Vector<Real>>::Type;

  static_assert(IsEigenObject<VectorSumType>::Value);
  static_assert(IsEigenObject<MatrixVectorMultType>::Value);
  static_assert(std::is_same_v<
    typename Dot<Math::Vector<Real>, Math::Vector<Real>>::Type,
    Real>);
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
