/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Rodin vector inner products, ownership and Eigen interoperability.
 */
#include <gtest/gtest.h>
#include "Rodin/Math/Traits.h"
#include "Rodin/Math/Vector.h"

using namespace Rodin;

namespace Rodin::Tests::Unit
{
  /// @brief Exercises the same algebra for dynamic and fixed Eigen-backed owners.
  template <class Vector>
  class VectorContractTest : public ::testing::Test
  {};

  using VectorTypes = ::testing::Types<Math::Vector<float>, Math::Vector<Real>,
    Math::Vector<Complex>, Math::Vector<std::complex<float>>,
    Math::Vector<std::complex<long double>>, Math::FixedSizeVector<Real, 3>,
    Math::FixedSizeVector<Complex, 3>, Math::FixedSizeVector<std::complex<float>, 3>,
    Math::FixedSizeVector<std::complex<long double>, 3>>;
  TYPED_TEST_SUITE(VectorContractTest, VectorTypes);

  /// @brief Member, expression and free contractions conjugate the right operand.
  TYPED_TEST(VectorContractTest, RightConjugatingProducts)
  {
    using Vector = TypeParam;
    using Scalar = typename Vector::Scalar;
    Vector left(3), right(3);
    for (int i = 0; i < 3; ++i)
    {
      if constexpr (Eigen::NumTraits<Scalar>::IsComplex)
      {
        left(i) = Scalar(i + 1, 2 - i);
        right(i) = Scalar(2 - i, -i - 3);
      }
      else
      {
        left(i) = Scalar(i + 1);
        right(i) = Scalar(2 - i);
      }
    }
    Scalar expected = 0;
    for (int i = 0; i < 3; ++i)
    {
      if constexpr (Eigen::NumTraits<Scalar>::IsComplex)
        expected += left(i) * std::conj(right(i));
      else
        expected += left(i) * right(i);
    }
    EXPECT_EQ(left.dot(right), expected);
    EXPECT_EQ(left.dot(right.head(3)), expected);
    EXPECT_EQ(left.dot(2 * right), Scalar(2) * expected);
    EXPECT_EQ(Math::dot(left, right), expected);
    EXPECT_EQ(Math::dot(left.head(3), right.head(3)), expected);
    EXPECT_EQ(Math::dot(left + left, right), Scalar(2) * expected);
    if constexpr (Eigen::NumTraits<Scalar>::IsComplex)
    {
      EXPECT_NE(expected, Scalar(std::conj(expected)));
      EXPECT_EQ(left.head(3).dot(right.head(3)), std::conj(expected));
    }
  }

  /// @brief Owners retain Eigen layout, traits, expression construction and solver use.
  TYPED_TEST(VectorContractTest, EigenExpressionsAndOwnership)
  {
    using Vector = TypeParam;
    using Scalar = typename Vector::Scalar;
    using Parent = typename Vector::Parent;
    static_assert(sizeof(Vector) == sizeof(Parent));
    static_assert(alignof(Vector) == alignof(Parent));
    static_assert(FormLanguage::IsEigenObject<const Vector&>::Value);
    static_assert(FormLanguage::IsVectorRange<Vector>::Value);
    static_assert(FormLanguage::IsEigenObject<Parent>::Value);
    static_assert(std::is_same_v<typename FormLanguage::Traits<Vector>::ScalarType,
      typename FormLanguage::Traits<Parent>::ScalarType>);
    Vector left = Vector::Ones(3);
    const Vector right = Vector::Constant(3, Scalar(2));
    Vector sum = left + right;
    EXPECT_EQ(sum(0), Scalar(3));
    sum = 2 * left + right;
    EXPECT_EQ(sum(2), Scalar(4));
    Vector copy(sum);
    copy(0) = Scalar(7);
    EXPECT_EQ(sum(0), Scalar(4));
    Vector moved(std::move(copy));
    EXPECT_EQ(moved(0), Scalar(7));
    left = std::move(moved);
    EXPECT_EQ(left(0), Scalar(7));
    left = sum;
    EXPECT_EQ(left(0), Scalar(4));
    left.head(2) = right.tail(2);
    EXPECT_EQ(left(0), Scalar(2));
    Eigen::Matrix<Scalar, 3, 3> matrix = 2 * Eigen::Matrix<Scalar, 3, 3>::Identity();
    Vector solution = matrix.fullPivLu().solve(right);
    EXPECT_EQ(solution(0), Scalar(1));
    solution.noalias() = matrix * left;
    EXPECT_EQ(solution(2), Scalar(8));
    Eigen::Map<Vector> mapped(left.data(), 3);
    EXPECT_EQ(Math::dot(mapped, right), Math::dot(left, right));
    Eigen::Ref<Vector> reference(left);
    reference(0) = Scalar(5);
    EXPECT_EQ(left(0), Scalar(5));
    const Eigen::Map<const Vector> readOnly(left.data(), 3);
    EXPECT_EQ(readOnly(0), Scalar(5));
    Parent copiedNative = Parent::Constant(3, Scalar(8));
    Vector copiedStorage(copiedNative);
    copiedNative(0) = Scalar(12);
    EXPECT_EQ(copiedStorage(0), Scalar(8));
    Parent native = Parent::Constant(3, Scalar(9));
    const auto* data = native.data();
    Vector converted(std::move(native));
    EXPECT_EQ(converted(0), Scalar(9));
    if constexpr (Vector::SizeAtCompileTime == Eigen::Dynamic)
      EXPECT_EQ(converted.data(), data);
    Parent assigned = Parent::Constant(3, Scalar(11));
    const auto* assignedData = assigned.data();
    converted = std::move(assigned);
    EXPECT_EQ(converted(0), Scalar(11));
    if constexpr (Vector::SizeAtCompileTime == Eigen::Dynamic)
      EXPECT_EQ(converted.data(), assignedData);
  }

  /// @brief Mixed scalar fields use a common precision without losing conjugation.
  TEST(VectorProducts, MixedPrecisionAndEmptyProducts)
  {
    Math::Vector<std::complex<float>> left(1);
    Math::FixedSizeVector<std::complex<long double>, 1> right;
    left(0) = std::complex<float>(1, 2);
    right(0) = std::complex<long double>(2, -3);
    static_assert(std::is_same_v<decltype(left.dot(right)), std::complex<long double>>);
    EXPECT_EQ(left.dot(right), std::complex<long double>(-4, 7));
    EXPECT_EQ(right.dot(left), std::complex<long double>(-4, -7));
    Math::Vector<double> real = Math::Vector<double>::Constant(1, 2);
    EXPECT_EQ(left.dot(real), Complex(2, 4));
    EXPECT_EQ(real.dot(left), Complex(2, -4));
    const Math::Vector<Complex> empty;
    const Math::FixedSizeVector<Complex, 0> fixedEmpty;
    EXPECT_EQ(empty.dot(empty), Complex(0));
    EXPECT_EQ(fixedEmpty.dot(empty), Complex(0));
    EXPECT_EQ(Math::dot(empty, fixedEmpty), Complex(0));
  }
}
