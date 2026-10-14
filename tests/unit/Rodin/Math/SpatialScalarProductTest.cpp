/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Exact real/complex scalar promotion for spatial objects.
 */

#include <limits>
#include <cstring>
#include <random>
#include <type_traits>
#include <vector>
#include <gtest/gtest.h>
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/SpatialMatrix.h"

using namespace Rodin;
using namespace Rodin::Math;

namespace Rodin::Tests::Unit
{
  /**
   * @brief Retains the original spatial dot reduction and ignores inactive storage.
   *
   * Real, complex and both mixed scalar orders are compared bit-for-bit with
   * the original dynamic Eigen reduction in dimensions zero through three.
   * Random non-dyadic values expose rounding changes; cancellation distinguishes
   * left-associated and reassociated sums. NaN padding detects inactive reads.
   */
  TEST(SpatialScalarProductTest, SpatialDotMatchesOriginalReductionExactly)
  {
    const auto check = []<class LHS, class RHS>() {
      SCOPED_TRACE((std::is_same_v<LHS, Complex> ? "Complex" : "Real"));
      SCOPED_TRACE((std::is_same_v<RHS, Complex> ? "Complex" : "Real"));
      std::mt19937 generator(3107);
      std::uniform_real_distribution<Real> sample(-10, 10);
      const auto value = [&]<class Scalar>() {
        if constexpr (std::is_same_v<Scalar, Complex>)
          return Scalar(sample(generator), sample(generator));
        else
          return sample(generator);
      };
      const auto sameBits = [](const auto& actual, const auto& expected) {
        const auto scalar = [](Real a, Real b) {
          if (std::isnan(b))
            EXPECT_TRUE(std::isnan(a));
          else
            EXPECT_EQ(std::memcmp(&a, &b, sizeof(Real)), 0) << a << " != " << b;
        };
        if constexpr (std::is_same_v<std::decay_t<decltype(actual)>, Complex>)
        {
          scalar(actual.real(), expected.real());
          scalar(actual.imag(), expected.imag());
        }
        else
          scalar(actual, expected);
      };
      for (std::uint8_t dim = 0; dim <= 3; ++dim)
      {
        SCOPED_TRACE(static_cast<unsigned>(dim));
        SpatialVector<LHS> lhs(dim);
        SpatialVector<RHS> rhs(dim);
        const auto& left = lhs;
        const auto& right = rhs;
        static_assert(std::is_const_v<std::remove_reference_t<decltype(left.getData())>>);
        static_assert(
          std::is_const_v<std::remove_reference_t<decltype(right.getData())>>);
        lhs.getData().setConstant(LHS(std::numeric_limits<Real>::quiet_NaN()));
        rhs.getData().setConstant(RHS(std::numeric_limits<Real>::quiet_NaN()));
        const auto compare = [&]() {
          const auto n = static_cast<Eigen::Index>(dim);
          const auto original =
            Math::dot(left.getData().head(n), right.getData().head(n));
          sameBits(Math::dot(lhs, rhs), original);
        };
        for (size_t trial = 0; trial < 1000; ++trial)
        {
          for (size_t i = 0; i < dim; ++i)
          {
            lhs(i) = value.template operator()<LHS>();
            rhs(i) = value.template operator()<RHS>();
          }
          compare();
        }
        if (dim == 3)
        {
          lhs(0) = LHS(1e16);
          lhs(1) = LHS(-1e16);
          lhs(2) = LHS(1);
          rhs.getData().setConstant(RHS(1));
          compare();
        }
        for (Real boundary : {Real(0), -Real(0), std::numeric_limits<Real>::infinity()})
        {
          for (size_t i = 0; i < dim; ++i)
          {
            lhs(i) = LHS(boundary);
            rhs(i) = RHS(1);
          }
          compare();
        }
      }
    };
    check.template operator()<Real, Real>();
    check.template operator()<Real, Complex>();
    check.template operator()<Complex, Real>();
    check.template operator()<Complex, Complex>();
  }

  /**
   * @brief Preserves the original complex reduction in pair-loop calling contexts.
   *
   * Non-dyadic operands with nonzero imaginary parts exercise both conjugation
   * and rounding. Growing scratch to three entries and shrinking to one leaves
   * poisoned inactive storage. Every full-batch result must match the original
   * const-view reduction bit-for-bit, and an unconjugated product must differ.
   */
  TEST(SpatialScalarProductTest, SingleComponentComplexBatchesMatchOriginalAfterResize)
  {
    for (const size_t count : {size_t(4), size_t(16), size_t(64)})
    {
      SCOPED_TRACE(count);
      std::vector<SpatialVector<Complex>> lhs(count, SpatialVector<Complex>(3));
      std::vector<SpatialVector<Complex>> rhs(count, SpatialVector<Complex>(3));
      for (size_t b = 0; b < count; ++b)
      {
        lhs[b].resize(1);
        rhs[b].resize(1);
        lhs[b].getData().setConstant(Complex(std::numeric_limits<Real>::quiet_NaN()));
        rhs[b].getData().setConstant(Complex(std::numeric_limits<Real>::quiet_NaN()));
        lhs[b](0) = Complex(Real(b + 1) / 7, Real(b + 2) / 11);
        rhs[b](0) = Complex(Real(b + 2) / 13, -Real(b + 1) / 17);
      }
      const auto& left = lhs;
      const auto& right = rhs;
      for (size_t b = 0; b < count; ++b)
      {
        for (size_t a = 0; a < count; ++a)
        {
          const auto n = static_cast<Eigen::Index>(left[b].size());
          const Complex expected =
            Math::dot(left[b].getData().head(n), right[a].getData().head(n));
          const Complex actual = Math::dot(left[b], right[a]);
          const auto sameBits = [](Real x, Real y) {
            return std::memcmp(&x, &y, sizeof(Real)) == 0;
          };
          ASSERT_TRUE(sameBits(actual.real(), expected.real()));
          ASSERT_TRUE(sameBits(actual.imag(), expected.imag()));
          EXPECT_NE(actual, left[b](0) * right[a](0));
        }
      }
    }
  }

  template <class Coefficient, class Value>
  concept ScalesSpatial = requires(const Coefficient& coefficient, const Value& value) {
    coefficient* value;
    value* coefficient;
  };

  static_assert(!ScalesSpatial<void*, SpatialVector<Real>>);
  static_assert(!ScalesSpatial<void*, SpatialMatrix<Complex>>);

  /// @brief Scalar/matrix promotion preserves dimensions, entries and zero padding.
  TEST(SpatialMatrixRegression, MixedRealComplexScaling)
  {
    static_assert(std::is_same_v<decltype(Real{} * SpatialMatrix<Complex>{}),
      SpatialMatrix<Complex>>);
    static_assert(std::is_same_v<decltype(SpatialMatrix<Complex>{} * Real{}),
      SpatialMatrix<Complex>>);
    static_assert(std::is_same_v<decltype(Complex{} * SpatialMatrix<Real>{}),
      SpatialMatrix<Complex>>);
    static_assert(std::is_same_v<decltype(SpatialMatrix<Real>{} * Complex{}),
      SpatialMatrix<Complex>>);
    static_assert(
      std::is_same_v<decltype(Real{} * SpatialMatrix<Real>{}), SpatialMatrix<Real>>);
    for (std::uint8_t rows = 0; rows <= 3; ++rows)
    {
      for (std::uint8_t cols = 0; cols <= 3; ++cols)
      {
        SpatialMatrix<Real> a(rows, cols);
        SpatialMatrix<Complex> z(rows, cols);
        for (std::uint8_t i = 0; i < rows; ++i)
        {
          for (std::uint8_t j = 0; j < cols; ++j)
          {
            a(i, j) = Real(i + 2 * j + 1);
            z(i, j) = Complex(a(i, j), Real(i + j + 1));
          }
        }
        for (const Real scale : {Real(0), Real(-0.5), Real(2)})
        {
          const Complex c(scale, 0.25);
          const auto left = scale * z;
          const auto right = z * scale;
          const auto promotedLeft = c * a;
          const auto promotedRight = a * c;
          EXPECT_EQ(left.rows(), rows);
          EXPECT_EQ(left.cols(), cols);
          EXPECT_EQ(promotedRight.rows(), rows);
          EXPECT_EQ(promotedRight.cols(), cols);
          for (std::uint8_t i = 0; i < 3; ++i)
          {
            for (std::uint8_t j = 0; j < 3; ++j)
            {
              const bool active = i < rows && j < cols;
              const Complex expected = active ? scale * z(i, j) : Complex(0);
              const Complex promoted = active ? c * a(i, j) : Complex(0);
              EXPECT_EQ(left.getData()(i, j), expected);
              EXPECT_EQ(right.getData()(i, j), expected);
              EXPECT_EQ(promotedLeft.getData()(i, j), promoted);
              EXPECT_EQ(promotedRight.getData()(i, j), promoted);
            }
          }
        }
      }
    }
  }

  template <class Scalar>
  Scalar baseValue()
  {
    if constexpr (std::is_same_v<Scalar, Complex>)
      return Complex(1, 2);
    else
      return Real(1);
  }

  template <class Scalar, class Coefficient>
  void checkVectorProduct(Coefficient coefficient)
  {
    using Factor = std::conditional_t<std::is_same_v<Scalar, Complex> &&
        std::is_arithmetic_v<Coefficient>,
      Real, Coefficient>;
    using Out = typename FormLanguage::Mult<Factor, Scalar>::Type;
    const Factor factor = coefficient;
    for (std::uint8_t dim = 0; dim <= 3; ++dim)
    {
      Math::SpatialVector<Scalar> value(dim);
      for (size_t i = 0; i < dim; ++i)
        value(i) = baseValue<Scalar>() + Real(i);
      const auto left = coefficient * value;
      const auto right = value * coefficient;
      static_assert(
        std::is_same_v<std::decay_t<decltype(left)>, Math::SpatialVector<Out>>);
      static_assert(
        std::is_same_v<std::decay_t<decltype(right)>, Math::SpatialVector<Out>>);
      EXPECT_EQ(left.size(), dim);
      EXPECT_EQ(right.size(), dim);
      for (size_t i = 0; i < dim; ++i)
      {
        EXPECT_EQ(left(i), factor * value(i));
        EXPECT_EQ(right(i), value(i) * factor);
        EXPECT_EQ(value(i), baseValue<Scalar>() + Real(i));
      }
      for (size_t i = dim; i < 3; ++i)
      {
        EXPECT_EQ(left.getData()(i), Out(0));
        EXPECT_EQ(right.getData()(i), Out(0));
      }
    }
  }

  template <class Scalar, class Coefficient>
  void checkMatrixProduct(Coefficient coefficient)
  {
    using Factor = std::conditional_t<std::is_same_v<Scalar, Complex> &&
        std::is_arithmetic_v<Coefficient>,
      Real, Coefficient>;
    using Out = typename FormLanguage::Mult<Factor, Scalar>::Type;
    const Factor factor = coefficient;
    for (std::uint8_t rows = 0; rows <= 3; ++rows)
    {
      for (std::uint8_t cols = 0; cols <= 3; ++cols)
      {
        Math::SpatialMatrix<Scalar> value(rows, cols);
        for (size_t i = 0; i < rows; ++i)
        {
          for (size_t j = 0; j < cols; ++j)
            value(i, j) = baseValue<Scalar>() + Real(i + j);
        }
        const auto left = coefficient * value;
        const auto right = value * coefficient;
        static_assert(
          std::is_same_v<std::decay_t<decltype(left)>, Math::SpatialMatrix<Out>>);
        static_assert(
          std::is_same_v<std::decay_t<decltype(right)>, Math::SpatialMatrix<Out>>);
        EXPECT_EQ(left.rows(), rows);
        EXPECT_EQ(left.cols(), cols);
        EXPECT_EQ(right.rows(), rows);
        EXPECT_EQ(right.cols(), cols);
        for (size_t i = 0; i < rows; ++i)
        {
          for (size_t j = 0; j < cols; ++j)
          {
            EXPECT_EQ(left(i, j), factor * value(i, j));
            EXPECT_EQ(right(i, j), value(i, j) * factor);
            EXPECT_EQ(value(i, j), baseValue<Scalar>() + Real(i + j));
          }
        }
      }
    }
  }

  TEST(SpatialScalarProductTest, VectorProductsPromoteRealAndComplex)
  {
    checkVectorProduct<Real>(Real(0.25));
    checkVectorProduct<Real>(Complex(0.25, 0.5));
    checkVectorProduct<Complex>(Real(0.25));
    checkVectorProduct<Complex>(Complex(0.25, 0.5));
  }

  TEST(SpatialScalarProductTest, MatrixProductsPromoteRealAndComplex)
  {
    checkMatrixProduct<Real>(Real(0.25));
    checkMatrixProduct<Real>(Complex(0.25, 0.5));
    checkMatrixProduct<Complex>(Real(0.25));
    checkMatrixProduct<Complex>(Complex(0.25, 0.5));
  }

  /// @brief Integer/float coefficients must compile without complex promotion loss.
  TEST(SpatialScalarProductTest, ArithmeticCoefficientsScaleComplexContainers)
  {
    checkVectorProduct<Complex>(Integer(2));
    checkVectorProduct<Complex>(0);
    checkVectorProduct<Complex>(-2);
    checkVectorProduct<Complex>(0.5f);
    checkVectorProduct<Complex>(0.5L);
    checkMatrixProduct<Complex>(Integer(2));
    checkMatrixProduct<Complex>(0);
    checkMatrixProduct<Complex>(-2);
    checkMatrixProduct<Complex>(0.5f);
    checkMatrixProduct<Complex>(0.5L);
    checkVectorProduct<Real>(2);
    checkMatrixProduct<Real>(2);
  }

  /// @brief Real scaling must not introduce spurious NaNs through complex zero terms.
  TEST(SpatialScalarProductTest, RealScalingPreservesInfiniteComponent)
  {
    const Real infinity = std::numeric_limits<Real>::infinity();
    Math::SpatialVector<Complex> v(1);
    Math::SpatialMatrix<Complex> a(1, 1);
    v(0) = Complex(1, infinity);
    a(0, 0) = v(0);
    for (const auto value :
      {(Real(2) * v)(0), (v * Real(2))(0), (Real(2) * a)(0, 0), (a * Real(2))(0, 0)})
    {
      EXPECT_EQ(value.real(), Real(2));
      EXPECT_EQ(value.imag(), infinity);
    }
  }

  template <class LHS, class RHS>
  void checkMixedMatrixProduct()
  {
    using Out = typename FormLanguage::Mult<LHS, RHS>::Type;
    for (std::uint8_t rows = 0; rows <= 3; ++rows)
    {
      for (std::uint8_t inner = 0; inner <= 3; ++inner)
      {
        for (std::uint8_t cols = 0; cols <= 3; ++cols)
        {
          Math::SpatialMatrix<LHS> a(rows, inner);
          Math::SpatialMatrix<RHS> b(inner, cols);
          for (size_t i = 0; i < rows; ++i)
          {
            for (size_t k = 0; k < inner; ++k)
              a(i, k) = baseValue<LHS>() + Real(i + k);
          }
          for (size_t k = 0; k < inner; ++k)
          {
            for (size_t j = 0; j < cols; ++j)
              b(k, j) = baseValue<RHS>() + Real(k + j);
          }
          const auto product = a * b;
          static_assert(
            std::is_same_v<std::decay_t<decltype(product)>, Math::SpatialMatrix<Out>>);
          EXPECT_EQ(product.rows(), rows);
          EXPECT_EQ(product.cols(), cols);
          for (size_t i = 0; i < rows; ++i)
          {
            for (size_t j = 0; j < cols; ++j)
            {
              Out expected(0);
              for (size_t k = 0; k < inner; ++k)
                expected += a(i, k) * b(k, j);
              EXPECT_EQ(product(i, j), expected);
            }
          }
        }
      }
    }
  }

  TEST(SpatialScalarProductTest, MatrixMatrixProductsPromoteRealAndComplex)
  {
    checkMixedMatrixProduct<Real, Real>();
    checkMixedMatrixProduct<Real, Complex>();
    checkMixedMatrixProduct<Complex, Real>();
    checkMixedMatrixProduct<Complex, Complex>();
  }
}
