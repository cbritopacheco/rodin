/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Explicit complex algebra contracts for spatial and Eigen counterparts.
 */
#include <limits>
#include <gtest/gtest.h>
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/SpatialTensor.h"

using namespace Rodin;
using namespace Rodin::Math;

namespace Rodin::Tests::Unit
{
  /**
   * @brief Compares independently accumulated real/complex contractions.
   * @tparam Left Left scalar field.
   * @tparam Right Right scalar field.
   */
  template <class Left, class Right>
  void checkSpatialInnerProducts()
  {
    const auto coefficient = []<class Scalar>(size_t index, bool right) {
      const Real real = Real(index + (right ? 2 : 1)) / 4;
      if constexpr (Eigen::NumTraits<Scalar>::IsComplex)
        return Scalar(real, Real(index + 1) * (right ? -0.5 : 0.25));
      else
        return Scalar(real);
    };
    for (std::uint8_t rows = 0; rows <= 3; ++rows)
    {
      SpatialVector<Left> u(rows);
      SpatialVector<Right> v(rows);
      u.getData().setConstant(Left(std::numeric_limits<Real>::quiet_NaN()));
      v.getData().setConstant(Right(std::numeric_limits<Real>::quiet_NaN()));
      using Result = decltype(Math::dot(Left{}, Right{}));
      Result expected(0);
      const auto tolerance =
        64 * std::numeric_limits<typename Eigen::NumTraits<Result>::Real>::epsilon();
      for (size_t i = 0; i < rows; ++i)
      {
        u(i) = coefficient.template operator()<Left>(i, false);
        v(i) = coefficient.template operator()<Right>(i, true);
        expected += Math::dot(u(i), v(i));
      }
      EXPECT_NEAR(std::abs(Math::dot(u, v) - expected), 0, tolerance);
      EXPECT_NEAR(
        std::abs(Math::dot(u, v.getData().head(rows)) - expected), 0, tolerance);
      EXPECT_NEAR(
        std::abs(Math::dot(u.getData().head(rows), v) - expected), 0, tolerance);
      EXPECT_NEAR(std::abs(u.dot(v) - expected), 0, tolerance);
      EXPECT_NEAR(std::abs(u.dot(v.getData().head(rows)) - expected), 0, tolerance);
      if constexpr (std::is_same_v<Left, Right>)
      {
        EXPECT_NEAR(
          std::abs(v.getData().head(rows).dot(u.getData().head(rows)) - expected), 0,
          tolerance);
      }
      for (std::uint8_t cols = 0; cols <= 3; ++cols)
      {
        SpatialMatrix<Left> a(rows, cols);
        SpatialMatrix<Right> b(rows, cols);
        a.getData().setConstant(Left(std::numeric_limits<Real>::quiet_NaN()));
        b.getData().setConstant(Right(std::numeric_limits<Real>::quiet_NaN()));
        expected = 0;
        for (size_t i = 0; i < rows; ++i)
        {
          for (size_t j = 0; j < cols; ++j)
          {
            a(i, j) = coefficient.template operator()<Left>(3 * i + j, false);
            b(i, j) = coefficient.template operator()<Right>(i + 3 * j, true);
            expected += Math::dot(a(i, j), b(i, j));
          }
        }
        EXPECT_NEAR(std::abs(Math::dot(a, b) - expected), 0, tolerance);
        EXPECT_NEAR(
          std::abs(Math::dot(a, b.getData().topLeftCorner(rows, cols)) - expected), 0,
          tolerance);
        EXPECT_NEAR(
          std::abs(Math::dot(a.getData().topLeftCorner(rows, cols), b) - expected), 0,
          tolerance);
        EXPECT_NEAR(std::abs(a.dot(b) - expected), 0, tolerance);
        EXPECT_NEAR(std::abs(a.dot(b.getData().topLeftCorner(rows, cols)) - expected), 0,
          tolerance);
      }
    }
  }

  TEST(SpatialComplexTest, InnerProductsUseTheirDeclaredConjugation)
  {
    checkSpatialInnerProducts<Real, Real>();
    checkSpatialInnerProducts<Real, Complex>();
    checkSpatialInnerProducts<Complex, Real>();
    checkSpatialInnerProducts<Complex, Complex>();
    const SpatialVector<Complex> u{Complex(1, 2), Complex(3, -1)};
    const SpatialVector<Complex> v{Complex(2, -3), Complex(-1, 4)};
    // A nonsymmetric complex pair distinguishes all three conventions exactly.
    EXPECT_EQ(u.dot(v), Math::dot(u, v));
    EXPECT_NE(Math::dot(u, v), u.getData().head(2).dot(v.getData().head(2)));
    EXPECT_EQ(Math::dot(u, v), std::conj(Math::dot(v, u)));
    EXPECT_EQ(Math::dot(Complex(1, 2), Complex(2, -3)), Complex(-4, 7));
  }

  /// @brief Complex contractions support each floating precision and mixed fields.
  TEST(SpatialComplexTest, InnerProductsPromoteScalarPrecision)
  {
    checkSpatialInnerProducts<std::complex<float>, std::complex<float>>();
    checkSpatialInnerProducts<std::complex<long double>, std::complex<long double>>();
    checkSpatialInnerProducts<float, std::complex<float>>();
    checkSpatialInnerProducts<std::complex<float>, double>();
    checkSpatialInnerProducts<double, std::complex<float>>();
    checkSpatialInnerProducts<std::complex<float>, std::complex<long double>>();
  }

  /// @brief Scalar conjugation and tensor contractions preserve extended precision.
  TEST(SpatialComplexTest, ConjugationPreservesScalarPrecision)
  {
    const long double value = 1 + 4 * std::numeric_limits<long double>::epsilon();
    const std::complex<long double> a(value, value), b(2, -3);
    static_assert(std::is_same_v<decltype(Math::conj(a)), std::complex<long double>>);
    static_assert(std::is_same_v<decltype(Math::conj(value)), long double>);
    static_assert(std::is_same_v<decltype(Math::dot(a, std::complex<float>())),
      std::complex<long double>>);
    EXPECT_EQ(Math::conj(a), std::conj(a));
    EXPECT_EQ(Math::conj(value), value);
    EXPECT_EQ(Math::dot(a, b), a * std::conj(b));
    EXPECT_EQ(Math::dot(value, b), value * std::conj(b));
    EXPECT_EQ(Math::dot(a, 2.0f), a * 2.0L);
    SpatialMatrix<std::complex<long double>> matrix(1, 1);
    matrix(0, 0) = a;
    EXPECT_EQ(matrix.conjugate()(0, 0), std::conj(a));
    EXPECT_EQ(matrix.adjoint()(0, 0), std::conj(a));
    SpatialTensor<std::complex<long double>> left({1, 1, 1});
    SpatialTensor<std::complex<float>> right({1, 1, 1});
    left[0] = a;
    right[0] = std::complex<float>(2, -3);
    EXPECT_EQ(left.conjugate()[0], std::conj(a));
    EXPECT_EQ(left.dot(right), a * std::conj(b));
    EXPECT_EQ(Math::dot(left, right), a * std::conj(b));
  }

  /// @brief The optimized second norm ignores padding and agrees with the norm API.
  TEST(SpatialComplexTest, SecondNormUsesActiveEuclideanNorm)
  {
    const auto check = []<class Scalar>() {
      using RealScalar = typename Eigen::NumTraits<Scalar>::Real;
      for (std::uint8_t size = 0; size <= 3; ++size)
      {
        SpatialVector<Scalar> value(size);
        value.getData().setConstant(Scalar(std::numeric_limits<RealScalar>::quiet_NaN()));
        for (size_t i = 0; i < size; ++i)
        {
          if constexpr (Eigen::NumTraits<Scalar>::IsComplex)
            value(i) = Scalar(i + 3, i + 4);
          else
            value(i) = Scalar(i + 3);
        }
        EXPECT_EQ(value.template lpNorm<2>(), value.norm());
        EXPECT_NEAR(value.template lpNorm<1>(),
          value.getData().head(size).template lpNorm<1>(),
          32 * std::numeric_limits<RealScalar>::epsilon());
        EXPECT_EQ(value.template lpNorm<2>(), value.getData().head(size).norm());
      }
    };
    check.template operator()<float>();
    check.template operator()<double>();
    check.template operator()<long double>();
    check.template operator()<std::complex<float>>();
    check.template operator()<Complex>();
    check.template operator()<std::complex<long double>>();
  }

  TEST(SpatialComplexTest, ConjugateTransposeAndAdjointExcludePadding)
  {
    for (std::uint8_t rows = 0; rows <= 3; ++rows)
    {
      SpatialVector<Complex> v(rows);
      v.getData().setConstant(Complex(std::numeric_limits<Real>::quiet_NaN()));
      for (size_t i = 0; i < rows; ++i)
        v(i) = Complex(Real(i + 1), -Real(i + 2));
      const auto vc = v.conjugate(), mc = Math::conj(v);
      const auto vt = v.transpose();
      EXPECT_EQ(vt.rows(), 1);
      EXPECT_EQ(vt.cols(), rows);
      for (size_t i = 0; i < rows; ++i)
      {
        EXPECT_EQ(vc(i), std::conj(v(i)));
        EXPECT_EQ(mc(i), std::conj(v(i)));
        EXPECT_EQ(vt(0, i), v(i));
      }
      for (std::uint8_t cols = 0; cols <= 3; ++cols)
      {
        SpatialMatrix<Complex> a(rows, cols);
        a.getData().setConstant(Complex(std::numeric_limits<Real>::quiet_NaN()));
        for (size_t i = 0; i < rows; ++i)
        {
          for (size_t j = 0; j < cols; ++j)
            a(i, j) = Complex(Real(2 * i + j + 1), -Real(i + 3 * j + 2));
        }
        const auto ac = a.conjugate(), fc = Math::conj(a);
        const auto at = a.transpose(), ah = a.adjoint();
        EXPECT_EQ(ah.rows(), cols);
        EXPECT_EQ(ah.cols(), rows);
        for (size_t i = 0; i < rows; ++i)
        {
          for (size_t j = 0; j < cols; ++j)
          {
            EXPECT_EQ(ac(i, j), std::conj(a(i, j)));
            EXPECT_EQ(fc(i, j), std::conj(a(i, j)));
            EXPECT_EQ(at(j, i), a(i, j));
            EXPECT_EQ(ah(j, i), std::conj(a(i, j)));
          }
        }
        EXPECT_NEAR(std::abs(Math::dot(a, a) - Complex(a.squaredNorm())), 0, 1e-13);
      }
      EXPECT_NEAR(std::abs(Math::dot(v, v) - Complex(v.squaredNorm())), 0, 1e-13);
      if (rows)
        EXPECT_NEAR(v.normalized().norm(), 1, 1e-13);
    }
  }

  TEST(SpatialComplexTest, CrossIsBilinearRatherThanEigenHermitianCross)
  {
    const SpatialVector<Complex> a{Complex(1, 2), Complex(-2, 1), Complex(3, -1)};
    const SpatialVector<Complex> b{Complex(2, -3), Complex(1, 4), Complex(-1, 2)};
    const auto spatial = a.cross(b);
    const auto interop = a.cross(b.getData());
    const auto eigen = a.getData().cross(b.getData()).eval();
    for (size_t i = 0; i < 3; ++i)
    {
      const auto expected =
        a((i + 1) % 3) * b((i + 2) % 3) - a((i + 2) % 3) * b((i + 1) % 3);
      EXPECT_EQ(spatial(i), expected);
      EXPECT_EQ(interop(i), expected);
      EXPECT_EQ(std::conj(eigen(i)), expected);
    }
    EXPECT_NEAR(std::abs(Math::dot(a, spatial.conjugate())), 0, 1e-13);
    EXPECT_GT(std::abs(Math::dot(a, spatial)), 1);
  }

  TEST(SpatialComplexTest, ComplexSolvesAndPseudoInversesUseAdjoints)
  {
    for (std::uint8_t rows = 1; rows <= 3; ++rows)
    {
      for (std::uint8_t cols = 1; cols <= 3; ++cols)
      {
        SpatialMatrix<Complex> a(rows, cols);
        for (size_t i = 0; i < rows; ++i)
        {
          for (size_t j = 0; j < cols; ++j)
          {
            a(i, j) =
              Complex((i == j ? 4 : 0) + Real(i + j + 1) / 8, Real(2 * i + j + 1) / 16);
          }
        }
        const auto p = a.pseudoInverse();
        EXPECT_LT((a * p * a - a).norm(), 1e-12);
        EXPECT_LT((p * a * p - p).norm(), 1e-12);
        const auto ap = a * p, pa = p * a;
        EXPECT_LT((ap - ap.adjoint()).norm(), 1e-12);
        EXPECT_LT((pa - pa.adjoint()).norm(), 1e-12);
        if (rows == cols)
        {
          SpatialVector<Complex> b(rows);
          for (size_t i = 0; i < rows; ++i)
            b(i) = Complex(Real(i + 1), -Real(i + 2));
          EXPECT_LT((a * a.solve(b) - b).norm(), 1e-12);
          EXPECT_LT((a.inverse() * b - a.solve(b)).norm(), 1e-12);
          const auto ref =
            a.getData().topLeftCorner(rows, cols).fullPivLu().determinant();
          EXPECT_NEAR(std::abs(a.determinant() - ref), 0, 1e-12);
        }
      }
    }
  }
}
