/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Exact spatial matrix-vector extent, promotion and padding regressions.
 */
#include <gtest/gtest.h>
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/SpatialMatrix.h"

namespace Rodin::Tests::Unit
{
  /**
   * @brief Checks all legal rectangular shapes with poisoned inactive storage.
   *
   * Dyadic inputs and at most three products make the independent row sums
   * exactly representable. Empty contractions must produce zero active values;
   * empty results and shrinking assignments must leave zero inactive capacity.
   * Complex operands have nonzero imaginary components and multiplication is
   * bilinear, not conjugating. No coordinate or numerical-tolerance predicate
   * is used to determine an extent.
   */
  TEST(Math_SpatialMatrixVector, RuntimeExtentsPromotionAndInactivePadding)
  {
    const auto check = []<class MatrixScalar, class VectorScalar>() {
      using Out = typename FormLanguage::Mult<MatrixScalar, VectorScalar>::Type;
      Math::SpatialVector<Out> result(3);
      for (const size_t rows : {size_t(3), size_t(0), size_t(2), size_t(1)})
      {
        for (const size_t columns : {size_t(3), size_t(0), size_t(1), size_t(2)})
        {
          SCOPED_TRACE(rows);
          SCOPED_TRACE(columns);
          Math::SpatialMatrix<MatrixScalar> matrix(rows, columns);
          Math::SpatialVector<VectorScalar> vector(columns);
          matrix.getData().setConstant(MatrixScalar(128));
          vector.getData().setConstant(VectorScalar(256));
          for (size_t i = 0; i < rows; ++i)
          {
            for (size_t j = 0; j < columns; ++j)
            {
              const Real real = -Real(i + 2 * j + 1) / 8;
              if constexpr (std::is_same_v<MatrixScalar, Complex>)
                matrix(i, j) = MatrixScalar(real, Real(2 * i + j + 1) / 16);
              else
                matrix(i, j) = real;
            }
          }
          for (size_t j = 0; j < columns; ++j)
          {
            const Real real = Real(j + 1) / 4;
            if constexpr (std::is_same_v<VectorScalar, Complex>)
              vector(j) = VectorScalar(real, -Real(j + 2) / 8);
            else
              vector(j) = real;
          }
          static_assert(
            std::is_same_v<decltype(matrix * vector), Math::SpatialVector<Out>>);
          result.getData().setConstant(Out(512));
          result = matrix * vector;
          ASSERT_EQ(result.size(), rows);
          for (size_t i = 0; i < 3; ++i)
          {
            Out expected = 0;
            if (i < rows)
            {
              for (size_t j = 0; j < columns; ++j)
                expected += Out(matrix(i, j)) * Out(vector(j));
            }
            EXPECT_EQ(result.getData()(i), expected);
          }
        }
      }
    };
    check.template operator()<Real, Real>();
    check.template operator()<Real, Complex>();
    check.template operator()<Complex, Real>();
    check.template operator()<Complex, Complex>();
  }
}

/// @brief Scalar factors are captured before writes when they alias an active entry.
TEST(SpatialArithmetic, AliasedScalarFactors)
{
  using namespace Rodin;
  for (size_t dimension = 1; dimension <= 3; ++dimension)
  {
    for (size_t factorIndex = 0; factorIndex < dimension; ++factorIndex)
    {
      Math::SpatialVector<Real> original(dimension);
      for (size_t i = 0; i < dimension; ++i)
        original[i] = Real(i + 2);
      auto multiplied = original;
      multiplied *= multiplied[factorIndex];
      auto divided = original;
      divided /= divided[factorIndex];
      for (size_t i = 0; i < dimension; ++i)
      {
        EXPECT_EQ(multiplied[i], original[i] * original[factorIndex]);
        EXPECT_EQ(divided[i], original[i] / original[factorIndex]);
      }
    }
    for (size_t columns = 1; columns <= 3; ++columns)
    {
      Math::SpatialMatrix<Real> original(dimension, columns);
      for (size_t row = 0; row < dimension; ++row)
      {
        for (size_t col = 0; col < columns; ++col)
          original(row, col) = Real(2 + row + 3 * col);
      }
      for (size_t row = 0; row < dimension; ++row)
      {
        for (size_t col = 0; col < columns; ++col)
        {
          auto scaled = original;
          scaled *= scaled(row, col);
          for (size_t i = 0; i < dimension; ++i)
          {
            for (size_t j = 0; j < columns; ++j)
              EXPECT_EQ(scaled(i, j), original(i, j) * original(row, col));
          }
        }
      }
    }
  }
}

/// @brief A zero inner extent produces active zeros instead of retained padding values.
TEST(SpatialArithmetic, InPlaceZeroInnerProduct)
{
  using namespace Rodin;
  for (size_t rows = 1; rows <= 3; ++rows)
  {
    for (size_t columns = 1; columns <= 3; ++columns)
    {
      Math::SpatialMatrix<Real> lhs(rows, columns);
      for (size_t row = 0; row < rows; ++row)
      {
        for (size_t col = 0; col < columns; ++col)
          lhs(row, col) = Real(5 + row + col);
      }
      lhs.resize(rows, 0);
      Math::SpatialMatrix<Real> rhs(0, columns);
      lhs *= rhs;
      EXPECT_EQ(lhs.rows(), rows);
      EXPECT_EQ(lhs.cols(), columns);
      for (size_t row = 0; row < rows; ++row)
      {
        for (size_t col = 0; col < columns; ++col)
          EXPECT_EQ(lhs(row, col), 0);
      }
    }
  }
}

/// @brief Complex full-rank pseudoinverses satisfy the Hermitian Penrose identities.
TEST(SpatialArithmetic, ComplexPseudoInverse)
{
  using namespace Rodin;
  for (size_t rows = 1; rows <= 3; ++rows)
  {
    for (size_t columns = 1; columns <= 3; ++columns)
    {
      Math::SpatialMatrix<Complex> matrix(rows, columns);
      for (size_t row = 0; row < rows; ++row)
      {
        for (size_t col = 0; col < columns; ++col)
        {
          matrix(row, col) = Complex((row == col ? 2 : 0) + 0.1 * (row + 1) * (col + 1),
            0.2 * (row + 1) - 0.3 * (col + 1));
        }
      }
      const auto inverse = matrix.pseudoInverse();
      const auto left = matrix * inverse;
      const auto right = inverse * matrix;
      EXPECT_LT((left * matrix - matrix).norm(), 1e-11);
      EXPECT_LT((right * inverse - inverse).norm(), 1e-11);
      EXPECT_LT((left - left.adjoint()).norm(), 1e-11);
      EXPECT_LT((right - right.adjoint()).norm(), 1e-11);
    }
  }
}

/// @brief Compound updates evaluate aliased Eigen expressions before writing components.
TEST(SpatialArithmetic, AliasedEigenCompoundUpdates)
{
  using namespace Rodin;
  for (const size_t dimension : {size_t(0), size_t(1), size_t(2), size_t(3)})
  {
    Math::SpatialVector<Real> original(dimension);
    for (size_t i = 0; i < dimension; ++i)
      original[i] = Real(i + 1);
    auto plus = original;
    plus += plus.getData().head(dimension).reverse();
    auto minus = original;
    minus -= minus.getData().head(dimension).reverse();
    for (size_t i = 0; i < dimension; ++i)
    {
      EXPECT_EQ(plus[i], original[i] + original[dimension - 1 - i]);
      EXPECT_EQ(minus[i], original[i] - original[dimension - 1 - i]);
    }
  }
  Math::SpatialVector<Real> product{1, 2, 3};
  const Eigen::Matrix3d matrix = Eigen::Matrix3d::Ones();
  product += matrix * product.getData();
  EXPECT_EQ(product[0], 7);
  EXPECT_EQ(product[1], 8);
  EXPECT_EQ(product[2], 9);
}
