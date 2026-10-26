/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Spatial algebra interoperability with Eigen vector and matrix expressions. */
#include <gtest/gtest.h>
#include <sstream>
#include <type_traits>
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialVector.h"

using namespace Rodin;
using namespace Rodin::Math;

namespace
{
  template <class Scalar>
  class SpatialEigenInterop : public ::testing::Test
  {};
  using Scalars = ::testing::Types<Float, Real, Complex>;
  TYPED_TEST_SUITE(SpatialEigenInterop, Scalars);

  /// @brief Both vector operand orders preserve active coordinates and column conversion.
  TYPED_TEST(SpatialEigenInterop, VectorSumsDifferencesAndColumnConversion)
  {
    using Scalar = TypeParam;
    for (std::uint8_t size = 0; size <= 3; ++size)
    {
      SpatialVector<Scalar> spatial(size);
      Vector<Scalar> eigen(size);
      for (size_t i = 0; i < size; ++i)
      {
        spatial[i] = Scalar(i + 1);
        eigen(i) = Scalar(5 + 2 * i);
      }
      const auto sumLeft = eigen + spatial;
      const auto sumRight = spatial + eigen;
      const auto differenceLeft = eigen - spatial;
      const auto differenceRight = spatial - eigen;
      const SpatialMatrix<Scalar> column(spatial);
      EXPECT_EQ(column.rows(), size);
      EXPECT_EQ(column.cols(), 1u);
      EXPECT_EQ(sumLeft.size(), size);
      EXPECT_EQ(sumRight.size(), size);
      EXPECT_EQ(differenceLeft.size(), size);
      EXPECT_EQ(differenceRight.size(), size);
      for (size_t i = 0; i < size; ++i)
      {
        EXPECT_EQ(sumLeft[i], Scalar(6 + 3 * i));
        EXPECT_EQ(sumRight[i], Scalar(6 + 3 * i));
        EXPECT_EQ(differenceLeft[i], Scalar(4 + i));
        EXPECT_EQ(differenceRight[i], -Scalar(4 + i));
        EXPECT_EQ(column(i, 0), Scalar(i + 1));
        EXPECT_EQ(spatial[i], Scalar(i + 1));
        EXPECT_EQ(eigen(i), Scalar(5 + 2 * i));
      }
      for (size_t i = size; i < 3; ++i)
        EXPECT_EQ(column.getData()(i, 0), Scalar(0));
    }
  }

  /// @brief Rectangular Eigen expressions preserve operand order, shape and scalar products.
  TYPED_TEST(SpatialEigenInterop, MatrixExpressionsAndVectorProducts)
  {
    using Scalar = TypeParam;
    for (std::uint8_t rows = 0; rows <= 3; ++rows)
    {
      for (std::uint8_t columns = 0; columns <= 3; ++columns)
      {
        SpatialMatrix<Scalar> spatial(rows, columns);
        Matrix<Scalar> eigen(rows, columns);
        SpatialVector<Scalar> vector(rows);
        for (size_t row = 0; row < rows; ++row)
        {
          vector[row] = Scalar(row + 1);
          for (size_t column = 0; column < columns; ++column)
          {
            spatial(row, column) = Scalar(1 + row + column);
            eigen(row, column) = Scalar(4 + 2 * row + column);
          }
        }
        const auto sum = spatial + eigen;
        const auto reversedSum = eigen + spatial;
        const auto difference = spatial - eigen;
        const auto reversedDifference = eigen - spatial;
        const auto product = vector * eigen;
        EXPECT_EQ(product.rows(), 1u);
        EXPECT_EQ(product.cols(), columns);
        const auto checkShape = [&](const auto& matrix) {
          EXPECT_EQ(matrix.rows(), rows);
          EXPECT_EQ(matrix.cols(), columns);
        };
        checkShape(sum);
        checkShape(reversedSum);
        checkShape(difference);
        checkShape(reversedDifference);
        for (size_t row = 0; row < rows; ++row)
        {
          for (size_t column = 0; column < columns; ++column)
          {
            EXPECT_EQ(sum(row, column), Scalar(5 + 3 * row + 2 * column));
            EXPECT_EQ(reversedSum(row, column), Scalar(5 + 3 * row + 2 * column));
            EXPECT_EQ(difference(row, column), -Scalar(3 + row));
            EXPECT_EQ(reversedDifference(row, column), Scalar(3 + row));
            EXPECT_EQ(spatial(row, column), Scalar(1 + row + column));
          }
        }
        for (size_t column = 0; column < columns; ++column)
        {
          Scalar expected = 0;
          for (size_t row = 0; row < rows; ++row)
            expected += Scalar((row + 1) * (4 + 2 * row + column));
          EXPECT_EQ(product(0, column), expected);
        }
        std::ostringstream actual, expected;
        EXPECT_EQ(&(actual << spatial), &actual);
        expected << spatial.getData().topLeftCorner(rows, columns);
        EXPECT_EQ(actual.str(), expected.str());
      }
    }
  }

  /// @brief Eigen-vector cross products preserve orientation and complex bilinearity.
  TYPED_TEST(SpatialEigenInterop, CrossProductWithEigenExpression)
  {
    using Scalar = TypeParam;
    SpatialVector<Scalar> a{Scalar(1), Scalar(2), Scalar(3)};
    Eigen::Matrix<Scalar, 3, 1> b;
    b << Scalar(3), Scalar(2), Scalar(1);
    if constexpr (std::is_same_v<Scalar, Complex>)
    {
      a[1] += Complex(0, 1);
      b(2) -= Complex(0, 1);
    }
    const auto product = a.cross(b);
    EXPECT_EQ(product.size(), 3u);
    if constexpr (std::is_same_v<Scalar, Complex>)
    {
      EXPECT_EQ(product[0], Complex(-3, -1));
      EXPECT_EQ(product[1], Complex(8, 1));
      EXPECT_EQ(product[2], Complex(-4, -3));
    }
    else
    {
      EXPECT_EQ(product[0], Scalar(-4));
      EXPECT_EQ(product[1], Scalar(8));
      EXPECT_EQ(product[2], Scalar(-4));
    }
    EXPECT_EQ(a[0], Scalar(1));
    EXPECT_EQ(b(0), Scalar(3));
  }
}
