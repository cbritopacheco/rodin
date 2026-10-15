/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Empty spatial algebra and boundary-size storage contracts. */
#include <gtest/gtest.h>
#include <sstream>
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialVector.h"

using namespace Rodin;
using namespace Rodin::Math;

namespace
{
  template <class Scalar>
  class SpatialBoundary : public ::testing::Test
  {};
  using SpatialScalars = ::testing::Types<Float, Real, Complex>;
  TYPED_TEST_SUITE(SpatialBoundary, SpatialScalars);

  /// @brief Binary-exact scaling, negation and array assignment preserve every active size.
  TYPED_TEST(SpatialBoundary, ScalingAndArrayAssignment)
  {
    for (std::uint8_t size = 0; size <= 3; ++size)
    {
      SCOPED_TRACE(unsigned(size));
      Math::Vector<TypeParam> entries(size);
      for (size_t i = 0; i < size; ++i)
      {
        entries[i] = TypeParam(i + 1);
        if constexpr (std::is_same_v<TypeParam, Complex>)
          entries[i] += Complex(0, Real(i + 2));
      }
      SpatialVector<TypeParam> values;
      values = entries.array();
      ASSERT_EQ(values.size(), size);
      const auto negative = -values;
      const auto divided = values / TypeParam(2);
      EXPECT_EQ(negative.size(), size);
      EXPECT_EQ(divided.size(), size);
      EXPECT_EQ(&(values *= TypeParam(2)), &values);
      for (size_t i = 0; i < size; ++i)
      {
        EXPECT_EQ(values[i], TypeParam(2) * entries[i]);
        EXPECT_EQ(negative[i], -entries[i]);
        EXPECT_EQ(divided[i], entries[i] / TypeParam(2));
      }
      EXPECT_EQ(&(values /= TypeParam(2)), &values);
      for (size_t i = 0; i < size; ++i)
        EXPECT_EQ(values[i], entries[i]);
    }
  }

  /// @brief Rectangular identity maps have exact Kronecker entries, including empty axes.
  TYPED_TEST(SpatialBoundary, RectangularIdentity)
  {
    for (std::uint8_t rows = 0; rows <= 3; ++rows)
    {
      for (std::uint8_t columns = 0; columns <= 3; ++columns)
      {
        const auto identity = SpatialMatrix<TypeParam>::Identity(rows, columns);
        EXPECT_EQ(identity.rows(), rows);
        EXPECT_EQ(identity.cols(), columns);
        for (size_t i = 0; i < rows; ++i)
        {
          for (size_t j = 0; j < columns; ++j)
            EXPECT_EQ(identity(i, j), TypeParam(i == j));
        }
      }
    }
  }

  /// @brief Zero-dimensional vectors preserve empty active storage and stream no entries.
  TYPED_TEST(SpatialBoundary, EmptyVectorsAndStreaming)
  {
    using V = SpatialVector<TypeParam>;
    const Math::Vector<TypeParam> eigenEmpty(0);
    V empty(eigenEmpty), assigned(3);
    assigned = eigenEmpty;
    EXPECT_EQ(empty.size(), 0u);
    EXPECT_EQ(assigned.size(), 0u);
    EXPECT_EQ(V::Zero(0).size(), 0u);
    EXPECT_EQ(empty.squaredNorm(), 0);
    EXPECT_EQ((empty + assigned).size(), 0u);
    EXPECT_EQ((TypeParam(2) * empty).size(), 0u);
    std::ostringstream emptyStream;
    emptyStream << empty;
    EXPECT_TRUE(emptyStream.str().empty());
    const V zero = V::Zero(3);
    for (size_t i = 0; i < 3; ++i)
      EXPECT_EQ(zero[i], TypeParam(0));
    V values(3);
    for (size_t i = 0; i < 3; ++i)
      values[i] = TypeParam(i + 1);
    std::ostringstream output;
    output << values;
    std::istringstream input(output.str());
    for (size_t i = 0; i < 3; ++i)
    {
      TypeParam value;
      ASSERT_TRUE(bool(input >> value));
      EXPECT_EQ(value, TypeParam(i + 1));
    }
    TypeParam extra;
    EXPECT_FALSE(bool(input >> extra));
  }

  /// @brief Empty determinant, trace and inverse use the conventional empty linear-map identities.
  TYPED_TEST(SpatialBoundary, EmptySquareMatrixIdentities)
  {
    const SpatialMatrix<TypeParam> empty(0, 0);
    EXPECT_EQ(empty.trace(), TypeParam(0));
    EXPECT_EQ(empty.determinant(), TypeParam(1));
    const auto inverse = empty.inverse();
    EXPECT_EQ(inverse.rows(), 0u);
    EXPECT_EQ(inverse.cols(), 0u);
    EXPECT_EQ(inverse.squaredNorm(), 0);
    EXPECT_EQ((empty * inverse).rows(), 0u);
    EXPECT_EQ((empty * inverse).cols(), 0u);
  }

  /// @brief Rectangular empty maps and zero input vectors have the correct output shape and exact zero entries.
  TYPED_TEST(SpatialBoundary, RectangularRowsColumnsAndZeroProducts)
  {
    for (std::uint8_t rows = 0; rows <= 3; ++rows)
    {
      for (std::uint8_t columns = 0; columns <= 3; ++columns)
      {
        SCOPED_TRACE(::testing::Message() << unsigned(rows) << 'x' << unsigned(columns));
        SpatialMatrix<TypeParam> matrix(rows, columns);
        for (size_t row = 0; row < rows; ++row)
        {
          for (size_t column = 0; column < columns; ++column)
            matrix(row, column) = TypeParam(1 + 3 * row + column);
        }
        for (size_t row = 0; row < rows; ++row)
        {
          const auto entries = matrix.row(row);
          EXPECT_EQ(entries.size(), columns);
          for (size_t column = 0; column < columns; ++column)
            EXPECT_EQ(entries[column], TypeParam(1 + 3 * row + column));
        }
        for (size_t column = 0; column < columns; ++column)
        {
          const auto entries = matrix.col(column);
          EXPECT_EQ(entries.size(), rows);
          for (size_t row = 0; row < rows; ++row)
            EXPECT_EQ(entries[row], TypeParam(1 + 3 * row + column));
        }
        const auto result = matrix * SpatialVector<TypeParam>::Zero(columns);
        EXPECT_EQ(result.size(), rows);
        for (size_t row = 0; row < rows; ++row)
          EXPECT_EQ(result[row], TypeParam(0));
      }
    }
  }
}
