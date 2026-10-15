/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Empty spatial algebra and boundary-size storage contracts. */
#include <gtest/gtest.h>
#include <array>
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

  /// @brief Aliased matrix expressions preserve values, rectangular shapes and zeroed inactive storage.
  TYPED_TEST(SpatialBoundary, MatrixExpressionsPreserveAliasedCoordinates)
  {
    using Scalar = TypeParam;
    for (std::uint8_t rows = 0; rows <= 3; ++rows)
    {
      for (std::uint8_t columns = 0; columns <= 3; ++columns)
      {
        for (const bool array : {false, true})
        {
          SCOPED_TRACE(::testing::Message()
            << unsigned(rows) << 'x' << unsigned(columns) << " array=" << array);
          SpatialMatrix<Scalar> values(rows, columns);
          std::array<Scalar, 9> original{};
          for (size_t row = 0; row < rows; ++row)
          {
            for (size_t column = 0; column < columns; ++column)
            {
              Scalar entry = Scalar(1 + 3 * row + column);
              if constexpr (std::is_same_v<Scalar, Complex>)
                entry += Complex(0, Real(2 + row + column));
              original[3 * row + column] = entry;
              values(row, column) = entry;
            }
          }
          if (array)
            EXPECT_EQ(
              &(values = values.getData().topLeftCorner(rows, columns).array()), &values);
          else
            EXPECT_EQ(&(values = values.getData().topLeftCorner(rows, columns)), &values);
          ASSERT_EQ(values.rows(), rows);
          ASSERT_EQ(values.cols(), columns);
          for (size_t row = 0; row < rows; ++row)
          {
            for (size_t column = 0; column < columns; ++column)
              EXPECT_EQ(values(row, column), original[3 * row + column]);
          }
          if (array)
            values = values.getData().topLeftCorner(rows, columns).array() + Scalar(2);
          else
            values = Scalar(2) * values.getData().topLeftCorner(rows, columns);
          std::array<Scalar, 9> transformed{};
          for (size_t row = 0; row < rows; ++row)
          {
            for (size_t column = 0; column < columns; ++column)
            {
              const Scalar entry = original[3 * row + column];
              transformed[3 * row + column] =
                array ? entry + Scalar(2) : Scalar(2) * entry;
              EXPECT_EQ(values(row, column), transformed[3 * row + column]);
            }
          }
          if (array)
            values = values.getData().topLeftCorner(rows, columns).array().transpose();
          else
            values = values.getData().topLeftCorner(rows, columns).transpose();
          ASSERT_EQ(values.rows(), columns);
          ASSERT_EQ(values.cols(), rows);
          for (size_t row = 0; row < columns; ++row)
          {
            for (size_t column = 0; column < rows; ++column)
              EXPECT_EQ(values(row, column), transformed[3 * column + row]);
          }
          const std::uint8_t shortRows = columns / 2, shortColumns = rows / 2;
          if (array)
            values = values.getData().topLeftCorner(shortRows, shortColumns).array();
          else
            values = values.getData().topLeftCorner(shortRows, shortColumns);
          ASSERT_EQ(values.rows(), shortRows);
          ASSERT_EQ(values.cols(), shortColumns);
          values.resize(3, 3);
          for (size_t row = 0; row < 3; ++row)
          {
            for (size_t column = 0; column < 3; ++column)
            {
              const Scalar expected = row < shortRows && column < shortColumns
                ? transformed[3 * column + row]
                : Scalar(0);
              EXPECT_EQ(values(row, column), expected);
            }
          }
        }
      }
    }
  }

  /// @brief Matrix and array expressions are evaluated before overwriting aliased coordinates.
  TYPED_TEST(SpatialBoundary, EigenExpressionAssignmentPreservesAliasedCoordinates)
  {
    using Scalar = TypeParam;
    for (std::uint8_t size = 0; size <= 3; ++size)
    {
      for (const bool array : {false, true})
      {
        SCOPED_TRACE(::testing::Message() << unsigned(size) << " array=" << array);
        SpatialVector<Scalar> values(size);
        std::array<Scalar, 3> original{};
        for (size_t i = 0; i < size; ++i)
        {
          original[i] = Scalar(i + 1);
          if constexpr (std::is_same_v<Scalar, Complex>)
            original[i] += Complex(0, Real(i + 2));
          values[i] = original[i];
        }
        if (array)
          EXPECT_EQ(&(values = values.getData().head(size).array()), &values);
        else
          EXPECT_EQ(&(values = values.getData().head(size)), &values);
        for (size_t i = 0; i < size; ++i)
          EXPECT_EQ(values[i], original[i]);
        if (array)
          values = values.getData().head(size).array() + Scalar(2);
        else
          values = Scalar(2) * values.getData().head(size);
        std::array<Scalar, 3> transformed{};
        for (size_t i = 0; i < size; ++i)
        {
          transformed[i] = array ? original[i] + Scalar(2) : Scalar(2) * original[i];
          EXPECT_EQ(values[i], transformed[i]);
        }
        if (array)
          values = values.getData().head(size).array().reverse();
        else
          values = values.getData().head(size).reverse();
        for (size_t i = 0; i < size; ++i)
          EXPECT_EQ(values[i], transformed[size - 1 - i]);
        const std::uint8_t shortened = size / 2;
        if (array)
          values = values.getData().head(size).head(shortened).array();
        else
          values = values.getData().head(size).head(shortened);
        ASSERT_EQ(values.size(), shortened);
        for (size_t i = 0; i < shortened; ++i)
          EXPECT_EQ(values[i], transformed[size - 1 - i]);
        values.resize(3);
        for (size_t i = shortened; i < 3; ++i)
          EXPECT_EQ(values[i], Scalar(0));
      }
    }
  }

  /// @brief Empty initializer lists and in-place addition retain zeroed inactive storage.
  TYPED_TEST(SpatialBoundary, EmptyInitializerAndAddition)
  {
    using Scalar = TypeParam;
    SpatialVector<Scalar> values(std::initializer_list<Scalar>{});
    const SpatialVector<Scalar> other(std::initializer_list<Scalar>{});
    EXPECT_EQ(&(values += other), &values);
    EXPECT_EQ(values.size(), 0u);
    values.resize(3);
    for (size_t i = 0; i < 3; ++i)
      EXPECT_EQ(values[i], Scalar(0));
  }

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
