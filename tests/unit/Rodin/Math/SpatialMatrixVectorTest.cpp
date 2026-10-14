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
