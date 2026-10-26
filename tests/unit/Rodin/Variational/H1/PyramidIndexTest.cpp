/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Checks the degree-only pyramid mode enumeration exactly.
 */
#include <gtest/gtest.h>
#include <limits>
#include "Rodin/Variational/H1/H1Element.h"
#include "Rodin/Utility/ForConstexpr.h"

namespace Rodin::Tests::Unit
{
  /**
   * @brief Certifies decoding independently by layer, row and column enumeration.
   *
   * Degrees zero through ten cover the constant mode, the benchmark
   * degrees and higher-order metadata. No floating-point tolerance is used.
   * The invalid-index behavior is a compatibility check, not an assertion
   * that an out-of-range mode is meaningful for basis evaluation.
   */
  TEST(H1PyramidIndex, EnumerationAndInverseAcrossDegrees)
  {
    Utility::ForIndex<11>([](auto degree) {
      using Indexing = Variational::PyramidIndex<degree.value>;
      constexpr size_t K = degree.value;
      SCOPED_TRACE(K);
      size_t index = 0;
      for (size_t k = 0; k <= K; ++k)
      {
        EXPECT_EQ(Indexing::getLayerOffset(k), index);
        for (size_t j = 0; j <= K - k; ++j)
        {
          for (size_t i = 0; i <= K - k; ++i, ++index)
          {
            size_t decodedI = Indexing::Count, decodedJ = Indexing::Count;
            size_t decodedK = Indexing::Count;
            Indexing::decode(index, decodedI, decodedJ, decodedK);
            EXPECT_EQ(decodedI, i);
            EXPECT_EQ(decodedJ, j);
            EXPECT_EQ(decodedK, k);
            EXPECT_EQ(Indexing::getIndex(decodedI, decodedJ, decodedK), index);
          }
        }
      }
      EXPECT_EQ(index, Indexing::Count);
      for (const size_t invalid : {Indexing::Count, std::numeric_limits<size_t>::max()})
      {
        size_t i = 1, j = 1, k = 1;
        Indexing::decode(invalid, i, j, k);
        EXPECT_EQ(i, 0u);
        EXPECT_EQ(j, 0u);
        EXPECT_EQ(k, 0u);
      }
    });
  }

  /// @brief Keeps the public index conversion usable in constant expressions.
  TEST(H1PyramidIndex, ConstantExpressionCompatibility)
  {
    constexpr bool correct = [] {
      size_t i = 0, j = 0, k = 0;
      Variational::PyramidIndex<6>::decode(
        Variational::PyramidIndex<6>::getIndex(2, 3, 1), i, j, k);
      return i == 2 && j == 3 && k == 1;
    }();
    static_assert(correct);
    EXPECT_TRUE(correct);
  }
}
