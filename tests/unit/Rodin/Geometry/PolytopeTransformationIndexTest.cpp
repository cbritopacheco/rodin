/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <array>
#include <atomic>
#include <memory>
#include <thread>
#include <vector>

#include <Rodin/Geometry/IdentityTransformation.h>
#include <Rodin/Geometry/PolytopeTransformationIndex.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  /// @brief Verifies that concurrent lookup constructs each transformation once.
  TEST(Geometry_PolytopeTransformationIndex, ConstructsEachSlotOnceConcurrently)
  {
    constexpr size_t count = 128;
    constexpr size_t threadCount = 8;
    constexpr size_t repetitions = 32;

    PolytopeTransformationIndex index;
    index.initialize(2);

    std::array<std::atomic<size_t>, count> factoryCalls{};
    std::vector<std::thread> threads;
    threads.reserve(threadCount);
    for (size_t thread = 0; thread < threadCount; ++thread)
    {
      threads.emplace_back([&, thread]() {
        for (size_t repetition = 0; repetition < repetitions; ++repetition)
          for (size_t offset = 0; offset < count; ++offset)
          {
            const Index i = static_cast<Index>((offset + thread) % count);
            const auto& transformation = index.get({2, i}, count, [&](size_t, Index idx) {
              ++factoryCalls[static_cast<size_t>(idx)];
              return std::make_unique<IdentityTransformation>(2);
            });
            EXPECT_EQ(transformation.getReferenceDimension(), 2u);
            EXPECT_EQ(transformation.getPhysicalDimension(), 2u);
          }
      });
    }

    for (auto& thread : threads)
      thread.join();

    for (const auto& calls : factoryCalls)
      EXPECT_EQ(calls.load(), 1u);
  }

  /// @brief Verifies that clearing the index invalidates cached transformations.
  TEST(Geometry_PolytopeTransformationIndex, ClearRebuildsTransformation)
  {
    PolytopeTransformationIndex index;
    index.initialize(2);

    size_t factoryCalls = 0;
    const auto factory = [&](size_t, Index) {
      ++factoryCalls;
      return std::make_unique<IdentityTransformation>(2);
    };

    index.get({2, 0}, 1, factory);
    index.get({2, 0}, 1, factory);
    EXPECT_EQ(factoryCalls, 1u);

    index.clear();
    index.get({2, 0}, 1, factory);
    EXPECT_EQ(factoryCalls, 2u);
  }
}
