/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <forward_list>
#include <memory>
#include <vector>
#include <set>

#include "Rodin/Geometry/IndexGenerator.h"

using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  // ============================================================================
  // EmptyIndexGenerator Tests
  // ============================================================================

  /// @brief Verifies default construction for geometry empty index generator by checking true predicates.
  TEST(Geometry_EmptyIndexGenerator, DefaultConstruction)
  {
    EmptyIndexGenerator gen;
    EXPECT_TRUE(gen.end());
  }

  /// @brief Verifies copy construction for geometry empty index generator by checking true predicates, copy semantics.
  TEST(Geometry_EmptyIndexGenerator, CopyConstruction)
  {
    EmptyIndexGenerator gen1;
    EmptyIndexGenerator gen2(gen1);
    
    EXPECT_TRUE(gen2.end());
  }

  /// @brief Verifies move construction for geometry empty index generator by checking true predicates, move semantics.
  TEST(Geometry_EmptyIndexGenerator, MoveConstruction)
  {
    EmptyIndexGenerator gen1;
    EmptyIndexGenerator gen2(std::move(gen1));
    
    EXPECT_TRUE(gen2.end());
  }

  /// @brief Verifies always at end for geometry empty index generator by checking true predicates.
  TEST(Geometry_EmptyIndexGenerator, AlwaysAtEnd)
  {
    EmptyIndexGenerator gen;
    EXPECT_TRUE(gen.end());
  }

  /// @brief Verifies copy for geometry empty index generator by checking true predicates, copy semantics.
  TEST(Geometry_EmptyIndexGenerator, Copy)
  {
    EmptyIndexGenerator gen;
    std::unique_ptr<EmptyIndexGenerator> copy(gen.copy());
    
    EXPECT_TRUE(copy->end());
  }

  /// @brief Verifies move for geometry empty index generator by checking true predicates, move semantics.
  TEST(Geometry_EmptyIndexGenerator, Move)
  {
    EmptyIndexGenerator gen;
    std::unique_ptr<EmptyIndexGenerator> moved(gen.move());
    
    EXPECT_TRUE(moved->end());
  }

  // ============================================================================
  // BoundedIndexGenerator Tests
  // ============================================================================

  /// @brief Verifies construction for geometry bounded index generator by checking false predicates.
  TEST(Geometry_BoundedIndexGenerator, Construction)
  {
    BoundedIndexGenerator gen(0, 5);
    EXPECT_FALSE(gen.end());
  }

  /// @brief Verifies empty range for geometry bounded index generator by checking true predicates.
  TEST(Geometry_BoundedIndexGenerator, EmptyRange)
  {
    BoundedIndexGenerator gen(5, 5);
    EXPECT_TRUE(gen.end());
  }

  /// @brief Copy construction preserves the first available bounded index.
  TEST(Geometry_BoundedIndexGenerator, CopyConstruction)
  {
    BoundedIndexGenerator gen1(0, 3);
    BoundedIndexGenerator gen2(gen1);

    EXPECT_FALSE(gen2.end());
    EXPECT_EQ(*gen2, 0);
  }

  /// @brief Verifies move construction for geometry bounded index generator by checking false predicates, move semantics.
  TEST(Geometry_BoundedIndexGenerator, MoveConstruction)
  {
    BoundedIndexGenerator gen1(0, 3);
    BoundedIndexGenerator gen2(std::move(gen1));
    
    EXPECT_FALSE(gen2.end());
  }

  /// @brief Verifies dereference for geometry bounded index generator by checking exact expected values.
  TEST(Geometry_BoundedIndexGenerator, Dereference)
  {
    BoundedIndexGenerator gen(10, 15);
    EXPECT_EQ(*gen, 10);
  }

  /// @brief Verifies increment for geometry bounded index generator by checking exact expected values, true predicates.
  TEST(Geometry_BoundedIndexGenerator, Increment)
  {
    BoundedIndexGenerator gen(0, 3);
    
    EXPECT_EQ(*gen, 0);
    ++gen;
    EXPECT_EQ(*gen, 1);
    ++gen;
    EXPECT_EQ(*gen, 2);
    ++gen;
    EXPECT_TRUE(gen.end());
  }

  /// @brief Verifies full iteration for geometry bounded index generator by checking exact expected values.
  TEST(Geometry_BoundedIndexGenerator, FullIteration)
  {
    BoundedIndexGenerator gen(5, 10);
    
    std::vector<Index> indices;
    while (!gen.end())
    {
      indices.push_back(*gen);
      ++gen;
    }
    
    EXPECT_EQ(indices.size(), 5);
    EXPECT_EQ(indices[0], 5);
    EXPECT_EQ(indices[1], 6);
    EXPECT_EQ(indices[2], 7);
    EXPECT_EQ(indices[3], 8);
    EXPECT_EQ(indices[4], 9);
  }

  /// @brief Polymorphic copying preserves the first available bounded index.
  TEST(Geometry_BoundedIndexGenerator, Copy)
  {
    BoundedIndexGenerator gen(0, 3);
    std::unique_ptr<BoundedIndexGenerator> copy(gen.copy());

    EXPECT_FALSE(copy->end());
    EXPECT_EQ(**copy, 0);
  }

  /// @brief Verifies move for geometry bounded index generator by checking exact expected values, false predicates, move semantics.
  TEST(Geometry_BoundedIndexGenerator, Move)
  {
    BoundedIndexGenerator gen(0, 3);
    std::unique_ptr<BoundedIndexGenerator> moved(gen.move());
    
    EXPECT_FALSE(moved->end());
    EXPECT_EQ(**moved, 0);
  }

  /// @brief Verifies single element for geometry bounded index generator by checking exact expected values, true predicates, false predicates.
  TEST(Geometry_BoundedIndexGenerator, SingleElement)
  {
    BoundedIndexGenerator gen(42, 43);
    
    EXPECT_FALSE(gen.end());
    EXPECT_EQ(*gen, 42);
    ++gen;
    EXPECT_TRUE(gen.end());
  }

  /// @brief Verifies large range for geometry bounded index generator by checking exact expected values.
  TEST(Geometry_BoundedIndexGenerator, LargeRange)
  {
    BoundedIndexGenerator gen(0, 1000);
    
    int count = 0;
    while (!gen.end())
    {
      ++gen;
      ++count;
    }
    
    EXPECT_EQ(count, 1000);
  }

  // ============================================================================
  // VectorIndexGenerator Tests
  // ============================================================================

  /// @brief Verifies empty vector for geometry vector index generator by checking true predicates, move semantics.
  TEST(Geometry_VectorIndexGenerator, EmptyVector)
  {
    std::vector<Index> indices;
    VectorIndexGenerator gen(std::move(indices));
    
    EXPECT_TRUE(gen.end());
  }

  /// @brief Verifies single element for geometry vector index generator by checking exact expected values, true predicates, false predicates.
  TEST(Geometry_VectorIndexGenerator, SingleElement)
  {
    std::vector<Index> indices{42};
    VectorIndexGenerator gen(std::move(indices));
    
    EXPECT_FALSE(gen.end());
    EXPECT_EQ(*gen, 42);
    ++gen;
    EXPECT_TRUE(gen.end());
  }

  /// @brief Verifies multiple elements for geometry vector index generator by checking exact expected values, true predicates, move semantics.
  TEST(Geometry_VectorIndexGenerator, MultipleElements)
  {
    std::vector<Index> indices{10, 20, 30, 40};
    VectorIndexGenerator gen(std::move(indices));
    
    EXPECT_EQ(*gen, 10);
    ++gen;
    EXPECT_EQ(*gen, 20);
    ++gen;
    EXPECT_EQ(*gen, 30);
    ++gen;
    EXPECT_EQ(*gen, 40);
    ++gen;
    EXPECT_TRUE(gen.end());
  }

  /// @brief Verifies full iteration for geometry vector index generator by checking exact expected values, move semantics.
  TEST(Geometry_VectorIndexGenerator, FullIteration)
  {
    std::vector<Index> original{1, 3, 5, 7, 9};
    VectorIndexGenerator gen(std::move(original));
    
    std::vector<Index> result;
    while (!gen.end())
    {
      result.push_back(*gen);
      ++gen;
    }
    
    EXPECT_EQ(result.size(), 5);
    EXPECT_EQ(result[0], 1);
    EXPECT_EQ(result[1], 3);
    EXPECT_EQ(result[2], 5);
    EXPECT_EQ(result[3], 7);
    EXPECT_EQ(result[4], 9);
  }

  /// @brief Verifies copy for geometry vector index generator by checking exact expected values, copy semantics, move semantics.
  TEST(Geometry_VectorIndexGenerator, Copy)
  {
    std::vector<Index> indices{1, 2, 3};
    VectorIndexGenerator gen(std::move(indices));
    
    std::unique_ptr<VectorIndexGenerator> copy(gen.copy());
    
    // Both should iterate the same way
    EXPECT_EQ(**copy, 1);
    ++(*copy);
    EXPECT_EQ(**copy, 2);
  }

  /// @brief Verifies move for geometry vector index generator by checking exact expected values, false predicates, move semantics.
  TEST(Geometry_VectorIndexGenerator, Move)
  {
    std::vector<Index> indices{5, 10, 15};
    VectorIndexGenerator gen(std::move(indices));
    
    std::unique_ptr<VectorIndexGenerator> moved(gen.move());
    
    EXPECT_FALSE(moved->end());
    EXPECT_EQ(**moved, 5);
  }

  // ============================================================================
  // SetIndexGenerator Tests
  // ============================================================================

  /// @brief Verifies empty set for geometry set index generator by checking true predicates, move semantics.
  TEST(Geometry_SetIndexGenerator, EmptySet)
  {
    IndexSet indices;
    SetIndexGenerator gen(std::move(indices));
    
    EXPECT_TRUE(gen.end());
  }

  /// @brief Verifies single element for geometry set index generator by checking exact expected values, true predicates, false predicates.
  TEST(Geometry_SetIndexGenerator, SingleElement)
  {
    IndexSet indices{42};
    SetIndexGenerator gen(std::move(indices));
    
    EXPECT_FALSE(gen.end());
    EXPECT_EQ(*gen, 42);
    ++gen;
    EXPECT_TRUE(gen.end());
  }

  /// @brief Verifies multiple elements sorted for geometry set index generator by checking exact expected values, move semantics.
  TEST(Geometry_SetIndexGenerator, MultipleElementsSorted)
  {
    IndexSet indices{30, 10, 20, 40};  // Set will sort these
    SetIndexGenerator gen(std::move(indices));
    
    // Sets are sorted, so we should get them in order
    std::vector<Index> result;
    while (!gen.end())
    {
      result.push_back(*gen);
      ++gen;
    }
    
    EXPECT_EQ(result.size(), 4);
    EXPECT_EQ(result[0], 10);
    EXPECT_EQ(result[1], 20);
    EXPECT_EQ(result[2], 30);
    EXPECT_EQ(result[3], 40);
  }

  /// @brief Verifies copy for geometry set index generator by checking exact expected values, false predicates, copy semantics.
  TEST(Geometry_SetIndexGenerator, Copy)
  {
    IndexSet indices{1, 2, 3};
    SetIndexGenerator gen(indices);
    
    std::unique_ptr<SetIndexGenerator> copy(gen.copy());
    
    EXPECT_FALSE(copy->end());
    EXPECT_EQ(**copy, 1);
  }

  /// @brief Verifies move for geometry set index generator by checking exact expected values, false predicates, move semantics.
  TEST(Geometry_SetIndexGenerator, Move)
  {
    IndexSet indices{5, 10, 15};
    SetIndexGenerator gen(std::move(indices));
    
    std::unique_ptr<SetIndexGenerator> moved(gen.move());
    
    EXPECT_FALSE(moved->end());
    EXPECT_EQ(**moved, 5);
  }

  /// @brief Verifies construct with const ref for geometry set index generator by checking exact expected values, false predicates.
  TEST(Geometry_SetIndexGenerator, ConstructWithConstRef)
  {
    IndexSet indices{7, 14, 21};
    SetIndexGenerator gen(indices);  // const ref constructor
    
    EXPECT_FALSE(gen.end());
    EXPECT_EQ(*gen, 7);
  }
}

namespace Rodin::Tests::Unit
{
  namespace
  {
    /// @brief Checks a finite suffix without dereferencing an exhausted generator.
    void expectRemainingIndices(
      IndexGeneratorBase& generator, const std::vector<Index>& expected)
    {
      for (const auto index : expected)
      {
        ASSERT_FALSE(generator.end());
        EXPECT_EQ(*generator, index);
        ++generator;
      }
      EXPECT_TRUE(generator.end());
    }

    /// @brief Owned clones preserve the cursor after every previous owner is destroyed.
    template <class Generator, class MakeGenerator>
    void checkOwnedIndexClones(
      MakeGenerator makeGenerator, const std::vector<Index>& sequence)
    {
      for (size_t offset = 0; offset <= sequence.size(); ++offset)
      {
        SCOPED_TRACE(offset);
        auto source = std::make_unique<Generator>(makeGenerator());
        for (size_t step = 0; step < offset; ++step)
          ++*source;
        const IndexGeneratorBase& borrowed = *source;
        std::unique_ptr<IndexGeneratorBase> copied(borrowed.copy());
        std::unique_ptr<IndexGeneratorBase> moved(source->move());
        source.reset();
        const std::vector<Index> suffix(sequence.begin() + offset, sequence.end());
        expectRemainingIndices(*moved, suffix);
        moved.reset();
        expectRemainingIndices(*copied, suffix);
        std::unique_ptr<IndexGeneratorBase> exhaustedCopy(copied->copy());
        std::unique_ptr<IndexGeneratorBase> exhaustedMove(copied->move());
        copied.reset();
        EXPECT_TRUE(exhaustedCopy->end());
        EXPECT_TRUE(exhaustedMove->end());
      }
    }
  }

  /// @brief Bounded clones retain their current position, including exhaustion.
  TEST(Geometry_IndexGeneratorLifecycle, BoundedCursorSurvivesCloning)
  {
    checkOwnedIndexClones<BoundedIndexGenerator>(
      [] { return BoundedIndexGenerator(7, 11); }, {7, 8, 9, 10});
    checkOwnedIndexClones<BoundedIndexGenerator>(
      [] { return BoundedIndexGenerator(7, 7); }, {});
  }

  /// @brief Vector-backed clones own their remaining sequence after source destruction.
  TEST(Geometry_IndexGeneratorLifecycle, VectorCursorOwnsCopiedStorage)
  {
    checkOwnedIndexClones<VectorIndexGenerator>(
      [] { return VectorIndexGenerator(std::vector<Index>{17, 4, 23, 4}); },
      {17, 4, 23, 4});
    checkOwnedIndexClones<VectorIndexGenerator>(
      [] { return VectorIndexGenerator(std::vector<Index>{}); }, {});
  }

  /// @brief Set-backed clones preserve sorted distinct values and the current cursor.
  TEST(Geometry_IndexGeneratorLifecycle, SetCursorOwnsCopiedStorage)
  {
    checkOwnedIndexClones<SetIndexGenerator>(
      [] { return SetIndexGenerator(IndexSet{23, 4, 17, 4}); }, {4, 17, 23});
    checkOwnedIndexClones<SetIndexGenerator>(
      [] { return SetIndexGenerator(IndexSet{}); }, {});
  }

  /// @brief Empty generators remain exhausted through polymorphic copies and moves.
  TEST(Geometry_IndexGeneratorLifecycle, EmptyClonesRemainExhausted)
  {
    checkOwnedIndexClones<EmptyIndexGenerator>([] { return EmptyIndexGenerator(); }, {});
  }

  /// @brief Iterator-backed clones borrow the range while preserving independent cursors.
  TEST(Geometry_IndexGeneratorLifecycle, IteratorClonesPreserveBorrowedSuffix)
  {
    const std::vector<Index> indices{31, 5, 19, 5};
    using Generator = IteratorIndexGenerator<std::vector<Index>::const_iterator>;
    checkOwnedIndexClones<Generator>(
      [&] { return Generator(indices.begin(), indices.end()); }, indices);
    checkOwnedIndexClones<Generator>(
      [&] { return Generator(indices.end(), indices.end()); }, {});
    using PointerGenerator = IteratorIndexGenerator<const Index*>;
    checkOwnedIndexClones<PointerGenerator>(
      [&] { return PointerGenerator(indices.data(), indices.data() + indices.size()); },
      indices);
    const std::forward_list<Index> forwardIndices{31, 5, 19, 5};
    using ForwardGenerator =
      IteratorIndexGenerator<std::forward_list<Index>::const_iterator>;
    checkOwnedIndexClones<ForwardGenerator>(
      [&] { return ForwardGenerator(forwardIndices.begin(), forwardIndices.end()); },
      indices);
  }
}
